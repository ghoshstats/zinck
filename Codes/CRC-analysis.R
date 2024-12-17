library(ggplot2)
library(rstan)
library(knockoff)
library(randomForest)
library(phyloseq)
library(dplyr)
library(zinck)
library(reshape2)
library(topicmodels)
library(VennDiagram)
library(grid)
library(RColorBrewer)

zinck_code <- "data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num words
  int<lower=0> D; // num docs
  int<lower=0> n[D, V]; // word counts for each doc

  // hyperparameters
  vector<lower=0, upper=1>[V] delta;
}

parameters {
  simplex[K] theta[D]; // topic mixtures
  vector<lower=0, upper=1>[V] zeta[K]; // zero-inflated betas
  vector<lower=0>[V] gamma1[K];
  vector<lower=0>[V] gamma2[K];
  vector<lower=0>[K] alpha;
}

transformed parameters {
  vector[V] beta[K];

  // Efficiently compute beta using vectorized operations
  for (k in 1:K) {
    vector[V] cum_log1m;
    cum_log1m[1:(V - 1)] = cumulative_sum(log1m(zeta[k, 1:(V - 1)]));
    cum_log1m[V] = 0;
    beta[k] = zeta[k] .* exp(cum_log1m);
    beta[k] = beta[k] / sum(beta[k]);
  }
}


model {
  for (k in 1:K) {
    alpha[k] ~ gamma(100,100);  // Change these hyperparameters as needed
  }
  for (d in 1:D) {
    theta[d] ~ dirichlet(alpha);
  }
  for (k in 1:K) {
    for (m in 1:V) {
        gamma1[k,m] ~ gamma(1,1);
        gamma2[k,m] ~ gamma(1,1);
    }
  }

  // Zero-inflated beta likelihood and data likelihood
  for (k in 1:K) {
    for (m in 1:V) {
      real lp_non_zero = bernoulli_lpmf(0 | delta[m]) + beta_lpdf(zeta[k, m] | gamma1[k, m], gamma2[k, m]);
      real lp_zero = bernoulli_lpmf(1 | delta[m]);
      target += log_sum_exp(lp_non_zero, lp_zero);
    }
  }

  // Compute the eta values and data likelihood more efficiently
  for (d in 1:D) {
    vector[V] eta = theta[d, 1] * beta[1];
    for (k in 2:K) {
      eta += theta[d, k] * beta[k];
    }
    eta = eta / sum(eta);
    n[d] ~ multinomial(eta);
  }
}
"

################## CRC data species level ##########################
####################################################################

load("meta.RData")
load("count.RData")

############ Filtering the count matrix in terms of abundances ############
norm_count <- count/rowSums(count)
col_means <- colMeans(norm_count>0)
indices <- which(col_means > 0.2)
sorted_indices <- indices[order(col_means[indices],decreasing = TRUE)]
dcount <- count[,sorted_indices]

X <- dcount
Y <- ifelse(meta$Group=="CRC",1,0)

###################### Method 1: Zinck ##########################################
dlt <- rep(0,ncol(X)) ## Initializing deltas with the individual column sparsities

for(t in (1:ncol(X)))
{
  dlt[t] <- 1-mean(X[,t]>0)
  if(dlt[t]==0)
  {
    dlt[t] = dlt[t]+0.01
  }
  if (dlt[t]==1)
  {
    dlt[t] = dlt[t]-0.01
  }
  
}

zinck_stan_data <- list(
  K = 20,
  V = ncol(X),
  D = nrow(X),
  n = X,
  delta = dlt
)

stan.model = stan_model(model_code = zinck_code)

set.seed(1)
fitCRC <- vb(stan.model, data=zinck_stan_data, algorithm="meanfield", iter=10000) ## Fitting the zinck model

theta <- fitCRC@sim[["est"]][["theta"]]
beta <- fitCRC@sim[["est"]][["beta"]]
X_tilde_CRC <- zinck::generateKnockoff(X,theta,beta,seed=1) ## Generating the knockoff copy

filter_zinck <- zinck.filter(as.matrix(X),as.matrix(X_tilde),as.factor(Y),
                model="Random Forest",offset = 1,seed=1,mtry=28,
                rftuning=TRUE,metric = "Gini")

selected_species <- filter_zinck$selected


## Importance scores ##
W <- filter_zinck$W

## Threshold ##
T <- filter_zinck$T

print(which(W>=T))
## Selected Taxa :: 36,76,150,158,164,204,236,240,255,259,280,311,331,374

################# Method 2: Model-X Knockoff Filter ############################

W_kf <- log_normalize(X)
set.seed(1)
W_tilde_kf <- create.second_order(W_kf) # Generating the second order Gaussian knockoff copy
set.seed(45)
W <- stat.random_forest(W_kf,W_tilde_kf,as.factor(Y))
T <- knockoff.threshold(W,fdr=0.1,offset=0)
print(which(W>=T))

## Selected Taxa :: 36, 76, 158, 331
names_KF <- colnames(X[,which(W>=T)])

################# Method 3: Vanilla LDA model ##################################
df.LDA <- as(as.matrix(X),"dgCMatrix")
vanilla.LDA.CRC <- LDA(df.LDA,k=16,method="VEM")
theta.LDA.CRC <- vanilla.LDA.CRC@gamma
beta.LDA.CRC <- vanilla.LDA.CRC@beta
beta.LDA.CRC <- t(apply(beta.LDA.CRC, 1, function(row) row/sum(row)))
X_tilde_LDA.CRC <- zinck::generateKnockoff(X,theta.LDA.CRC,beta.LDA.CRC,seed=1) ## Generating the vanilla LDA knockoff copy
X_aug <- cbind(X,X_tilde_LDA.CRC) ## Creating the augmented matrix

set.seed(6)
W <- stat.random_forest(X,X_tilde_LDA.CRC,Y)
T <- knockoff.threshold(W,fdr=0.1,offset=1)
print(which(W>=T))

## Selected Taxa :: 5,  36, 150, 158, 164,204, 236, 255, 311, 331, 374
names_LDA <- colnames(X[,which(W>=T)])


################## Method 4: DeepLINK ##########################################

```{python}
import numpy as np
import DeepLINK as dl
import pandas as pd
from PCp1_numFactors import PCp1 as PCp1
import keras
from keras.layers import Dense, Dropout
from keras.models import Sequential
from pairwise_connected_layer import PairwiseConnected
from itertools import combinations
from keras.callbacks import EarlyStopping
import tensorflow as tf
import random
import numpy as np

aut_epoch = 100 # number of autoencoder training epochs
aut_loss = 'mean_squared_error' # loss function used in autoencoder training
aut_verb = 0 # verbose level of autoencoder
mlp_epoch = 100 # number of mlp training epochs
mlp_loss = 'binary_crossentropy'
#mlp_loss = 'mean_squared_error' # loss function used in mlp training
dnn_loss = 'binary_crossentropy'
dnn_verb = 0
aut_met = 'relu'
dnn_met = 'elu'
mlp_verb = 0 # verbose level of mlp
l1 = 0.001 # l1 regularization factor in mlp
lr = 0.001 # learning rate for mlp training

X -= np.mean(X, axis=0)

# Prevent division by zero: Replace zero std with 1 before dividing
std_dev = np.std(X, axis=0, ddof=1)
std_dev[std_dev == 0] = 1  # Avoid division by zero
X /= std_dev

# Normalize X1 while ensuring no division by zero in normalization
norms = np.sqrt(np.sum(X ** 2, axis=0))
norms[norms == 0] = 1  # Avoid division by zero
X1 = X / norms
r_hat = PCp1(X1, 15)
random.seed(4)
np.random.seed(4)
tf.random.set_seed(4)
Xnew = dl.knockoff_construct(X, r_hat, 'relu', aut_epoch, aut_loss, aut_verb)

  # compute knockoff statistics

p = Xnew.shape[1] // 2

 # implement DeepPINK
es = EarlyStopping(monitor='val_loss', patience=30, verbose=2)
dp = Sequential()
dp.add(PairwiseConnected(input_shape=(2 * p,)))
dp.add(Dense(p, activation='elu', kernel_regularizer=keras.regularizers.l1(l1=l1)))
dp.add(Dense(1, activation='relu',
                    kernel_regularizer=keras.regularizers.l1_l2(l1=0.001, l2=0.001)))
dp.compile(loss=mlp_loss, optimizer=keras.optimizers.Adam(learning_rate=lr))
dp.fit(Xnew, y, epochs=mlp_epoch, batch_size=32, verbose=mlp_verb, validation_split=0.1, callbacks=[es])

 # calculate knockoff statistics W_j
weights = dp.get_weights()
w = weights[1] @ weights[3]
w = w.reshape(p, )
z = weights[0][:p]
z_tilde = weights[0][p:]
W = (w * z) ** 2 - (w * z_tilde) ** 2
  # feature selection
np.random.seed(1)

selected = dl.knockoff_select(W, 0.1, ko_plus=False)
```
                        
## Selected Taxa :: 6,  19, 36, 106, 126, 161, 164, 234, 240, 241, 255, 294, 349
names_DL <- colnames(X[,which(W>=T)])

################################################################################
################ Plot the Feature importance Statistics ########################

######### (1) Zinck #############
#################################

names_zinck <- c("Ruminococcus torques species (1376)","Bacteroides caccae species (1382)","Clostridium symbiosum species (1475)","Hungatella hathewayi species (0882)","Unknown Ruminococcus species (6664)","Unknown Clostridiales species (6105)",  "Parvimonas micra species (1145)","Unknown Clostridium species (5413)","Gemella morbillorum species (4513)","Unknown Clostridiales species (6009)","Clostridium clostridioforme species (0979)", "Peptostreptococcus stomatis species (4614)","Solobacterium moorei species (0531)","Fusobacterium nucleatum species (0776)")

data.species <- data.frame(             ### Creating the data frame with Feature Importance Statistics
  impscores = sort(W[which(W>=T)], decreasing=FALSE) , 
  name = factor(names_zinck, levels = names_zinck),
  y = seq(length(names_zinck)) * 0.9
)

norm_count <- count/rowSums(count)
col_means <- colMeans(norm_count>0)
indices <- which(col_means > 0)
sorted_indices <- indices[order(col_means[indices],decreasing = TRUE)]
Xnorm <-norm_count[,sorted_indices]

# Calculate the column sums for cases and controls
case_sums <- colMeans(Xnorm[Y == 1, which(W>=T)])
control_sums <- colMeans(Xnorm[Y == 0, which(W>=T)])

# Determine colors based on the sum comparison
colors <- ifelse(case_sums > control_sums, "red", "blue")

# Create a vector to indicate which labels should be bold
bold_indices <- c(10,11)
bold_labels <- ifelse(seq_along(names_zinck) %in% bold_indices, "bold", "plain")

# Create a data frame for plotting
data.species <- data.frame(
  impscores = sort(W[which(W >= T)], decreasing = FALSE),
  name = factor(names_zinck, levels = names_zinck),
  y = seq(length(names_zinck)) * 0.9,
  color = colors,
  fontface = bold_labels
)

tiff("Zinck_CRC_feature_scores_not_log_normalized.tiff", units="in", width=10, height=6, res=600)

plt.species <- ggplot(data.species) +
  geom_col(aes(impscores, name, fill = color), width = 0.6) +
  scale_fill_identity() +
  theme_bw() +
  ylab("Species") +
  xlab("Feature Statistic") +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, face = data.species$fontface) # Adjust fontface based on the condition
  )
print(plt.species)
dev.off()

######### (2) Model-X Knockoffs #######
#######################################

data.species <- data.frame(
  impscores = sort(W[which(W>=T)], decreasing=FALSE) , 
  name = factor(names, levels = names),
  y = seq(length(names)) * 0.9
)

plt.species <- ggplot(data.species) +
  geom_col(aes(impscores, name), fill = "black", width = 0.6)+theme_bw()+ylab("Species")+xlab("Feature Statistic")

plt.species

######### (3) vanilla LDA Knockoffs #########
#############################################

data.species <- data.frame(
  impscores = sort(W[which(W>=T)], decreasing=FALSE) , 
  name = factor(names, levels = names),
  y = seq(length(names)) * 0.9
)


plt.species <- ggplot(data.species) +
  geom_col(aes(impscores, name), fill = "black", width = 0.6)+theme_bw()+ylab("Species")+xlab("Feature Statistic")

plt.species

######### (4) DeepLINK knockoffs ###########
############################################

data.species <- data.frame(
  impscores = sort(W[which(W>=T)], decreasing=FALSE) , 
  name = factor(names, levels = names),
  y = seq(length(names)) * 0.9
)

plt.species <- ggplot(data.species) +
  geom_col(aes(impscores, name), fill = "black", width = 0.6)+theme_bw()+ylab("Species")+xlab("Feature Statistic")

plt.species

############## Draw the Venn Diagram ########################
#############################################################

####### Set of selected taxa #######
zinCK <- c(36,76,150,158,164,204,236,240,255,259,280,311,331,374)
DL <- c(6,  19, 36, 106, 126, 161, 164, 234, 240, 241, 255, 294, 349)
KF <- c(36, 76, 158, 331)
vanilla_LDA <- c(5,  36, 150, 158, 164,204, 236, 255, 311, 331, 374)

venn.plot <- venn.diagram(
  x = list(zinCK = zinCK, KF = KF, DeepLINK = DL, vanilla_LDA = vanilla_LDA),
  category.names = c("Zinck" , "MX-KF" , "DeepLINK", "LDA-KF"),
  filename = 'CRC_venn_diagram_not_lognormalized.png',
  output = TRUE,  # Set to FALSE so that we can draw it with grid.draw
  imagetype="png",
  height = 600,
  width = 600,
  resolution = 400,
  lwd = 1,
  col=c("red", 'blue', 'green', "orange"),
  fill = c(alpha("red", 0.3), alpha('blue', 0.3), alpha('green', 0.3), alpha('orange', 0.3)),
  cex = 0.7,
  fontfamily = "sans",
  cat.cex = 0,  # Modified for readability
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0, 0),  # Modified for 4 categories
  cat.dist = c(0.055, 0.055, 0.055, 0.055),  # Modified for 4 categories
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', 'green', 'orange')  # Modified for 4 categories
)

