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

################################################################################
###################### IBD data genus level ####################################

load("genera.RData") ## Loading the meta IBD studies
combined_studies <- as.data.frame(t(physeq_genera@otu_table))
study_names <- physeq_genera@sam_data[["dataset_name"]]

## RISK ##
risk_indices <- which(study_names == "RISK")
risk_otu <- combined_studies[risk_indices, ]

## CS-PRISM ##
prism_indices <- which(study_names == "CS-PRISM")
prism_otu <- combined_studies[prism_indices, ]

## HMP2 ##
HMP_indices <- which(study_names == "HMP2")
hmp_otu <- combined_studies[HMP_indices, ]

## Pouchitis ##
pouchitis_indices <- which(study_names == "Pouchitis")
pouchitis_otu <- combined_studies[pouchitis_indices, ]

## MucosalIBD ##
mi_indices <- which(study_names == "MucosalIBD")
mi_otu <- combined_studies[mi_indices, ]

metadata <- physeq_genera@sam_data

# Get the indices of the last occurrence of each unique subject_accession
latest_indices <- sapply(unique(metadata$subject_accession), function(x) {
  max(which(metadata$subject_accession == x))
})

# Subset the metadata to keep only the latest samples
latest_metadata <- metadata[latest_indices, ]

study_names <- c("CS-PRISM","HMP2","MucosalIBD","Pouchitis","RISK")
meta_IBD <- latest_metadata[latest_metadata$dataset_name %in% study_names,]

# Extract row names for each study
risk_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "RISK", ])
prism_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "CS-PRISM", ])
hmp_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "HMP2", ])
mi_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "MucosalIBD", ])
pouchitis_row_names <- row.names(latest_metadata[latest_metadata$dataset_name == "Pouchitis", ])

# Subset the OTU matrices and Y vectors based on these row names
risk_otu_latest <- risk_otu[row.names(risk_otu) %in% risk_row_names, ]
prism_otu_latest <- prism_otu[row.names(prism_otu) %in% prism_row_names, ]
hmp_otu_latest <- hmp_otu[row.names(hmp_otu) %in% hmp_row_names, ]
mi_otu_latest <- mi_otu[row.names(mi_otu) %in% mi_row_names, ]
pouchitis_otu_latest <- pouchitis_otu[row.names(pouchitis_otu) %in% pouchitis_row_names, ]

combined_otu <- rbind(prism_otu_latest,hmp_otu_latest,mi_otu_latest,
                      pouchitis_otu_latest,risk_otu_latest)
combined_otu <- combined_otu[ match(rownames(meta_IBD), rownames(combined_otu)), ] 
## To make sure samples in OTU correspond to the samples in meta

IBD_resp <- meta_IBD$disease
Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0)

X <- combined_otu

########### (1) Zinck ############
dlt <- rep(0,ncol(X)) ## Initializing the deltas with the sparsity of each column

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

zinLDA_stan_data <- list(
  K = 27,
  V = ncol(X),
  D = nrow(X),
  n = X,
  delta = dlt
)

## Fitting the zinck model ##
set.seed(1)
fitIBD <- vb(stan.model, data=zinLDA_stan_data, algorithm="meanfield", importance_resampling=TRUE, iter=10000,tol_rel_obj=0.01,elbo_samples=500)

theta <- fitIBD@sim[["est"]][["theta"]]
beta <- fitIBD@sim[["est"]][["beta"]]
X_tilde <- zinck::generateKnockoff(X,theta,beta,seed=1) ## Generating the kncokoff copy

filter_zinck <- zinck.filter(as.matrix(X),as.matrix(X_tilde),as.factor(Y),
                             model="Random Forest",offset=0,seed=312)

selected_genera <- filter_zinck$selected


## Importance scores ##
W <- filter_zinck$W

## Threshold ##
T <- filter_zinck$T

print(which(W>=T))

# 2, 30,35, 36, 39, 55,58,69,77,79,90, 94 selected

############## Feature Importance Plot #############

names_zinck <- c("Bilophila genus (f. Desulfovibrionaceae)","Unknown genus (f. Ruminococcaceae)","Holdemania genus (f. Erysipelotrichaceae)","Parabacteroides genus (f. Porphyromonadaceae)","Unknown genus (f. Rikenellaceae)", "Bacteroides genus (f. Bacteroidaceae)","Oscillospira genus (f. Ruminococcaceae)","Phascolarctobacterium genus (f. Veillonellaceae)","Coprococcus genus (f. Lachnospiraceae)","Blautia genus (f. Lachnospiraceae)", "Lachnospira genus (f. Lachnospiraceae)","Faecalibacterium genus (f. Ruminococcaceae)")
Xnorm <- X / rowSums(X)

# Calculate the column sums for cases and controls
case_sums <- colMeans(Xnorm[Y == 1, which(W>=T)])
control_sums <- colMeans(Xnorm[Y == 0, which(W>=T)])

colors <- ifelse(case_sums > control_sums, "red", "blue")


# Create a vector to indicate which labels should be bold
bold_indices <- c(6,10)
bold_labels <- ifelse(seq_along(names_zinck) %in% bold_indices, "bold", "plain")

# Create a data frame for plotting
data.species <- data.frame(
  impscores = sort(W[which(W >= T)], decreasing = FALSE),
  name = factor(names_zinck, levels = names_zinck),
  y = seq(length(names_zinck)) * 0.9,
  color = colors,
  fontface = bold_labels
)

tiff("Zinck_IBD_feature_scores_notlognormalized.tiff", units="in", width=10, height=6, res=600)

# Create the plot
plt.species <- ggplot(data.species) +
  geom_col(aes(impscores, name, fill = color), width = 0.6) +
  scale_fill_identity() +
  theme_bw() +
  ylab("Genera") +
  xlab("Feature Statistic") +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18, face = data.species$fontface) # Adjust fontface based on the condition
  )
print(plt.species)
dev.off()

######### (2) MX-KF #########
W <- log_normalize(X)
set.seed(3)
W_tilde.KF <- create.second_order(W)
set.seed(46)
imp <- stat.random_forest(W,W_tilde.KF,as.factor(Y))
T <- knockoff.threshold(imp,fdr=0.1,offset=0)
print(which(imp>=T))

######### (3) LDA-KF ############
df.LDA <- as(as.matrix(X),"dgCMatrix") ## Converting the matrix X into a dgCMatrix
vanilla.LDA.comb <- LDA(df.LDA,k=20,method="VEM") ## Fitting LDA()
theta.LDA.comb <- vanilla.LDA.comb@gamma
beta.LDA.comb <- vanilla.LDA.comb@beta
beta.LDA.comb <- t(apply(beta.LDA.comb, 1, function(row) row/sum(row))) ## Normalizing the posterior beta estimates
X_tilde.LDA.comb <- zinck::generateKnockoff(X,theta.LDA.comb,beta.LDA.comb,seed=1) ## Standard LDA Knockoff matrix for X

set.seed(2)
W <- stat.random_forest(X,X_tilde.LDA.comb,as.factor(Y))
T <- knockoff.threshold(W,fdr=0.1,offset=0)
print(which(W>=T))

########## (4) DeepLINK ##########

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
random.seed(3)
np.random.seed(3)
tf.random.set_seed(3)
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


############################ Plotting the Venn Diagram ################################
#######################################################################################
                      
zinCK <- c(2, 30,35, 36, 39, 55,58,69,77,79,90, 94)
DL <- c(12,17,18,30,36,44,49,52,70,88,90,94,96)
KF <- c(2, 19, 39, 58, 77, 90)
vanilla_LDA <- c(2, 30,35,36, 39,58,69,77, 90,94)
venn.plot <- venn.diagram(
  x = list(zinCK = zinCK, KF = KF, DeepLINK = DL, vanilla_LDA = vanilla_LDA),
  category.names = c("Zinck" , "MX-KF" , "DeepLINK", "LDA-KF"),
  filename = 'IBD_venn_diagram_not_lognormalized.png',
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




