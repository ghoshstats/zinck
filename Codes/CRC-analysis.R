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

##### CRC data species level ######
###################################

load("meta.RData")
load("count.RData")
dcount <- count[,order(decreasing=T,colSums(count,na.rm=T),apply(count,2L,paste,collapse=''))][,1:300] ## ordering the columns w/ decreasing abundance
X <- dcount
Y <- as.factor(meta$Group)
lookup <- c("CTR" = 0, "CRC" = 1)
Y <- lookup[Y]    ## Converting into 0/1 data

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

zinLDA_stan_data <- list(
  K = 19,
  V = ncol(X),
  D = nrow(X),
  n = X,
  delta = dlt
)

stan.model = stan_model(model_code = zinck_code)

set.seed(1)
fitCRC <- vb(stan.model, data=zinLDA_stan_data, algorithm="meanfield", iter=10000) ## Fitting the zinck model

theta <- fitCRC@sim[["est"]][["theta"]]
beta <- fitCRC@sim[["est"]][["beta"]]
X_tilde_CRC <- zinck::generateKnockoff(X,theta,beta,seed=1) ## Generating the knockoff copy
X_aug <- cbind(X,X_tilde_CRC) ## Creating the augmented matrix

######### Tuning the Random Forest model ####################

bestmtry <- tuneRF(X_aug,as.factor(Y),stepFactor = 1.5,improve=1e-5,ntree=1000) ## Getting the best mtry hyperparameter
m <- bestmtry[as.numeric(which.min(bestmtry[,"OOBError"])),1]
df_X <- as.data.frame(X_aug)
colnames(df_X) <- NULL
rownames(df_X) <- NULL
df_X$Y <- Y
model_rf <- randomForest(formula=as.factor(Y)~.,ntree=1000,mtry=m,importance=TRUE,data=as.matrix(df_X)) ## Fitting the tuned Random Forest model
cf <- as.data.frame(importance(model_rf))[,3] ## Extracting the Mean Decrease in Impurities for each variable
W <- abs(cf[1:300])-abs(cf[301:600])
T <- knockoff.threshold(W,fdr = 0.1, offset = 0) ## This is the knockoff threshold
print(which(W>=T))
names_zinck <- colnames(X[,which(W>=T)])
## Selected Taxa :: 8,  33,  47,  54,  81, 124, 130, 136, 146, 193, 215, 245, 255, 263, 264, 268, 283

################# Method 2: Model-X Knockoff Filter ############################

X_tilde_kf <- create.second_order(X) ## Generating the second order Gaussian knockoff copy
set.seed(2)
W <- stat.random_forest(X,X_tilde_kf,as.factor(Y)) ## Fitting the Random Forest model
T <- knockoff.threshold(W,fdr=0.1,offset=0)
print(which(W>=T))
names_KF <- colnames(X[,which(W>=T)])
## Selected Taxa :: 245, 263

################# Method 3: Vanilla LDA model ##################################

df.LDA <- as(as.matrix(X),"dgCMatrix")
vanilla.LDA.CRC <- LDA(df.LDA,k=16,method="VEM")
theta.LDA.CRC <- vanilla.LDA.CRC@gamma
beta.LDA.CRC <- vanilla.LDA.CRC@beta
beta.LDA.CRC <- t(apply(beta.LDA.CRC, 1, function(row) row/sum(row)))
X_tilde_LDA.CRC <- zinck::generateKnockoff(X,theta.LDA.CRC,beta.LDA.CRC,seed=1) ## Generating the vanilla LDA knockoff copy

set.seed(47)
W <- stat.random_forest(X,X_tilde_LDA.CRC,Y)
T <- knockoff.threshold(W,fdr=0.1,offset=0)
print(which(W>=T))
names_LDA <- colnames(X[,which(W>=T)])
## Selected Taxa :: 8, 47, 136, 245, 263 

################## Method 4: DeepLINK ##########################################

## Selected Taxa :: 6, 47, 54, 63, 84, 89, 114, 124, 136, 137, 161, 173, 174, 176, 194, 207, 223, 240, 245, 247
names_DL <- colnames(X[,which(W>=T)])

################################################################################
################ Plot the Feature importance Statistics ########################

######### (1) Zinck #############
#################################

data.species <- data.frame(             ### Creating the data frame with Feature Importance Statistics
  impscores = sort(W[which(W>=T)], decreasing=FALSE) , 
  name = factor(names_zinck, levels = names_zinck),
  y = seq(length(names)) * 0.9
)


plt.species <- ggplot(data.species) +
  geom_col(aes(impscores, name), fill = "black", width = 0.6)+theme_bw()+ylab("Species")+xlab("Feature Statistic")

plt.species

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
zinCK <- c(8,  33,  47,  54,  81, 124, 130, 136, 146, 193, 215, 245, 255, 263, 264, 268, 283)
DL <- c(6, 47, 54, 63, 84, 89, 114, 124, 136, 137, 161, 173, 174, 176, 194, 207, 223, 240, 245, 247)
KF <- c(245, 263)
vanilla_LDA <- c(8, 47, 136, 245, 263 )

venn.plot <- venn.diagram(
  x = list(zinCK = zinCK, KF = KF, DeepLINK = DL, vanilla_LDA = vanilla_LDA),
  category.names = c("Zinck" , "KF" , "DeepLINK", "vanilla_LDA"),
  filename = 'CRC_species_venn_diagram.png',
  output = TRUE,  # Set to FALSE so that we can draw it with grid.draw
  imagetype="png",
  height = 600,
  width = 600,
  resolution = 300,
  lwd = 1,
  col=c("red", 'blue', 'green', "orange"),
  fill = c(alpha("red", 0.3), alpha('blue', 0.3), alpha('green', 0.3), alpha('orange', 0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.6,  # Modified for readability
  cat.default.pos = "outer",
  cat.pos = c(0, 0, 0, 0),  # Modified for 4 categories
  cat.dist = c(0.055, 0.055, 0.055, 0.055),  # Modified for 4 categories
  cat.fontfamily = "sans",
  cat.col = c("red", 'blue', 'green', 'orange')  # Modified for 4 categories
)


