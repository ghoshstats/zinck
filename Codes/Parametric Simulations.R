###################################### Replicating the simulation settings for CKF with varying p #############################################
###############################################################################################################################################

library(knockoff)
library(glmnet)
library(MethylCapSig)
library(GUniFrac)
library(energy)
library(dirmult)
library(zinck)
library(kosel)
library(rstan)
library(topicmodels)
library(dplyr)
library(reshape2)

load("DirMultOutput.RData")

n = 500
p = 100 ## Change to 200, 300, 400
#niter = 100
#################### Stan code for zinck ####################

zinck_code <- "data {
    int<lower=1> K; // num clusters
    int<lower=1> V; // num words
    int<lower=0> D; // num docs
    int<lower=0> n[D, V]; // word counts for each doc

    // hyperparameters
    vector<lower=0, upper=1>[V] delta;
    }

    parameters {
      simplex[K] theta[D]; // cluster mixtures
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

stan.model <- stan_model(model_code = zinck_code)

mvlognormal<-function(n , Mu, Sigma, R, seed){
  ## Generate Lognormal random variables with zeros. 
  ## Mu - Mean of the actual lognormal variables
  ## Sigma - Diagonal of the covariance matrix of the actual variables.
  ## R - Correlation matrix for the log-transformed normal variables. This is to ensure the normal distribution corresponding to the lognormals exists
  if (dim(R)[1] != dim(R)[2]){
    stop("Correlation matrix is not square.");
  }
  p <- length(Mu);
  if (length(Sigma) != p){
    stop("Mean and covariance matrix are not of same dimension")
  }
  Alpha <- matrix(0,p,1);
  Beta <- matrix(0,p,p);
  
  ## Alpha and Beta are the converted mean and covariance-diagonal matrix of the log-transformed normal variable.
  for (i in 1:p){
    if (abs(Mu[i]) >= .Machine$double.eps){
      Alpha[i] = log(Mu[i]) - (1/2)*log(1 + Sigma[i]/Mu[i]^2);
      Beta[i,i] = log(1 + Sigma[i]/(Mu[i]^2));
    }
    if (abs(Mu[i]) < .Machine$double.eps){
      Alpha[i] = -Inf;
      Beta[i,i] = 0;
    }
  }
  Delta = sqrt(Beta)%*%R%*%sqrt(Beta);
  Delta.Eigen = eigen(Delta);
  Delta.EigenVal <- Delta.Eigen$values*(abs(Delta.Eigen$values) >= 1e-12);
  RootDelta = (Delta.Eigen$vectors)%*%diag(sqrt(Delta.EigenVal))%*%t(Delta.Eigen$vectors);
  set.seed(seed)
  X <- matrix(rnorm(p*n,0,1), nrow = n, ncol = p);
  for (i in 1:n){
    X[i,] <- exp(Alpha + RootDelta%*%X[i,]);
  }      
  return(X);
}

########################## DM data generation ############################
power1_DM_cts_list <- c()
power2_DM_cts_list <- c()
power3_DM_cts_list <- c()
power4_DM_cts_list <- c()

power1_DM_bin_list <- c()
power2_DM_bin_list <- c()
power3_DM_bin_list <- c()
power4_DM_bin_list <- c()


fdr1_DM_cts_list <- c()
fdr2_DM_cts_list <- c()
fdr3_DM_cts_list <- c()
fdr4_DM_cts_list <- c()

fdr1_DM_bin_list <- c()
fdr2_DM_bin_list <- c()
fdr3_DM_bin_list <- c()
fdr4_DM_bin_list <- c()

data(throat.tree)
data(throat.otu.tab)
tree = throat.tree 
tree.dist = cophenetic(tree)
otu.ids <- tree$tip.label

for(i in 1:500)
{
set.seed(i)
p.est = dd$pi
names(p.est) <- names(dd$pi)
theta <- dd$theta
gplus <- (1 - theta) / theta
p.est <- p.est[otu.ids]
g.est <- p.est * gplus
comm <- matrix(0, n , length(g.est))  
rownames(comm) <- 1:nrow(comm)
colnames(comm) <- names(g.est)
comm.p <- comm   # comm.p hold the underlying proportions
nSeq <- rnbinom(n , mu = 25000, size = 25)
for (i in 1:n ) {
    comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
    comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
}
colindex=order(p.est, decreasing = T)[1:p]
X = comm[,colindex]
Pi = comm.p[,colindex]
rowPi = rowSums(Pi)
Pi = sweep(Pi, 1, rowPi, "/")
#W = log_normalize(X)
##### Generating continuous responses ######
# rowX = rowSums(X)
# Pi = sweep(X, 1, rowX, "/")
col_abundances = colMeans(Pi)
set.seed(1)
signals = (2 * rbinom(30, 1, 0.5) - 1) * runif(30, 1.5, 3)
kBeta = c(signals / sqrt(col_abundances[1:30]), rep(0, p - 30))
eps=rnorm(n,mean = 0, sd=1)
Y <- Pi^2 %*% (kBeta/2) + Pi %*% kBeta + eps
# X <- matrix(0, nrow = nrow(Pi), ncol = ncol(Pi))
 
# # Loop over each row to generate the new counts based on the multinomial distribution
 
# nSeq = 25000
# set.seed(1)
# for (i in 1:nrow(Pi)) {
#   X[i, ] <- rmultinom(1, size = nSeq, prob = Pi[i, ])
# }
 
colnames(X) <- colnames(Pi)
rownames(X) <- rownames(Pi)
dlt <- rep(0,ncol(X))
for(t in (1:ncol(X)))
{
dlt[t] <- 1-mean(X[,t]>0)
}
    
zinck_stan_data <- list(
K = 16,
V = ncol(X),
D = nrow(X),
n = X,
delta = dlt
)
stan.model <- stan_model(model_code = zinck_code)
    
set.seed(1)
    
fit_zinck <- vb(
stan.model,
data = zinck_stan_data,
algorithm = "meanfield",
importance_resampling = TRUE,
iter = 10000,
tol_rel_obj = 0.01,
elbo_samples = 500
)
beta <- fit_zinck@sim[["est"]][["beta"]]
theta <- fit_zinck@sim[["est"]][["theta"]]
    
X_tilde <- generateKnockoff(X,theta,beta,seed=1)
index <- 1:30 ## Index set of the true non-zero signals
################### Varying the target FDR thresholds #########################
############## FDR = 0.05 ###############
index_est <- zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.05)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr1_DM_cts_list[i] <- FP/(FP+TP)
power1_DM_cts_list[i] <- TP/(TP+FN)
    
############### FDR = 0.1 ################
index_est <- zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.1)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr2_DM_cts_list[i] <- FP/(FP+TP)
power2_DM_cts_list[i] <- TP/(TP+FN)
    
################ FDR = 0.15 ###############
index_est <- zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.15)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr3_DM_cts_list[i] <- FP/(FP+TP)
power3_DM_cts_list[i] <- TP/(TP+FN)
    
################ FDR = 0.2 ################
index_est <- zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.2)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr4_DM_cts_list[i] <- FP/(FP+TP)
power4_DM_cts_list[i] <- TP/(TP+FN)

########## Generating binary response ############
set.seed(1)
signals = (2 * rbinom(30, 1, 0.5) - 1) * runif(30, 3, 10)
kBeta = c(signals / sqrt(col_abundances[1:30]), rep(0, p - 30))
pr = 1/(1+exp(-(Pi^2 %*% (kBeta/2) + Pi %*% kBeta)))
Y_bin = rbinom(n,1,pr)

################### Varying the target FDR thresholds #########################
############## FDR = 0.05 ###############
index_est <- zinck.filter(X,X_tilde,as.factor(Y_bin),model="Random Forest",fdr=0.05)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr1_DM_bin_list[i] <- FP/(FP+TP)
power1_DM_bin_list[i] <- TP/(TP+FN)
    
############### FDR = 0.1 ################
index_est <- zinck.filter(X,X_tilde,as.factor(Y_bin),model="Random Forest",fdr=0.1)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr2_DM_bin_list[i] <- FP/(FP+TP)
power2_DM_bin_list[i] <- TP/(TP+FN)
    
################ FDR = 0.15 ###############
index_est <- zinck.filter(X,X_tilde,as.factor(Y_bin),model="Random Forest",fdr=0.15)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr3_DM_bin_list[i] <- FP/(FP+TP)
power3_DM_bin_list[i] <- TP/(TP+FN)
    
################ FDR = 0.2 ################
index_est <- zinck.filter(X,X_tilde,as.factor(Y_bin),model="Random Forest",fdr=0.2)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr4_DM_bin_list[i] <- FP/(FP+TP)
power4_DM_bin_list[i] <- TP/(TP+FN)
}





########################## LN data generation ############################
power1_LN_cts_list <- c()
power2_LN_cts_list <- c()
power3_LN_cts_list <- c()
power4_LN_cts_list <- c()

power1_LN_bin_list <- c()
power2_LN_bin_list <- c()
power3_LN_bin_list <- c()
power4_LN_bin_list <- c()


fdr1_LN_cts_list <- c()
fdr2_LN_cts_list <- c()
fdr3_LN_cts_list <- c()
fdr4_LN_cts_list <- c()

fdr1_LN_bin_list <- c()
fdr2_LN_bin_list <- c()
fdr3_LN_bin_list <- c()
fdr4_LN_bin_list <- c()

sigma = matrix(0,p,p)
gamma = 0.5
for(i in 1:nrow(sigma)){
    for(j in 1:nrow(sigma)){
        sigma[i,j] = gamma^(abs(i-j))   
    }
}
### Generate mean ###
means = rep(1,p)

for(i in 1:200)
{
W = mvlognormal(n, means, Sigma = diag(sigma), R = sigma,seed=i)
rowW = rowSums(W)
Pi = sweep(W, 1, rowW, "/")
colnames(Pi) <-  c(paste0("Taxa", 1:p))

X <- matrix(0, nrow = nrow(Pi), ncol = ncol(Pi))
    
# Loop over each row to generate the new counts based on the multinomial distribution
nSeq = 25000
set.seed(1)
for (i in 1:nrow(Pi)) {
    X[i, ] <- rmultinom(1, size = nSeq, prob = Pi[i, ])
}
    
colnames(X) <- colnames(Pi)
rownames(X) <- rownames(Pi)
##### Generating continuous responses ######
col_abundances = colMeans(Pi)
set.seed(1)
signals = (2 * rbinom(30, 1, 0.5) - 1) * runif(30, 1.5, 3)
kBeta = c(signals / sqrt(col_abundances[1:30]), rep(0, p - 30))
eps=rnorm(n,mean = 0, sd=1)
Y <- Pi^2 %*% (kBeta/2) + Pi %*% kBeta + eps

dlt <- rep(0,ncol(X))
for(t in (1:ncol(X)))
{
dlt[t] <- 1-mean(X[,t]>0)
}
    
zinck_stan_data <- list(
K = 16,
V = ncol(X),
D = nrow(X),
n = X,
delta = dlt
)
stan.model <- stan_model(model_code = zinck_code)
    
set.seed(1)
    
fit_zinck <- vb(
stan.model,
data = zinck_stan_data,
algorithm = "meanfield",
importance_resampling = TRUE,
iter = 10000,
tol_rel_obj = 0.01,
elbo_samples = 500
)
beta <- fit_zinck@sim[["est"]][["beta"]]
theta <- fit_zinck@sim[["est"]][["theta"]]
    
X_tilde <- generateKnockoff(X,theta,beta,seed=1)
index <- 1:30 ## Index set of the true non-zero signals
################### Varying the target FDR thresholds #########################
############## FDR = 0.05 ###############
index_est <- zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.05)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr1_LN_cts_list[i] <- FP/(FP+TP)
power1_LN_cts_list[i] <- TP/(TP+FN)
    
############### FDR = 0.1 ################
index_est <- zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.1)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr2_LN_cts_list[i] <- FP/(FP+TP)
power2_LN_cts_list[i] <- TP/(TP+FN)
    
################ FDR = 0.15 ###############
index_est <- zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.15)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr3_LN_cts_list[i] <- FP/(FP+TP)
power3_LN_cts_list[i] <- TP/(TP+FN)
    
################ FDR = 0.2 ################
index_est <- zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.2)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr4_LN_cts_list[i] <- FP/(FP+TP)
power4_LN_cts_list[i] <- TP/(TP+FN)

########## Generating binary response ############
set.seed(1)
signals = (2 * rbinom(30, 1, 0.5) - 1) * runif(30, 3, 10)
kBeta = c(signals / sqrt(col_abundances[1:30]), rep(0, p - 30))
pr = 1/(1+exp(-(Pi^2 %*% (kBeta/2) + Pi %*% kBeta)))
Y_bin = rbinom(n,1,pr)

################### Varying the target FDR thresholds #########################
############## FDR = 0.05 ###############
index_est <- zinck.filter(X,X_tilde,as.factor(Y_bin),model="Random Forest",fdr=0.05)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr1_LN_bin_list[i] <- FP/(FP+TP)
power1_LN_bin_list[i] <- TP/(TP+FN)
    
############### FDR = 0.1 ################
index_est <- zinck.filter(X,X_tilde,as.factor(Y_bin),model="Random Forest",fdr=0.1)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr2_LN_bin_list[i] <- FP/(FP+TP)
power2_LN_bin_list[i] <- TP/(TP+FN)
    
################ FDR = 0.15 ###############
index_est <- zinck.filter(X,X_tilde,as.factor(Y_bin),model="Random Forest",fdr=0.15)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr3_LN_bin_list[i] <- FP/(FP+TP)
power3_LN_bin_list[i] <- TP/(TP+FN)
    
################ FDR = 0.2 ################
index_est <- zinck.filter(X,X_tilde,as.factor(Y_bin),model="Random Forest",fdr=0.2)
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
    
fdr4_LN_bin_list[i] <- FP/(FP+TP)
power4_LN_bin_list[i] <- TP/(TP+FN)
}

########## Mean Empirical FDRs for different thresholds #############
print(mean(fdr1_DM_cts_list)) 
print(mean(fdr1_LN_cts_list)) 
print(mean(fdr1_DM_bin_list)) 
print(mean(fdr1_LN_bin_list)) 
print(mean(fdr2_DM_cts_list)) 
print(mean(fdr2_LN_cts_list)) 
print(mean(fdr2_DM_bin_list)) 
print(mean(fdr2_LN_bin_list)) 
print(mean(fdr3_DM_cts_list)) 
print(mean(fdr3_LN_cts_list)) 
print(mean(fdr3_DM_bin_list)) 
print(mean(fdr3_LN_bin_list)) 
print(mean(fdr4_DM_cts_list)) 
print(mean(fdr4_LN_cts_list)) 
print(mean(fdr4_DM_bin_list)) 
print(mean(fdr4_LN_bin_list)) 

########## Mean detection powers for different thresholds #############
print(mean(power1_DM_cts_list)) 
print(mean(power1_LN_cts_list)) 
print(mean(power1_DM_bin_list)) 
print(mean(power1_LN_bin_list)) 
print(mean(power2_DM_cts_list)) 
print(mean(power2_LN_cts_list)) 
print(mean(power2_DM_bin_list)) 
print(mean(power2_LN_bin_list)) 
print(mean(power3_DM_cts_list)) 
print(mean(power3_LN_cts_list)) 
print(mean(power3_DM_bin_list)) 
print(mean(power3_LN_bin_list)) 
print(mean(power4_DM_cts_list)) 
print(mean(power4_LN_cts_list)) 
print(mean(power4_DM_bin_list)) 
print(mean(power4_LN_bin_list)) 

########################################### Comparing other methods ########################################################
############################################################################################################################

## The methods under comparison are MX-KF (Model-X Knockoffs) [Candes et al. (2018)], LDA-KF (vanilla LDA Knockoffs), 
## DeepLINK [Zhu et al. (2021)], and CKF (Compositional Knockoff Filter) [Srinivasan et al. (2021)].

## You can run the same simulation settings by just replacing the "Fitting the zinck model" chunk using the corresponding code blocks --

######################## Fitting MX-KF ##########################

##### For both DM & LN generation ######
W <- log_normalize(X)
set.seed(1)
kfp = knockoff.filter(X=W,y=Y,fdr = FDR,statistic = stat.random_forest,offset=0) ## Change FDR accordingly!
kfStat = kfp$statistic
t = knockoff.threshold(kfStat, fdr = FDR,offset=1)
kfSelect = sort(which(kfStat >= t))
index_est = kfSelect

######################## Fitting LDA-KF ##########################

##### For both DM & LN generation #####
df.LDA = as(as.matrix(X),"dgCMatrix")
set.seed(1)
vanilla.LDA <- LDA(df.LDA,k=8,method="VEM") ## k = 8 (for p = 100), k = 10 (for p = 200), k = 11 (for p = 300) and k = 12 (for p = 400) 
theta.LDA <- vanilla.LDA@gamma
beta.LDA <- vanilla.LDA@beta
beta.LDA <- t(apply(beta.LDA,1,function(row) row/sum(row)))
X_tilde <- zinck::generateKnockoff(X,theta.LDA,beta.LDA,seed=1) ## Generating vanilla LDA knockoff copy
W_tilde <- log_normalize(X_tilde) ## Making the knockoff copy compositional
W <- stat.random_forest(W,W_tilde,Y)
T <- knockoff.threshold(W,fdr=FDR,offset=1) ## Vary FDR accordingly!
index_est <- which(W>=T) 


#################### Fitting DeepLINK #########################

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
q = FDR
		    
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

############## Continuous Outcomes ############

Xnew = dl.knockoff_construct(X1, r_hat, 'elu', aut_epoch, aut_loss, aut_verb)

  # compute knockoff statistics

p = Xnew.shape[1] // 2

 # implement DeepPINK
es = EarlyStopping(monitor='val_loss', patience=30, verbose=2)
dp = Sequential()
dp.add(PairwiseConnected(input_shape=(2 * p,)))
dp.add(Dense(p, activation='elu', kernel_regularizer=keras.regularizers.l1(l1=l1)))
dp.add(Dense(1, activation=None))
dp.compile(loss=mlp_loss, optimizer=keras.optimizers.Adam(learning_rate=lr))
dp.fit(Xnew, Y, epochs=mlp_epoch, batch_size=32, verbose=mlp_verb, validation_split=0.1, callbacks=[es])

 # calculate knockoff statistics W_j
weights = dp.get_weights()
w = weights[1] @ weights[3]
w = w.reshape(p, )
z = weights[0][:p]
z_tilde = weights[0][p:]
W = (w * z) ** 2 - (w * z_tilde) ** 2
  # feature selection

selected = dl.knockoff_select(W, q, ko_plus=False)
selected_plus = dl.knockoff_select(W, q, ko_plus=True)

############## Binary Outcomes ################
mlp_loss = 'binary_crossentropy'
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
dp.fit(Xnew, Y_bin, epochs=mlp_epoch, batch_size=32, verbose=mlp_verb, validation_split=0.1, callbacks=[es])

 # calculate knockoff statistics W_j
weights = dp.get_weights()
w = weights[1] @ weights[3]
w = w.reshape(p, )
z = weights[0][:p]
z_tilde = weights[0][p:]
W = (w * z) ** 2 - (w * z_tilde) ** 2
  # feature selection
selected = dl.knockoff_select(W, q, ko_plus=False)
selected_plus = dl.knockoff_select(W, q, ko_plus=True)
```
##################### Fitting CKF ####################

#source("utilityFnCKF.R")  ## Obtained from the CKF software file
#source("functionsCKF.R")
library(gurobi)

n0 = 100
n1 = n - n0
nScreen = 40

W=X
Z = log_normalize(W)
noPsudeoCountW = W
while(anotherZ == T){      
relativeZ = Z
### sample the sets
sampleCol = sample(1:n, size = n0, replace = F)
relativeZ0 = relativeZ[sampleCol,]
relativeZ1 = relativeZ[-sampleCol,]
Xscreen = X[sampleCol,]
Yscreen = Y[sampleCol]
X0 = Xscreen
Y0 = Yscreen
X1 = X[-sampleCol,]
Y1 = Y[-sampleCol]
W0 = W[sampleCol,]
W1 = W[-sampleCol,]
Z0dupe = sum(duplicated(t(relativeZ0)))
Z1dupe = sum(duplicated(t(relativeZ1)))
duplicates = Z0dupe+Z1dupe
if(duplicates == 0){
anotherZ = F
}
}
########## MIO #########
mioFit = screeningStep(relativeZ0, Y0, nScreen, n0 = n0, n1 = n1, 60*10)
mioSel = mioFit[[1]]
    
##### MIO Fits ####
mSelect = selectionStep(X,Y,X0,Y0,X1,Y1,W0,W1,mioSel,FDR)
index_est <- mSelect[[1]]
# For DeepLINK, please refer to the software written in Python 3.7.6 named DeepLINK publicly available in the GitHub repository : https://github.com/zifanzhu/DeepLINK
# For CKF, please refer to the software published in the supplementary material of the paper: "Compositional knockoff filter for high-dimensional regression analysis of microbiome data" 
# by Srinivasan et al. (2021)
# Note that CKF requires a Gurobi license.
