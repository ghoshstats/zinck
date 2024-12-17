############################ Empirical Power and FDR for different target FDR levels ####################################
#########################################################################################################################

load("count.Rdata")
load("meta.RData")
library(rstan)
library(glmnet)
library(knockoff)
library(zinck)
library(kosel)
library(dplyr)
library(reshape2)
library(topicmodels)

######################### Stan code for fitting our zinck model ###########################

zinck_code <- "data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num words
  int<lower=0> D; // num docs
  int<lower=0> n[D, V]; // word counts for each doc

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] gamma1;
  vector<lower=0>[V] gamma2;
  vector<lower=0, upper=1>[V] delta;
}
parameters {
  simplex[K] theta[D]; // topic mixtures
  vector<lower=0,upper=1>[V] zeta[K]; // zero-inflated betas
}


transformed parameters {
  vector<lower=0>[V] beta[K];
  for (k in 1:K) {
	beta[k,1] =  zeta[k,1];
  for (m in 2:V) {
    beta[k,m] = zeta[k,m]*prod(1 - zeta[k,1:(m - 1)]);  // stick breaking
  }
  }
  for (k in 1:K) {
      beta[k]=beta[k]/sum(beta[k,1:V]);  // GD construction
  }
}


model {
  for (d in 1:D) {
    theta[d] ~ dirichlet(alpha);  
  }

  for (k in 1:K) {
    for (m in 1:V) {
      if (zeta[k,m]==0){  // Zero-inflated beta likelihood
        target += bernoulli_lpmf(1 | delta[m]);
      }else{
        target += bernoulli_lpmf(0 | delta[m]) + beta_lpdf(zeta[k,m] | gamma1[m], gamma2[m]);
      }
		}
  }

  for (d in 1:D) {
    vector[V] eta;
    eta = beta[1] * theta[d, 1];
    for (k in 2:K) {
      eta = eta + beta[k] * theta[d, k];
    }
    eta = eta/sum(eta[1:V]);
    n[d] ~ multinomial(eta);  // generation of each sample
  }
}
"
stan.model = stan_model(model_code = zinck_code)  ## Reading the stan model

######################################### Replicating the simulation setup ##########################################################

generate_data <- function(p,seed){
  dcount <- count[,order(decreasing=T,colSums(count,na.rm=T),apply(count,2L,paste,collapse=''))] ## ordering the columns w/ decreasing abundance
  ####### Randomly sampling patients from 574 observations #######
  set.seed(seed)
  norm_count <- count/rowSums(count)
  col_means <- colMeans(norm_count > 0)
  indices <- which(col_means > 0.2)
  sorted_indices <- indices[order(col_means[indices], decreasing=TRUE)]
  if(p %in% c(100,200,300,400)){
    dcount <- count[,sorted_indices][,1:p]
    sel_index <- sort(sample(1:nrow(dcount), 500))
    dcount <- dcount[sel_index,]
    original_OTU <- dcount + 0.5
    seq_depths <- rowSums(original_OTU)
    Pi = sweep(original_OTU, 1, seq_depths, "/")
    n = nrow(Pi)
    
    col_abundances = colMeans(Pi)
    
    ##### Generating continuous responses ######
    set.seed(1)
    signals = (2 * rbinom(30, 1, 0.5) - 1) * runif(30, 1.5, 3)
    kBeta = c(signals / sqrt(col_abundances[1:30]), rep(0, p - 30))
    eps=rnorm(n,mean = 0, sd=1)
    Y <- Pi^2 %*% (kBeta/2) + Pi %*% kBeta + eps
    
    ##### Generating binary responses #####
    set.seed(1)
    signals = (2 * rbinom(30, 1, 0.5) - 1) * runif(30, 3, 10)
    kBeta = c(signals / sqrt(col_abundances[1:30]), rep(0, p - 30))
    pr = 1/(1+exp(-(Pi^2 %*% (kBeta/2) + Pi %*% kBeta)))
    Y_bin = rbinom(n,1,pr)
    
    ######### Generate a copy of X #########
    X <- matrix(0, nrow = nrow(Pi), ncol = ncol(Pi))
    nSeq <- seq_depths
    # Loop over each row to generate the new counts based on the multinomial distribution
    
    set.seed(1)
    
    for (i in 1:nrow(Pi)) {
      X[i, ] <- rmultinom(1, size = nSeq, prob = Pi[i, ])
    }  
  } else if (p %in% c(200,300,400)) {
    print("Enter p within 100 to 400")
  }
  colnames(X) <- colnames(Pi)
  return(list(Y = Y, X = X, Y_bin = Y_bin, index = 1:30))
}

ntaxa = 100 # Change to p = 200, 300, 400 accordingly.

### Index 1 corresponds to FDR=0.05, 2 corresponds to FDR=0.1, 3 corresponds to FDR=0.15, and 4 corresponds to FDR=0.2 ###

power1_cts_list <- c()
power2_cts_list <- c()
power3_cts_list <- c()
power4_cts_list <- c()

power1_bin_list <- c()
power2_bin_list <- c()
power3_bin_list <- c()
power4_bin_list <- c()


fdr1_cts_list <- c()
fdr2_cts_list <- c()
fdr3_cts_list <- c()
fdr4_cts_list <- c()

fdr1_bin_list <- c()
fdr2_bin_list <- c()
fdr3_bin_list <- c()
fdr4_bin_list <- c()

for(i in 1:niter)
{
  tryCatch({
    X1 <- generate_data(p=ntaxa, seed=i)$X
    Y1 <- generate_data(p=ntaxa, seed=i)$Y

    ####################### Continuous Outcomes ##############################

    ################################### Fitting the zinck model ################################################
    ####### Initializing Delta for ADVI ########
    
   for(t in (1:ncol(X1)))
    {
    dlt[t] <- 1-mean(X1[,t]>0)
    }

    zinck_stan_data <- list(
    K = 15,
    V = ncol(X1),
    D = nrow(X1),
    n = X1,
    alpha = rep(0.1, 15),
    gamma1 = rep(0.5, ncol(X1)),
    gamma2 = rep(10,ncol(X1)),
    delta = dlt
    )
    
    set.seed(seed_list[i])
    fit1 <- vb(
      stan.model,
      data = zinck_stan_data,
      algorithm = "meanfield",
      importance_resampling = TRUE,
      iter = 10000,
      tol_rel_obj = 0.01,
      elbo_samples = 500
    )
    
    theta <- fit1@sim[["est"]][["theta"]]
    beta <- fit1@sim[["est"]][["beta"]]
    
    X1_tilde <- generateKnockoff(X1,theta,beta,seed=1) ## Generating the knockoff copy
  
    #############################################################################################################
    index <- 1:30 ## Index set of the true non-zero signals
    
    ############################ Fitting the Random Forest models with Y against X and X_tilde ####################
    
    ################### Varying the target FDR thresholds #########################

    ############## FDR = 0.05 ###############
    cts_rf <- suppressWarnings(zinck.filter(X1,X1_tilde,Y1,model="Random Forest",fdr=0.05,offset=0,mtry=200,seed=12,rftuning = TRUE))
    index_est <- cts_rf[["selected"]]
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr1_cts_list[i] <- FP/(FP+TP)
    power1_cts_list[i] <- TP/(TP+FN)
    
    ############### FDR = 0.1 ################
    cts_rf <- suppressWarnings(zinck.filter(X1,X1_tilde,Y1,model="Random Forest",fdr=0.1,offset=0,mtry=200,seed=12,rftuning = TRUE))
    index_est <- cts_rf[["selected"]]    
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr2_cts_list[i] <- FP/(FP+TP)
    power2_cts_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.15 ###############
    cts_rf <- suppressWarnings(zinck.filter(X1,X1_tilde,Y1,model="Random Forest",fdr=0.15,offset=0,mtry=200,seed=12,rftuning = TRUE))
    index_est <- cts_rf[["selected"]]    
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr3_cts_list[i] <- FP/(FP+TP)
    power3_cts_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.2 ################
    cts_rf <- suppressWarnings(zinck.filter(X1,X1_tilde,Y1,model="Random Forest",fdr=0.2,offset=0,mtry=200,seed=12,rftuning = TRUE))
    index_est <- cts_rf[["selected"]]     
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr4_cts_list[i] <- FP/(FP+TP)
    power4_cts_list[i] <- TP/(TP+FN)
    
    
    ############################ Binary Outcomes #################################
    
    Y1_bin <- generate_data(p = ntaxa, seed = i)$Y_bin
    ############## FDR = 0.05 ###############
    bin_rf <- suppressWarnings(zinck.filter(X1,X1_tilde,as.factor(Y1_bin),model="Random Forest",fdr=0.05,offset=0,mtry=14,seed=15,rftuning = TRUE))
    index_est <- bin_rf[["selected"]]
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr1_bin_list[i] <- FP/(FP+TP)
    power1_bin_list[i] <- TP/(TP+FN)
    
    ############### FDR = 0.1 ################
    bin_rf <- suppressWarnings(zinck.filter(X1,X1_tilde,as.factor(Y1_bin),model="Random Forest",fdr=0.1,offset=0,mtry=14,seed=15,rftuning = TRUE))
    index_est <- bin_rf[["selected"]]
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr2_bin_list[i] <- FP/(FP+TP)
    power2_bin_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.15 ###############
    bin_rf <- suppressWarnings(zinck.filter(X1,X1_tilde,as.factor(Y1_bin),model="Random Forest",fdr=0.15,offset=0,mtry=14,seed=15,rftuning = TRUE))
    index_est <- bin_rf[["selected"]]
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr3_bin_list[i] <- FP/(FP+TP)
    power3_bin_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.2 ################
    bin_rf <- suppressWarnings(zinck.filter(X1,X1_tilde,as.factor(Y1_bin),model="Random Forest",fdr=0.2,offset=0,mtry=14,seed=15,rftuning = TRUE))
    index_est <- bin_rf[["selected"]]
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr4_bin_list[i] <- FP/(FP+TP)
    power4_bin_list[i] <- TP/(TP+FN)
  }, error = function(e) {
    cat("An error occurred in iteration", i, "\n")
  })
}
########## Mean Empirical FDRs for different thresholds #############
print(mean(fdr1_cts_list)) 
print(mean(fdr1_bin_list)) 
print(mean(fdr2_cts_list)) 
print(mean(fdr2_bin_list)) 
print(mean(fdr3_cts_list)) 
print(mean(fdr3_bin_list)) 
print(mean(fdr4_cts_list)) 
print(mean(fdr4_bin_list)) 

########## Mean detection powers for different thresholds #############
print(mean(power1_cts_list)) 
print(mean(power1_bin_list)) 
print(mean(power2_cts_list)) 
print(mean(power2_bin_list)) 
print(mean(power3_cts_list)) 
print(mean(power3_bin_list)) 
print(mean(power4_cts_list)) 
print(mean(power4_bin_list)) 

########################################### Comparing other methods ########################################################
############################################################################################################################

## The methods under comparison are MX-KF (Model-X Knockoffs) [Candes et al. (2018)], LDA-KF (vanilla LDA Knockoffs), 
## DeepLINK [Zhu et al. (2021)], and CKF (Compositional Knockoff Filter) [Srinivasan et al. (2021)].

## You can run the same simulation settings by just replacing the "Fitting the zinck model" chunk using the corresponding code blocks --

######################## Fitting MX-KF ##########################

##### For continuous outcomes ######
W1 = log_normalize(X1)
W1_tilde = create.second_order(W1)
set.seed(3)
imp <- stat.random_forest(W1,W1_tilde,Y1)
T <- knockoff.threshold(imp,fdr=FDR)
index_est <- which(imp >=T)

##### For binary outcomes ######
set.seed(21)
imp <- stat.random_forest(W1,W1_tilde,as.factor(Y1_bin))
T <- knockoff.threshold(imp,fdr=FDR)
index_est <- which(imp >=T)

######################## Fitting LDA-KF ##########################

##### For continuous outcomes #####
df.LDA = as(as.matrix(X1),"dgCMatrix")
set.seed(1)
vanilla.LDA <- LDA(df.LDA,k=8,method="VEM") 
theta.LDA <- vanilla.LDA@gamma
beta.LDA <- vanilla.LDA@beta
beta.LDA <- t(apply(beta.LDA,1,function(row) row/sum(row)))
X1_tilde.LDA <- zinck::generateKnockoff(X1,theta.LDA,beta.LDA,seed=1) ## Generating vanilla LDA knockoff copy

set.seed(34)
imp <- stat.random_forest(X1,X1_tilde.LDA,Y1)
T <- knockoff.threshold(imp,fdr=FDR,offset=0)
index_est <- which(imp>=T)


##### For binary outcomes #####
set.seed(34)
imp <- stat.random_forest(X1,X1_tilde.LDA,as.factor(Y1_bin))
T <- knockoff.threshold(imp,fdr=FDR,offset=0)
index_est <- which(imp >=T)


 
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
		    
X1 -= np.mean(X1, axis=0)

# Prevent division by zero: Replace zero std with 1 before dividing
std_dev = np.std(X1, axis=0, ddof=1)
std_dev[std_dev == 0] = 1  # Avoid division by zero
X1 /= std_dev

# Normalize X1 while ensuring no division by zero in normalization
norms = np.sqrt(np.sum(X1 ** 2, axis=0))
norms[norms == 0] = 1  # Avoid division by zero
X2 = X1 / norms
r_hat = PCp1(X2, 15)

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
dp.fit(Xnew, Y1, epochs=mlp_epoch, batch_size=32, verbose=mlp_verb, validation_split=0.1, callbacks=[es])

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
dp.fit(Xnew, Y2, epochs=mlp_epoch, batch_size=32, verbose=mlp_verb, validation_split=0.1, callbacks=[es])

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
library(knockoff)
library(glmnet)
library(MethylCapSig)
library(energy)
library(dirmult)
library(gurobi)

n0 = 100
n1 = n - n0
nScreen = 40

W=X1
Y=Y1
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
              

