############################ Taxa Detection Power Simulation Setup ####################################
#######################################################################################################

load("count.Rdata")
load("meta.RData")
library(rstan)
library(glmnet)
library(knockoff)
library(zinck)
library(kosel)
library(dplyr)
library(reshape2)
############################# Single iteration of our simulation setup ###################################

####################### Working with CRC species level data ########################
dcount <- count[,order(decreasing=T,colSums(count,na.rm=T),apply(count,2L,paste,collapse=''))] ## ordering the columns w/ decreasing abundance
##################################################################################

####### Randomly sampling patients from 574 observations #######
set.seed(1)
sel_index <- rbinom(nrow(meta),size=1,prob=0.5)
selected_samples <- which(sel_index==1)
meta_selected <- meta[selected_samples,]
X <- dcount[selected_samples,]
mean(X>0) ### sparsity of X = 72%

##################################################################################
######################## Random Signal Injection #############################

p = 100
X1 <- X[,1:p]
n = nrow(X1)

Five1 = c(-3,3,2.5,-1, -1.5)
Five2 = c(3,3,-2,-2,-2)
Five3 = c(1,-1,3,-2,-1)
Five4 = c(-1,1,2,-1,-1)
Five5 = c(3,3,-3,-2,-1)
Five6 = c(-1,1,-2,1,1)
Five_all <- c(Five1,Five2,Five3,Five4,Five5,Five6)  ## Signals satisfy sum-to-zero constraint ##
randBeta <- rep(0,p)
set.seed(1)
rand_indices <- sample(1:p,size=30,replace=FALSE)  ## Randomly injecting 30 signals out of p=100 ##
set.seed(1)
randBeta[rand_indices] <- sample(Five_all,size=30,replace=FALSE)

W1 <- log_normalize(X1) ## log-normalizing counts to make it compositional

##### Generating continuous responses ######
set.seed(1)
eps=rnorm(n,mean = 0, sd=1)
Y1 <- W1 %*% randBeta + eps

##### Generating binary responses #####
set.seed(1)
pr = 1/(1+exp(-W1 %*% randBeta))
Y1_bin = rbinom(n,1,pr)


############################## p = 200 #############################
X2 <- X[,1:p]
n = nrow(X2)
Five1 = c(-3,3,2.5,-1, -1.5)
Five2 = c(3,3,-2,-2,-2)
Five3 = c(1,-1,3,-2,-1)
Five4 = c(-1,1,2,-1,-1)
Five5 = c(3,3,-3,-2,-1)
Five6 = c(-1,1,-2,1,1)
Five_all <- c(Five1,Five2,Five3,Five4,Five5,Five6)
randBeta <- rep(0,p)
set.seed(1)
rand_indices <- sample(1:200,size=30,replace=FALSE) ## Randomly inject 30 signals only among the first 200 taxa ##
set.seed(1)
randBeta[rand_indices] <- sample(Five_all,size=30,replace=FALSE)

W2 <- log_normalize(X2) ## log-normalizing counts to make them compositional

###### Generating continuous responses #######
set.seed(1)
eps=rnorm(n,mean = 0, sd=1)
Y2 <- W2 %*% randBeta + eps

####### Generating binary responses ########
set.seed(4)
pr = 1/(1+exp(-W2 %*% randBeta))
Y2_bin = rbinom(n,1,pr)

############################### p=300 ###################################
p = 300
X3 <- X[,1:p]
n = nrow(X3)
Five1 = c(-3,3,2.5,-1, -1.5)
Five2 = c(3,3,-2,-2,-2)
Five3 = c(1,-1,3,-2,-1)
Five4 = c(-1,1,2,-1,-1)
Five5 = c(3,3,-3,-2,-1)
Five6 = c(-1,1,-2,1,1)
Five_all <- c(Five1,Five2,Five3,Five4,Five5,Five6)
randBeta <- rep(0,p)
set.seed(1)
rand_indices <- sample(1:200,size=30,replace=FALSE) ## Randomly inject 30 signals only among the first 200 taxa ##
set.seed(1)
randBeta[rand_indices] <- sample(Five_all,size=30,replace=FALSE)

W3 <- log_normalize(X3) ## log-normalizing counts to make them compositional

###### Generating continuous responses #######
set.seed(1)
eps=rnorm(n,mean = 0, sd=1)
Y3 <- W3 %*% randBeta + eps

####### Generating binary responses ########
set.seed(1)
pr = 1/(1+exp(-W3 %*% randBeta))
Y3_bin = rbinom(n,1,pr)

########################### p=400 ################################

p = 400
X4 <- X[,1:p]
n = nrow(X4)
Five1 = c(-3,3,2.5,-1, -1.5)
Five2 = c(3,3,-2,-2,-2)
Five3 = c(1,-1,3,-2,-1)
Five4 = c(-1,1,2,-1,-1)
Five5 = c(3,3,-3,-2,-1)
Five6 = c(-1,1,-2,1,1)
Five_all <- c(Five1,Five2,Five3,Five4,Five5,Five6)
randBeta <- rep(0,p)
set.seed(1)
rand_indices <- sample(1:200,size=30,replace=FALSE) ## Randomly inject 30 signals only among the first 200 taxa ##
set.seed(1)
randBeta[rand_indices] <- sample(Five_all,size=30,replace=FALSE)

W4 <- log_normalize(X4) ## log-normalizing counts to make them compositional

###### Generating continuous responses #######
set.seed(1)
eps=rnorm(n,mean = 0, sd=1)
Y4 <- W4 %*% randBeta + eps

####### Generating binary responses ########
set.seed(1)
pr = 1/(1+exp(-W4 %*% randBeta))
Y4_bin = rbinom(n,1,pr)

############################### Fitting the zinck model ####################################
############################################################################################
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

####################### Fitting zinck for p=100 ##############################

#### Case (1): Continuous responses #####
######## Initializing Delta for ADVI #########
dlt <- rep(0,ncol(X1))
for(t in (1:ncol(X1)))
{
  dlt[t] <- 1-mean(X1[,t]>0)
}
######### Input data for variational Bayes ###########
zinLDA_stan_data <- list(
  K = 8,
  V = ncol(X1),
  D = nrow(X1),
  n = X1,
  alpha = rep(0.1, 8),
  gamma1 = rep(0.5, ncol(X1)),
  gamma2 = rep(10,ncol(X1)),
  delta = dlt
)

set.seed(1)  ## Very sensitive to starting values!
fit1 <- vb(stan.model, data=zinLDA_stan_data, algorithm="meanfield", iter=10000,tol_rel_obj=0.01)  ### Fitting the zinck model ###

#### Posterior Estimates #####
theta <- fit1@sim[["est"]][["theta"]] 
beta <- fit1@sim[["est"]][["beta"]]  

####### Generating the knockoff copy #########
X1_tilde <- generateKnockoff(X1,theta,beta,seed=1)
W_tilde1 <- log_normalize(X1_tilde)

selected_species = zinck.filter(W1,W_tilde1,Y1,model="glmnet",fdr=0.1,offset=1) ## Selected taxa

###### Estimated Power and FDR #######

index <- rand_indices
index_est <- selected_species
neg_index <- (1:100)[-c(index)]

if(length(index_est)==0){
  neg_index_est <- 1:100
} else{
  neg_index_est <- (1:100)[-c(index_est)]
}

### Evaluation metrics ###
TP <- sum(index_est %in% index==TRUE) ## True Positives
FP <- sum(index_est %in% index==FALSE) ## False Positives
FN <- sum(index %in% index_est == FALSE) ## False Negatives
TN <- sum(neg_index_est %in% neg_index == TRUE) ## True Negatives
FDR <- FP/(FP+TP)
Power <- TP/(TP+FN)

###### Case (2): Binary Responses ######
########################################

W <- stat.lasso_coefdiff_bin(W1,W_tilde1,Y1_bin)
msel <- ko.sel(W, print = FALSE)
index_est <- which(W>=msel$threshold)

neg_index <- (1:100)[-c(rand_indices)]
if(length(index_est) == 0) {
  neg_index_est <- 1:100
} else {
  neg_index_est <- (1:100)[-c(index_est)]
}
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
TN <- sum(neg_index_est %in% neg_index == TRUE)
FDR <- FP/(FP+TP)
Power <- TP/(TP+FN)

####################### Fitting zinck for p=200 ##############################

#### Case (1): Continuous responses #####
######## Initializing Delta for ADVI #########
dlt <- rep(0,ncol(X2))
for(t in (1:ncol(X2)))
{
  dlt[t] <- 1-mean(X2[,t]>0)
}
######### Input data for variational Bayes ###########
zinLDA_stan_data <- list(
  K = 16,
  V = ncol(X2),
  D = nrow(X2),
  n = X2,
  alpha = rep(0.1, 16),
  gamma1 = rep(0.5, ncol(X2)),
  gamma2 = rep(10,ncol(X2)),
  delta = dlt
)

set.seed(1)  ## Very sensitive to starting values!
fit2 <- vb(stan.model, data=zinLDA_stan_data, algorithm="meanfield", iter=10000,tol_rel_obj=0.01)  ### Fitting the zinck model ###

#### Posterior Estimates #####
theta <- fit2@sim[["est"]][["theta"]] 
beta <- fit2@sim[["est"]][["beta"]]  

####### Generating the knockoff copy #########
X2_tilde <- generateKnockoff(X2,theta,beta,seed=1)
W_tilde2 <- log_normalize(X2_tilde)

selected_species = zinck.filter(W2,W_tilde2,Y2,model="glmnet",fdr=0.1,offset=1) ## Selected taxa

###### Estimated Power and FDR #######

index <- rand_indices
index_est <- selected_species
neg_index <- (1:200)[-c(index)]

if(length(index_est)==0){
  neg_index_est <- 1:200
} else{
  neg_index_est <- (1:200)[-c(index_est)]
}

### Evaluation metrics ###
TP <- sum(index_est %in% index==TRUE) ## True Positives
FP <- sum(index_est %in% index==FALSE) ## False Positives
FN <- sum(index %in% index_est == FALSE) ## False Negatives
TN <- sum(neg_index_est %in% neg_index == TRUE) ## True Negatives
FDR <- FP/(FP+TP)
Power <- TP/(TP+FN)

###### Case (2): Binary Responses ######
########################################

W <- stat.lasso_coefdiff_bin(W2,W_tilde2,Y2_bin)
msel <- ko.sel(W, print = FALSE, method="gaps")
index_est <- which(W>=msel$threshold)

neg_index <- (1:200)[-c(rand_indices)]
if(length(index_est) == 0) {
  neg_index_est <- 1:200
} else {
  neg_index_est <- (1:200)[-c(index_est)]
}
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
TN <- sum(neg_index_est %in% neg_index == TRUE)
FDR <- FP/(FP+TP)
Power <- TP/(TP+FN)

####################### Fitting zinck for p=300 ##############################

#### Case (1): Continuous responses #####
######## Initializing Delta for ADVI #########
dlt <- rep(0,ncol(X3))
for(t in (1:ncol(X3)))
{
  dlt[t] <- 1-mean(X3[,t]>0)
}
######### Input data for variational Bayes ###########
zinLDA_stan_data <- list(
  K = 12,
  V = ncol(X3),
  D = nrow(X3),
  n = X3,
  alpha = rep(0.1, 12),
  gamma1 = rep(0.5, ncol(X3)),
  gamma2 = rep(10,ncol(X3)),
  delta = dlt
)

set.seed(12)  ## Very sensitive to starting values!
fit3 <- vb(stan.model, data=zinLDA_stan_data, algorithm="meanfield", iter=10000,tol_rel_obj=0.01)  ### Fitting the zinck model ###

#### Posterior Estimates #####
theta <- fit3@sim[["est"]][["theta"]] 
beta <- fit3@sim[["est"]][["beta"]]  

####### Generating the knockoff copy #########
X3_tilde <- generateKnockoff(X3,theta,beta,seed=1)
W_tilde3 <- log_normalize(X3_tilde)

selected_species = zinck.filter(W3,W_tilde3,Y3,model="glmnet",fdr=0.1,offset=1) ## Selected taxa

###### Estimated Power and FDR #######

index <- rand_indices
index_est <- selected_species
neg_index <- (1:300)[-c(index)]

if(length(index_est)==0){
  neg_index_est <- 1:300
} else{
  neg_index_est <- (1:300)[-c(index_est)]
}

### Evaluation metrics ###
TP <- sum(index_est %in% index==TRUE) ## True Positives
FP <- sum(index_est %in% index==FALSE) ## False Positives
FN <- sum(index %in% index_est == FALSE) ## False Negatives
TN <- sum(neg_index_est %in% neg_index == TRUE) ## True Negatives
FDR <- FP/(FP+TP)
Power <- TP/(TP+FN)

###### Case (2): Binary Responses ######
########################################

W <- stat.lasso_coefdiff_bin(W3,W_tilde3,Y3_bin)
msel <- ko.sel(W, print = FALSE, method="gaps")
index_est <- which(W>=msel$threshold)

neg_index <- (1:300)[-c(rand_indices)]
if(length(index_est) == 0) {
  neg_index_est <- 1:300
} else {
  neg_index_est <- (1:300)[-c(index_est)]
}
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
TN <- sum(neg_index_est %in% neg_index == TRUE)
FDR <- FP/(FP+TP)
Power <- TP/(TP+FN)

####################### Fitting zinck for p=400 ##############################

#### Case (1): Continuous responses #####
######## Initializing Delta for ADVI #########
dlt <- rep(0,ncol(X4))
for(t in (1:ncol(X4)))
{
  dlt[t] <- 1-mean(X4[,t]>0)
}
######### Input data for variational Bayes ###########
zinLDA_stan_data <- list(
  K = 15,
  V = ncol(X4),
  D = nrow(X4),
  n = X4,
  alpha = rep(0.1, 15),
  gamma1 = rep(0.5, ncol(X4)),
  gamma2 = rep(10,ncol(X4)),
  delta = dlt
)

set.seed(12)  ## Very sensitive to starting values!
fit4 <- vb(stan.model, data=zinLDA_stan_data, algorithm="meanfield", iter=10000,tol_rel_obj=0.01)  ### Fitting the zinck model ###

#### Posterior Estimates #####
theta <- fit4@sim[["est"]][["theta"]] 
beta <- fit4@sim[["est"]][["beta"]]  

####### Generating the knockoff copy #########
X4_tilde <- generateKnockoff(X4,theta,beta,seed=1)
W_tilde4 <- log_normalize(X4_tilde)

selected_species = zinck.filter(W4,W_tilde4,Y4,model="glmnet",fdr=0.1,offset=1) ## Selected taxa

###### Estimated Power and FDR #######

index <- rand_indices
index_est <- selected_species
neg_index <- (1:400)[-c(index)]

if(length(index_est)==0){
  neg_index_est <- 1:400
} else{
  neg_index_est <- (1:400)[-c(index_est)]
}

### Evaluation metrics ###
TP <- sum(index_est %in% index==TRUE) ## True Positives
FP <- sum(index_est %in% index==FALSE) ## False Positives
FN <- sum(index %in% index_est == FALSE) ## False Negatives
TN <- sum(neg_index_est %in% neg_index == TRUE) ## True Negatives
FDR <- FP/(FP+TP)
Power <- TP/(TP+FN)

###### Case (2): Binary Responses ######
########################################

W <- stat.lasso_coefdiff_bin(W4,W_tilde4,Y4_bin)
msel <- ko.sel(W, print = FALSE, method="gaps")
index_est <- which(W>=msel$threshold)

neg_index <- (1:400)[-c(rand_indices)]
if(length(index_est) == 0) {
  neg_index_est <- 1:400
} else {
  neg_index_est <- (1:400)[-c(index_est)]
}
FN <- sum(index %in% index_est == FALSE) 
FP <- sum(index_est %in% index == FALSE)
TP <- sum(index_est %in% index == TRUE)
TN <- sum(neg_index_est %in% neg_index == TRUE)
FDR <- FP/(FP+TP)
Power <- TP/(TP+FN)















