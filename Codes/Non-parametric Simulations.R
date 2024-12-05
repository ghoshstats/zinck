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
#################################### Working with CRC species level data #############################

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
return(list(Y = Y, X = X, Y_bin = Y_bin, index = 1:30))
}

#niter = 100
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
    X1 <- generate_data(p=ntaxa, seed=i)$X1
    Y1 <- generate_data(p=ntaxa, seed=i)$Y1

    ####################### Continuous Outcomes ##############################

    ################################### Fitting the zinck model ################################################
    ####### Initializing Delta for ADVI ########
    
    dlt <- rep(0,ncol(X1))
    for(t in (1:ncol(X1)))
    {
      dlt[t] <- 1-mean(X1[,t]>0)
    }
    
    zinck_stan_data <- list(
      K = 15,    ### NOTE: K = 15 is optimal for p = 100, K = 16 is optimal for p = 200, K = 17 is optimal for p = 300 and K = 18 is optimal for p = 400.
      V = ncol(X1),
      D = nrow(X1),
      n = X1,
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
    
    W_tilde1 <- log_normalize(X1_tilde) ## Making the knockoff copy compositional
    #############################################################################################################
    index <- generate_data(p=ntaxa,seed=i)$index ## Index set of the true non-zero signals
    
    ################### Varying the target FDR thresholds #########################
    
    ############## FDR = 0.05 ###############
    index_est <- zinck.filter(W1,W_tilde1,Y1,model="glmnet",fdr=0.05)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr1_cts_list[i] <- FP/(FP+TP)
    power1_cts_list[i] <- TP/(TP+FN)
    
    ############### FDR = 0.1 ################
    index_est <- zinck.filter(W1,W_tilde1,Y1,model="glmnet",fdr=0.1)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr2_cts_list[i] <- FP/(FP+TP)
    power2_cts_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.15 ###############
    index_est <- zinck.filter(W1,W_tilde1,Y1,model="glmnet",fdr=0.15)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr3_cts_list[i] <- FP/(FP+TP)
    power3_cts_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.2 ################
    index_est <- zinck.filter(W1,W_tilde1,Y1,model="glmnet",fdr=0.2)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr4_cts_list[i] <- FP/(FP+TP)
    power4_cts_list[i] <- TP/(TP+FN)
    
    
    ############################ Binary Outcomes #################################
    
   Y1_bin <- generate_data(p = ntaxa, seed = i)$Y1_bin
    ############## FDR = 0.05 ###############
    index_est <- zinck.filter(W1,W_tilde1,Y1_bin,model="glmnet",fdr=0.05)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr1_bin_list[i] <- FP/(FP+TP)
    power1_bin_list[i] <- TP/(TP+FN)
    
    ############### FDR = 0.1 ################
    index_est <- zinck.filter(W1,W_tilde1,Y1_bin,model="glmnet",fdr=0.1)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr2_bin_list[i] <- FP/(FP+TP)
    power2_bin_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.15 ###############
    index_est <- zinck.filter(W1,W_tilde1,Y1_bin,model="glmnet",fdr=0.15)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr3_bin_list[i] <- FP/(FP+TP)
    power3_bin_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.2 ################
    index_est <- zinck.filter(W1,W_tilde1,Y1_bin,model="glmnet",fdr=0.2)
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
set.seed(1)
W1 = log_normalize(X1)
kfp = knockoff.filter(X=W1,y=Y1,fdr = FDR,statistic = stat.glmnet_lambdasmax) ## Change FDR accordingly!
kfStat = kfp$statistic
t = knockoff.threshold(kfStat, fdr = FDR)
kfSelect = sort(which(kfStat >= t))
index_est = kfSelect

##### For binary outcomes ######
set.seed(1)
kfp = knockoff.filter(X=W1,y=Y1_bin,fdr = FDR,statistic = stat.lasso_coefdiff_bin) ## Change FDR accordingly!
kfStat = kfp$statistic
t = knockoff.threshold(kfStat, fdr = FDR)
kfSelect = sort(which(kfStat >= t))
index_est = kfSelect

######################## Fitting LDA-KF ##########################

##### For continuous outcomes #####
df.LDA = as(as.matrix(X1),"dgCMatrix")
set.seed(1)
vanilla.LDA <- LDA(df.LDA,k=6,method="VEM") ## k = 6 (for p = 100, 200, 300) and k = 8 (for p = 400) 
theta.LDA <- vanilla.LDA@gamma
beta.LDA <- vanilla.LDA@beta
beta.LDA <- t(apply(beta.LDA,1,function(row) row/sum(row)))
X1_tilde <- zinck::generateKnockoff(X1,theta.LDA,beta.LDA,seed=1) ## Generating vanilla LDA knockoff copy
W_tilde1 <- log_normalize(X1_tilde) ## Making the knockoff copy compositional
W <- stat.glmnet_coefdiff(W1,W_tilde1,Y1)
T <- knockoff.threshold(W,fdr=FDR,offset=1) ## Vary FDR accordingly!
index_est <- which(W>=T) 

##### For binary outcomes #####
W <- stat.lasso_coefdiff_bin(W1,W_tilde1,Y1_bin)
T <- knockoff.threshold(W,fdr=FDR,offset=1) ## Vary FDR accordingly!
index_est <- which(W>=T) 

# For DeepLINK, please refer to the software written in Python 3.7.6 named DeepLINK publicly available in the GitHub repository : https://github.com/zifanzhu/DeepLINK
# For CKF, please refer to the software published in the supplementary material of the paper: "Compositional knockoff filter for high-dimensional regression analysis of microbiome data" 
# by Srinivasan et al. (2021)
# Note that CKF requires a Gurobi license.
              

                
              

