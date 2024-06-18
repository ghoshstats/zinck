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

#################################### Working with CRC species level data #############################

generate_data <- function(p,seed){
  dcount <- count[,order(decreasing=T,colSums(count,na.rm=T),apply(count,2L,paste,collapse=''))] ## ordering the columns w/ decreasing abundance
  ####### Randomly sampling patients from 574 observations #######
  set.seed(seed)
  sel_index <- rbinom(nrow(meta),size=1,prob=0.5)
  selected_samples <- which(sel_index==1)
  meta_selected <- meta[selected_samples,]
  X <- dcount[selected_samples,]
  if(p == 100){
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
  } else if (p %in% c(200,300,400)) {
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
    rand_indices <- sample(1:200,size=30,replace=FALSE)  ## Randomly injecting 30 signals out of p=200 ##
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
  }

  return(list(Y1 = Y1, X1 = X1,  W1 = W1, Y1_bin = Y1_bin, index = rand_indices))
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
    W1 <- generate_data(p=ntaxa, seed=i)$W1
    Y1 <- generate_data(p=ntaxa, seed=i)$Y1
    
    ################ Fitting the zinck model ####################
    
    ####################### Continuous Outcomes ##############################
    
    ####### Initializing Delta for ADVI ########
    
    dlt <- rep(0,ncol(X1))
    for(t in (1:ncol(X1)))
    {
      dlt[t] <- 1-mean(X1[,t]>0)
    }
    
    zinck_stan_data <- list(
      K = 15,
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







