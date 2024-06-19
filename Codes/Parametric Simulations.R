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

load("DirMultOutput.RData")

n = 250
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

########## Function to generate Dirichlet Multinomial (DM) and Logistic Normal (LN) data ############

generate_data_parametric <- function(p,type,seed){
  ### Beta Vector
  Five1 = c(-3,3,2.5,-1, -1.5)
  Five2 = c(3,3,-2,-2,-2)
  Five3 = c(1,-1,3,-2,-1)
  Five4 = c(-1,1,2,-1,-1)
  Five5 = c(3,3,-3,-2,-1)
  Five6 = c(-1,1,-2,1,1)
  kBeta = c(Five1, Five2, Five3, Five4, Five5, Five6,rep(0,p - 30))
  if(type=="DM"){
    data(throat.tree)
    data(throat.otu.tab)
    tree = throat.tree 
    tree.dist = cophenetic(tree)
    otu.ids <- tree$tip.label

    set.seed(seed)
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
    nSeq <- rnbinom(n , mu = 10000, size = 25)
    for (i in 1:n ) {
      comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
      comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
    }
    colindex=order(p.est, decreasing = T)[1:p]
    X = comm[,colindex]
    W = log_normalize(X)
    ##### Generating continuous responses ######
    set.seed(1)
    eps=rnorm(n,mean = 0, sd=1)
    Y <- W %*% kBeta + eps
  } else if (type=="LN"){
    ### Generate covariance matrix ###
    sigma = matrix(0,p,p)
    gamma = 0.5
    for(i in 1:nrow(sigma)){
      for(j in 1:nrow(sigma)){
        sigma[i,j] = gamma^(abs(i-j))   
      }
    }
    ### Generate mean ###
    means = rep(1,p)

    W = mvlognormal(n, means, Sigma = diag(sigma), R = sigma,seed=seed)
    Z = matrix(0, nrow(W), ncol(W))
    X = matrix(0, nrow(W), ncol(W))
    
    set.seed(1)
    N <- rnbinom(n,size=25,mu=8000) ## Library sizes
    
    rowW = rowSums(W)
    colW = colSums(W)
    Z <- sweep(W, 1, rowW, "/")
    X <- floor(N * Z)
    colnames(X) <-  c(paste0("Taxa", 1:p))
    colnames(W) <-  c(paste0("Taxa", 1:p))
    ##### Generating continuous responses ######
    set.seed(1)
    eps=rnorm(n,mean = 0, sd=1)
    Y <- W %*% kBeta + eps
  }
  return(list(Y = Y, X = X,  W = W))
}

power1_DM_list <- c()
power2_DM_list <- c()
power3_DM_list <- c()
power4_DM_list <- c()

fdr1_DM_list <- c()
fdr2_DM_list <- c()
fdr3_DM_list <- c()
fdr4_DM_list <- c()

power1_LN_list <- c()
power2_LN_list <- c()
power3_LN_list <- c()
power4_LN_list <- c()

fdr1_LN_list <- c()
fdr2_LN_list <- c()
fdr3_LN_list <- c()
fdr4_LN_list <- c()

for(i in 1:niter)
{
  tryCatch({
    
    ################### DM generation ######################
    X <- generate_data_parametric(p=100, type="DM", seed=i)$X
    W <- generate_data_parametric(p=100, type="DM", seed=i)$W
    Y <- generate_data_parametric(p=100, type="DM", seed=i)$Y
    
    ################ Fitting the zinck model ####################
    
    ####### Initializing Delta for ADVI ########
    
    dlt <- rep(0,ncol(X))
    
    for(t in (1:ncol(X)))
    {
      dlt[t] <- 1-mean(X[,t]>0)
    }
    
    zinck_stan_data <- list(
      K = 8, ## k = 8 (for p = 100), k = 10 (for p = 200), k = 12 (for p = 300) and k = 15 (for p = 400) 
      V = ncol(X),
      D = nrow(X),
      n = X,
      delta = dlt
    )
  
    set.seed(seed_list_DM[i])
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
    X_tilde <- zinck::generateKnockoff(X,theta,beta,seed=1)
    W_tilde <- log_normalize(X_tilde)

    index <- 1:30 ## Index set of the true non-zero signals
    
    ################### Varying the target FDR thresholds #########################
    
    ############## FDR = 0.05 ###############
    index_est <- zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.05)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr1_DM_list[i] <- FP/(FP+TP)
    power1_DM_list[i] <- TP/(TP+FN)
    
    ############### FDR = 0.1 ################
    index_est <- zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.1)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr2_DM_list[i] <- FP/(FP+TP)
    power2_DM_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.15 ###############
    index_est <- zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.15)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr3_DM_list[i] <- FP/(FP+TP)
    power3_DM_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.2 ################
    index_est <- zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.2)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr4_DM_list[i] <- FP/(FP+TP)
    power4_DM_list[i] <- TP/(TP+FN)
    
    
    ############################ LN generation #################################
    
    
    X <- generate_data_parametric(p=100, type="LN", seed=i)$X
    W <- generate_data_parametric(p=100, type="LN", seed=i)$W
    Y <- generate_data_parametric(p=100, type="LN", seed=i)$Y
    
    ################ Fitting the zinck model ####################
    
    ####### Initializing Delta for ADVI ########
    
    dlt <- rep(0,ncol(X))
    
    for(t in (1:ncol(X)))
    {
      dlt[t] <- 1-mean(X[,t]>0)
    }
    
    zinck_stan_data <- list(
      K = 8,      ## k = 8 (for p = 100), k = 10 (for p = 200), k = 11 (for p = 300) and k = 12 (for p = 400) 
      V = ncol(X),
      D = nrow(X),
      n = X,
      delta = dlt
    )
    
    stan.model <- stan_model(model_code = zinck_code)
    set.seed(seed_list_LN[i])
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
    X_tilde <- zinck::generateKnockoff(X,theta,beta,seed=1)
    W_tilde <- log_normalize(X_tilde)
    
    index <- 1:30 ## Index set of the true non-zero signals
    
    ################### Varying the target FDR thresholds #########################
    
    ############## FDR = 0.05 ###############
    index_est <- zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.05)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr1_LN_list[i] <- FP/(FP+TP)
    power1_LN_list[i] <- TP/(TP+FN)
    
    ############### FDR = 0.1 ################
    index_est <- zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.1)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr2_LN_list[i] <- FP/(FP+TP)
    power2_LN_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.15 ###############
    index_est <- zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.15)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr3_LN_list[i] <- FP/(FP+TP)
    power3_LN_list[i] <- TP/(TP+FN)
    
    ################ FDR = 0.2 ################
    index_est <- zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.2)
    FN <- sum(index %in% index_est == FALSE) 
    FP <- sum(index_est %in% index == FALSE)
    TP <- sum(index_est %in% index == TRUE)
    
    fdr4_LN_list[i] <- FP/(FP+TP)
    power4_LN_list[i] <- TP/(TP+FN)
    
    
  }, error = function(e) {
    cat("An error occurred in iteration", i, "\n")
  })
}
########## Mean Empirical FDRs for different thresholds #############
print(mean(fdr1_DM_list)) 
print(mean(fdr1_LN_list)) 
print(mean(fdr2_DM_list)) 
print(mean(fdr2_LN_list)) 
print(mean(fdr3_DM_list)) 
print(mean(fdr3_LN_list)) 
print(mean(fdr4_DM_list)) 
print(mean(fdr4_LN_list)) 

########## Mean detection powers for different thresholds #############
print(mean(power1_DM_list)) 
print(mean(power1_LN_list)) 
print(mean(power2_DM_list)) 
print(mean(power2_LN_list)) 
print(mean(power3_DM_list)) 
print(mean(power3_LN_list)) 
print(mean(power4_DM_list)) 
print(mean(power4_LN_list)) 

########################################### Comparing other methods ########################################################
############################################################################################################################

## The methods under comparison are MX-KF (Model-X Knockoffs) [Candes et al. (2018)], LDA-KF (vanilla LDA Knockoffs), 
## DeepLINK [Zhu et al. (2021)], and CKF (Compositional Knockoff Filter) [Srinivasan et al. (2021)].

## You can run the same simulation settings by just replacing the "Fitting the zinck model" chunk using the corresponding code blocks --

######################## Fitting MX-KF ##########################

##### For both DM & LN generation ######
set.seed(1)
kfp = knockoff.filter(X=W,y=Y,fdr = FDR,statistic = stat.glmnet_lambdasmax,offset=0) ## Change FDR accordingly!
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
W <- stat.lasso_coefdiff(W,W_tilde,Y)
T <- knockoff.threshold(W,fdr=FDR,offset=1) ## Vary FDR accordingly!
index_est <- which(W>=T) 


# For DeepLINK, please refer to the software written in Python 3.7.6 named DeepLINK publicly available in the GitHub repository : https://github.com/zifanzhu/DeepLINK
# For CKF, please refer to the software published in the supplementary material of the paper: "Compositional knockoff filter for high-dimensional regression analysis of microbiome data" 
# by Srinivasan et al. (2021)
# Note that CKF requires a Gurobi license.
              
