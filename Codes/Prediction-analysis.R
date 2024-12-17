################################## Prediction Analysis for CRC & IBD data ##########################################
####################################################################################################################
library(caret)
library(pROC)
library(randomForest)
library(phyloseq)
library(dplyr)
library(reshape2)
library(pROC)
library(GUniFrac)
library(ggplot2)
library(knockoff)
library(glmnet)
library(rstan)
library(topicmodels)

load("meta.RData")
load("count.RData")
load("genera.RData")


##################################### AUROC values for CRC #############################################
norm_count <- count/rowSums(count)
col_means <- colMeans(norm_count > 0)
indices <- which(col_means > 0.2)
sorted_indices <- indices[order(col_means[indices], decreasing=TRUE)]
dcount <- count[,sorted_indices][,1:300]
X <- dcount
Y <- as.factor(meta$Group)
Y <- ifelse(meta$Group=="CRC",1,0)

############ Leave-one study out design ################
left_out_study <- "AT-CRC" ### Also other studies -- "US-CRC", "CN-CRC", "DE-CRC", "FR-CRC"
study_names <- unique(meta$Study)
remaining_studies <- setdiff(study_names, left_out_study)

roc_zinck_AT <- rep(0,100) ## vector for storing AUROC values 
nspecies_zinck_AT <- rep(0,100) ## vector for storing species level counts

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

for(i in 1:100)
{
  tryCatch({
  set.seed(i)  
  train_indices_list <- list()
  val_indices_list <- list()
  for(study in unique(meta$Study)) {
    study_indices <- which(meta$Study == study)
    ## 80:20 train and validation sets
    train_indices_temp <- createDataPartition(meta$Group[study_indices], p = 0.8, list = TRUE, times = 1)[[1]]     
    # Convert to actual indices from the subset of study_indices
    train_indices <- study_indices[train_indices_temp]
    # Determine validation indices as the set difference
    val_indices <- setdiff(study_indices, train_indices)
    # Store indices in their respective lists
    train_indices_list[[study]] <- train_indices
    val_indices_list[[study]] <- val_indices
  }
  
  training_indices <- unlist(train_indices_list[setdiff(study_names, left_out_study)])
  validation_indices <- unlist(val_indices_list[setdiff(study_names, left_out_study)])
  
  X_train <- dcount[training_indices,1:300]
  Y_train <- as.factor(meta$Group[training_indices])
  Y_train <- ifelse(meta$Group[training_indices]=="CRC",1,0)

  ########### Initializing ADVI ##############  
  dlt <- rep(0,ncol(X_train))

  for(t in (1:ncol(X_train)))
  {
    dlt[t] <- 1-mean(X_train[,t]>0)
  }
   
  zinck_stan_data <- list(
    K = 15,
    V = ncol(X_train),
    D = nrow(X_train),
    n = X_train,
    alpha = rep(0.1, 15),
    gamma1 = rep(0.5, ncol(X_train)),
    gamma2 = rep(10,ncol(X_train)),
    delta = dlt
  )

  set.seed(12) ## Set seeds carefully since vb is sensitive to starting points. If there is an error for iteration i switch to seed = 11 ##
  fitCRC_train <- vb(stan.model, data=zinck_stan_data, algorithm="meanfield", iter=10000)
  theta <- fitCRC_train@sim[["est"]][["theta"]]
  beta <- fitCRC_train@sim[["est"]][["beta"]]
  X_tilde.zinck <- zinck::generateKnockoff(X_train,theta,beta,seed=1) ## getting the knockoff copy
  X_aug <- cbind(X_train,X_tilde.zinck)
  
  ######  Random Forest ####
  bestmtry <- tuneRF(X_aug,as.factor(Y_train),stepFactor = 1.5,improve=1e-5,ntree=1000,plot=FALSE) ## tuning mtry parameter
  m <- bestmtry[as.numeric(which.min(bestmtry[,"OOBError"])),1]
  df_X <- as.data.frame(X_aug)
  colnames(df_X) <- NULL
  rownames(df_X) <- NULL
  df_X$Y_train <- Y_train
  model_rf <- randomForest(formula=as.factor(Y_train)~.,ntree=1000,mtry=m,importance=TRUE,data=as.matrix(df_X))
  
  cf <- as.data.frame(importance(model_rf))[,3]  ## Gini impurities 
  #cf <- as.data.frame(importance(model_rf))[,1]  ## % Node impurities
     
  W <- abs(cf[1:300])-abs(cf[301:600]) ## feature statistics
  T <- knockoff.threshold(W,fdr = 0.1, offset = 0) ## knockoff threshold
  selected_taxa <- which(W>=T) 
  
  all_fprs <- list()
  all_tprs <- list()
  X_val <- dcount[validation_indices,1:300]
  Y_val <- as.factor(meta$Group[validation_indices])
  Y_val <- ifelse(meta$Group[validation_indices]=="CRC",1,0)

  ########## Normalizing the training and validation sets ############
  norm_train <- as.data.frame(t(apply(X_train,1, function(x) x/sum(x))))
  norm_val <- as.data.frame(t(apply(X_val,1, function(x) x/sum(x))))
  auroc.all <- c()
  df <- as.data.frame(norm_train[,selected_taxa])
  df1 <- as.data.frame(norm_val[,selected_taxa])
  
  #### Train and validation data to tune the random forest ####
  xtrain = df
  ytrain = Y_train
  xval = df1
  yval = Y_val
  
  ######### Tuning the hyperparameters #############
  ## ntrees range from 500 to 1000.
  ## mytrys have 3 options.
  ntrees <- c(50:100)*10
  #ntrees <- c(500,1000)
  mtrys <- c(round(sqrt(5)) / 2, round(sqrt(5)), round(sqrt(5)) * 2)
  
  aucall <- matrix(0,3,51)
  for(l in 1:length(ntrees)){
    for(ll in 1:length(mtrys)){
      set.seed(2022)
      output.forest <- randomForest::randomForest(x = xtrain, y = as.factor(ytrain),
                                                  xtest = xval, ytest = as.factor(yval),
                                                  ntree = ntrees[l], mytry = mtrys[ll] ,
                                                  importance= TRUE, keep.forest=TRUE)
      
      prob <- predict(output.forest, newdata = norm_val, type = "prob")
      roc_object <- roc(response = Y_val, predictor = prob[,"1"], verbose = 0)
      aucall[ll,l] <- roc_object$auc
    }
  }
  ######### plug-in the best hyperparameters #############
  
  best_ntree <- ntrees[which.max(apply(aucall, 2, max))]  ## ntree
  best_mtry <- mtrys[which.max(apply(aucall, 1, max))]  ## mtry
  # Training the final model on all remaining data with best hyperparameters
  all_training_indices <- which(meta$Study != left_out_study)
  X_all_train <- dcount[all_training_indices,1:300]
  norm_all_train <- as.data.frame(t(apply(X_all_train,1, function(x) x/sum(x))))
  Y_all_train <- as.factor(meta$Group[all_training_indices])
  Y_all_train <- ifelse(meta$Group[all_training_indices]=="CRC",1,0)
  
  df_all_train <- as.data.frame(norm_all_train[,selected_taxa])
  
  #### Train and test data ####
  xtrain = df_all_train
  ytrain = Y_all_train
  
  test_indices <- which(meta$Study == left_out_study)
  X_test <- dcount[test_indices,1:300]
  norm_test <- as.data.frame(t(apply(X_test,1, function(x) x/sum(x))))
  Y_test <- as.factor(meta$Group[test_indices])
  Y_test <- ifelse(meta$Group[test_indices]=="CRC",1,0)
  
  df_test <- as.data.frame(norm_test[,selected_taxa])
  xtest <- df_test
  ytest <- Y_test
  
  set.seed(2022)
  output.forest <- randomForest(x = xtrain, y = as.factor(ytrain), 
                                xtest = xtest, ytest = as.factor(ytest),
                                ntree = best_ntree, mytry = best_mtry,
                                importance= TRUE, keep.forest=TRUE)
  
  y_pred <- predict(output.forest, newdata = norm_test, type = "prob")
  
  roc_object <- roc(response = Y_test, predictor = y_pred[,"1"])
  
  # Store FPR, TPR, and Thresholds
  all_fprs[[left_out_study]] <- (1-roc_object$specificities)
  all_tprs[[left_out_study]] <- roc_object$sensitivities
  common_fprs <- seq(0, 1, by=0.01)
  mean_tprs <- sapply(common_fprs, function(fpr) {
  interpolated_tprs <- approx(x = all_fprs[[1]], y = all_tprs[[1]], xout = fpr)$y
  mean(interpolated_tprs, na.rm = TRUE)
  })
  roc_df <- data.frame(fpr = common_fprs, tpr = mean_tprs)
  roc_df <- roc_df[order(roc_df$fpr), ]
  
  # Calculating AUC  
  roc_zinck_AT[i]  <- sum(with(roc_df, diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2))
  nspecies_zinck_AT[i] <- length(selected_taxa)
  }, error = function(e) {
    cat("An error occurred in iteration", i, "\n")
  })
}

print(roc_zinck_AT)
print(nspecies_zinck_AT)

##############################################################################################################
###################################### AUROC values for IBD data #############################################

combined_studies <- as.data.frame(t(physeq_genera@otu_table))
study_names <- physeq_genera@sam_data[["dataset_name"]]

## RISK ##
risk_indices <- which(study_names == "RISK")
# Subset combined_studies using these indices
risk_otu <- combined_studies[risk_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][risk_indices]
risk_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0)

## CS-PRISM ##
prism_indices <- which(study_names == "CS-PRISM")
# Subset combined_studies using these indices
prism_otu <- combined_studies[prism_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][prism_indices]
prism_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0)

## HMP2 ##
HMP_indices <- which(study_names == "HMP2")
# Subset combined_studies using these indices
hmp_otu <- combined_studies[HMP_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][HMP_indices]
hmp_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0)

## Pouchitis ##
pouchitis_indices <- which(study_names == "Pouchitis")
# Subset combined_studies using these indices
pouchitis_otu <- combined_studies[pouchitis_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][pouchitis_indices]
pouchitis_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0)

## MucosalIBD ##
mi_indices <- which(study_names == "MucosalIBD")
# Subset combined_studies using these indices
mi_otu <- combined_studies[mi_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][mi_indices]
mi_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0)

metadata <- physeq_genera@sam_data

# Get the indices of the last occurrence of each unique subject_accession
latest_indices <- sapply(unique(metadata$subject_accession), function(x) {
  max(which(metadata$subject_accession == x))
})

# Subset the metadata to keep only the latest samples
latest_metadata <- metadata[latest_indices, ]

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

risk_Y_latest <- risk_Y[row.names(risk_otu) %in% risk_row_names]
prism_Y_latest <- prism_Y[row.names(prism_otu) %in% prism_row_names]
hmp_Y_latest <- hmp_Y[row.names(hmp_otu) %in% hmp_row_names]
mi_Y_latest <- mi_Y[row.names(mi_otu) %in% mi_row_names]
pouchitis_Y_latest <- pouchitis_Y[row.names(pouchitis_otu) %in% pouchitis_row_names]

combined_otu <- rbind(risk_otu_latest,prism_otu_latest,hmp_otu_latest,mi_otu_latest,pouchitis_otu_latest)
dcount <- combined_otu[,order(decreasing=T,colSums(combined_otu,na.rm=T),apply(combined_otu,2L,paste,collapse=''))] ## ordering the columns w/ decreasing abundance

Y <- c(prism_Y_latest, hmp_Y_latest, mi_Y_latest, pouchitis_Y_latest,risk_Y_latest)

study_names <- c("RISK","CS-PRISM","HMP2","MucosalIBD","Pouchitis")
meta_IBD <- latest_metadata[latest_metadata$dataset_name %in% study_names,]

IBD_resp <- meta_IBD$disease
meta_IBD$disease <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0)

ngenera_zinck_risk <- rep(0,100) ## detected number of taxa for each left out study 
roc_zinck_risk <- rep(0,100) ## vector of AUROC values

###################### The zinck model ######################
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
stan.model = stan_model(model_code = zinck_code) ## Defining the stan model
                               
############## Leave - one - study - out analysis ############
                               
left_out_study <- "RISK"  ## Try with other studies like "CS-PRISM", "HMP2", "MucosalIBD", and "Pouchitis"
# Determine the indices for the remaining studies
remaining_studies <- setdiff(study_names, left_out_study)

for(i in 1:100)
{
tryCatch({
set.seed(i)  
train_indices_list <- list()
val_indices_list <- list()
for(study in unique(meta_IBD$dataset_name)) {
  study_indices <- which(meta_IBD$dataset_name == study)
  
  train_indices_temp <- createDataPartition(meta_IBD$disease[study_indices], p = 0.8, list = TRUE, times = 1)[[1]] ## 80:20 train, validation split
  
  # Convert to actual indices from the subset of study_indices
  train_indices <- study_indices[train_indices_temp]
  # Determine validation indices as the set difference
  val_indices <- setdiff(study_indices, train_indices)
  # Store indices in their respective lists
  train_indices_list[[study]] <- train_indices
  val_indices_list[[study]] <- val_indices
}

training_indices <- unlist(train_indices_list[setdiff(study_names, left_out_study)])
validation_indices <- unlist(val_indices_list[setdiff(study_names, left_out_study)])

# Nested cross-validation for hyperparameter tuning
# Assuming a simple approach: iterate over splits for validation and use the rest for training
# Map indices to row names from `meta_IBD`
training_row_names <- rownames(meta_IBD)[training_indices]
validation_row_names <- rownames(meta_IBD)[validation_indices]

X_train <- dcount[training_row_names, ]
Y_train <- meta_IBD$disease[training_indices]

dlt <- rep(0,ncol(X_train))

for(t in (1:ncol(X_train)))
{
  dlt[t] <- 1-mean(X_train[,t]>0)
  if(dlt[t]==0)
  {
    dlt[t] = dlt[t]+0.01
  }
  if (dlt[t]==1)
  {
    dlt[t] = dlt[t]-0.01
  }

}
# 
zinck_stan_data <- list(
  K = 18,
  V = ncol(X_train),
  D = nrow(X_train),
  n = X_train,
  delta = dlt
)

set.seed(1) ## Very sensitive to initialization
fitIBD_train <- vb(stan.model, data=zinck_stan_data, algorithm="meanfield", iter=10000)

theta <- fitIBD_train@sim[["est"]][["theta"]]
beta <- fitIBD_train@sim[["est"]][["beta"]]
X_tilde_IBD <- zinck::generateKnockoff(X_train,theta,beta,seed=1)

set.seed(1) 
W <- stat.random_forest(X_train,X_tilde_IBD,as.factor(Y_train)) ## Fitting a random forest model
T <- knockoff.threshold(W,fdr=0.1,offset=0) ## knockoff threshold
print(which(W>=T))

selected_taxa <- which(W>=T)

all_fprs <- list()
all_tprs <- list()

X_val <- dcount[validation_row_names, ]
Y_val <- meta_IBD$disease[validation_indices]

######### Normalizing the training and validation sets ############
norm_train <- as.data.frame(t(apply(X_train,1, function(x) x/sum(x))))
norm_val <- as.data.frame(t(apply(X_val,1, function(x) x/sum(x))))

auroc.all <- c()
df <- as.data.frame(norm_train[,selected_taxa])
df1 <- as.data.frame(norm_val[,selected_taxa])

xtrain = df
ytrain = Y_train
xval = df1
yval = Y_val

######### Tuning the hyperparameters #############
## ntrees range from 500 to 1000.
## mytrys have 3 options.
ntrees <- c(50:100)*10
#ntrees <- c(500,1000)
mtrys <- c(round(sqrt(5)) / 2, round(sqrt(5)), round(sqrt(5)) * 2)

aucall <- matrix(0,3,51)
for(l in 1:length(ntrees)){
  for(ll in 1:length(mtrys)){
    set.seed(2022)
    output.forest <- randomForest::randomForest(x = xtrain, y = as.factor(ytrain),
                                                xtest = xval, ytest = as.factor(yval),
                                                ntree = ntrees[l], mytry = mtrys[ll] ,
                                                importance= TRUE, keep.forest=TRUE)
    
    prob <- predict(output.forest, newdata = norm_val, type = "prob")
    roc_object <- roc(response = Y_val, predictor = prob[,"1"], verbose = 0)
    aucall[ll,l] <- roc_object$auc
  }
}
######### plug-in the best hyperparameters #############

best_ntree <- ntrees[which.max(apply(aucall, 2, max))]  ## ntree
best_mtry <- mtrys[which.max(apply(aucall, 1, max))]  ## mtry
# Training the final model on all remaining data with best hyperparameters
all_training_indices <- which(meta_IBD$dataset_name != left_out_study)
training_row_names_all <- rownames(meta_IBD)[all_training_indices]
X_all_train <- dcount[training_row_names_all, ]

norm_all_train <- as.data.frame(t(apply(X_all_train,1, function(x) x/sum(x))))
Y_all_train <- as.factor(meta_IBD$disease[all_training_indices])

df_all_train <- as.data.frame(norm_all_train[,selected_taxa])

#### Train and test data ####
xtrain = df_all_train
ytrain = Y_all_train

test_indices <- which(meta_IBD$dataset_name == left_out_study)
test_row_names <- rownames(meta_IBD)[test_indices]
X_test <- dcount[test_row_names, ]
norm_test <- as.data.frame(t(apply(X_test,1, function(x) x/sum(x))))
Y_test <- as.factor(meta_IBD$disease[test_indices])

df_test <- as.data.frame(norm_test[,selected_taxa])
xtest <- df_test
ytest <- Y_test

set.seed(2022)
output.forest <- randomForest(x = xtrain, y = as.factor(ytrain), 
                              xtest = xtest, ytest = as.factor(ytest),
                              ntree = best_ntree, mytry = best_mtry,
                              importance= TRUE, keep.forest=TRUE)

y_pred <- predict(output.forest, newdata = norm_test, type = "prob")

roc_object <- roc(response = Y_test, predictor = y_pred[,"1"])

# Store FPR, TPR, and Thresholds
all_fprs[[left_out_study]] <- (1-roc_object$specificities)
all_tprs[[left_out_study]] <- roc_object$sensitivities

common_fprs <- seq(0, 1, by=0.01)

mean_tprs <- sapply(common_fprs, function(fpr) {
  interpolated_tprs <- approx(x = all_fprs[[1]], y = all_tprs[[1]], xout = fpr)$y
  mean(interpolated_tprs, na.rm = TRUE)
})

roc_df <- data.frame(fpr = common_fprs, tpr = mean_tprs)
roc_df <- roc_df[order(roc_df$fpr), ]
# Calculating AUC
roc_zinck_risk[i] <- sum(with(roc_df, diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2))
ngenera_zinck_risk[i] <- length(selected_taxa)
}, error = function(e) {
    cat("An error occurred in iteration", i, "\n")
  })
}

                               

                               
