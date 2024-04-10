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

#### Study : RISK ####
risk_indices <- which(study_names == "RISK")
risk_otu <- combined_studies[risk_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][risk_indices]
risk_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0) ## Labelling "CD" or "UC" as 1, rest as 0

#### Study : CS-PRISM ####
prism_indices <- which(study_names == "CS-PRISM")
prism_otu <- combined_studies[prism_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][prism_indices]
prism_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0) 

#### Study : HMP2 ######

hmp_indices <- which(study_names == "HMP2")
hmp_otu <- combined_studies[hmp_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][hmp_indices]
hmp_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0) 

##### Study : MucosalIBD #####

mi_indices <- which(study_names == "MucosalIBD")
mi_otu <- combined_studies[mi_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][mi_indices]
mi_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0) 

##### Study : Pouchitis #####

pouchitis_indices <- which(study_names == "Pouchitis")
pouchitis_otu <- combined_studies[pouchitis_indices, ]
IBD_resp <- physeq_genera@sam_data[["disease"]][pouchitis_indices]
pouchitis_Y <- ifelse(IBD_resp %in% c("CD", "UC"), 1, 0) 


######### Combining all the 5 studies together #############
X <- rbind(risk_otu,prism_otu,hmp_otu,mi_otu,pouchitis_otu)
Y <- c(risk_Y,prism_Y,hmp_Y,mi_Y,pouchitis_Y)

X <- X[,order(decreasing=T,colSums(X,na.rm=T),apply(X,2L,paste,collapse=''))] ## ordering the columns w/ decreasing abundance

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

set.seed(5)
W <- stat.random_forest(X,X_tilde,Y)
T <- knockoff.threshold(W,fdr=0.1,offset = 0) ## This is the knockoff filter threshold
print(which(W>=T))

############## Feature Importance Plot #############
names <- colnames(X[,which(W>=T)]) ## Extracting the names of the important genera
data.genus <- data.frame(
  impscores = sort(W[which(W>=T)], decreasing=FALSE) , 
  name = factor(names, levels = names),
  y = seq(length(names)) * 0.9
)

plt.genus <- ggplot(data.genus) +
  geom_col(aes(impscores, name), fill = "black", width = 0.6)+theme_bw()+ylab("Genera")+xlab("Feature Statistic")

plt.genus


