## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE,message=FALSE)
knitr::opts_knit$set(root.dir = "/Users/Patron/Documents/zinLDA research")
library(zinck)
library(knockoff)
library(dplyr)
library(zinLDA)
library(ggplot2)
library(reshape2)
library(glmnet)
library(gridExtra)
library(randomForest)
library(cowplot)
library(kosel)

## ----toy-matrix, message=FALSE------------------------------------------------
set.seed(1) ##reproducibility

N.d=zinLDA::rdu(n=20,min=400,max=500) # Drawing random library sizes between 400, 500

sim_data = zinLDA::simulateZINLDA(D=20,V=30,N=N.d,K=5,Alpha=0.1,Pi=0.8,a=0.5,b=10)

X_original <- sim_data$sampleTaxaMatrix ## The original sample taxa count matrix


## ----knockoff-toy-matrix, message=FALSE, warning=FALSE, results='hide'--------
model_zinck <- fit.zinck(X_original, num_clusters=5, method="Gibbs", seed=1)

Theta <- model_zinck$theta

Beta <- model_zinck$beta

X_zinck <- generateKnockoff(X_original,Theta,Beta,seed=1)

rownames(X_zinck) <- rownames(X_original)


## ----heatmap-comparison, fig.width=6, fig.height=4----------------------------

heat1 <- draw_heatmap(X_original, "Original")

heat2 <- draw_heatmap(X_zinck, "Knockoff")

plot_grid(heat1, heat2, ncol = 2, align="v")


## ----CRC-example--------------------------------------------------------------
load("count.Rdata")

dcount <- count[,order(decreasing=T,colSums(count,na.rm=T),
                       apply(count,2L,paste,collapse=''))][,1:300] 
## ordering the columns w/ decreasing abundance

## Random Subject Selection
set.seed(1)

sel_index <- rbinom(nrow(dcount),size=1,prob=0.5)

selected_samples <- which(sel_index==1)

X <- dcount[selected_samples,]


## ----signal-injection---------------------------------------------------------
Five1 <- c(-3,3,2.5,-1,-1.5)

Five2 <- c(3,3,-2,-2,-2)

Five3 <- c(1,-1,3,-2,-1)

Five4 <- c(-1,1,2,-1,-1)

Five5 <- c(3,3,-3,-2,-1)

Five6 <- c(-1,1,-2,1,1)

Five_all <- c(Five1,Five2,Five3,Five4,Five5,Five6)

randBeta <- rep(0,300)

set.seed(1)

rand_indices <- sample(1:200,size=30,replace=FALSE)

set.seed(1)

randBeta[rand_indices] <- sample(Five_all, size=30, replace=FALSE)


## ----response-generation------------------------------------------------------
n = nrow(X)

W <- log_normalize(X)

set.seed(1)

eps=rnorm(n,mean = 0, sd=1)

Y <- W %*% randBeta + eps


## ----learn-parameters, message=FALSE, warning=FALSE, echo=TRUE, results='hide'----
species_fit <- fit.zinck(X,num_clusters = 6,method = "ADVI", seed=123, alpha_param = 0.1)

theta <- species_fit$theta

beta <- species_fit$beta


## ----knockoffgen--------------------------------------------------------------
X_tilde <- generateKnockoff(X,theta,beta,seed=1)

## ----glmnet-model-------------------------------------------------------------
W_tilde <- log_normalize(X_tilde)

selected_species = zinck.filter(W,W_tilde,Y,model="glmnet",fdr=0.1,offset=1)


## ----evaluation---------------------------------------------------------------
index <- rand_indices

index_est <- selected_species

neg_index <- (1:300)[-c(index)]

if(length(index_est)==0){
  neg_index_est <- 1:300
} else{
  neg_index_est <- (1:300)[-c(index_est)]
}

# ### Evaluation metrics ###
TP <- sum(index_est %in% index==TRUE) ## True Positives
FP <- sum(index_est %in% index==FALSE) ## False Positives
FN <- sum(index %in% index_est == FALSE) ## False Negatives
TN <- sum(neg_index_est %in% neg_index == TRUE) ## True Negatives

# ## False Discovery Proportion ##
estimated_FDR <- FP/(FP+TP)
print(estimated_FDR)

## Power ##
estimated_power <- TP/(TP+FN)
print(estimated_power)

## ----CRCgenus-----------------------------------------------------------------
load("count.genus.RData")
load("meta.RData")

X <- count.genus[,order(decreasing = T,colSums(count.genus,na.rm=T),apply(count.genus,2L,paste,collapse=''))]

# ## Case-control status (Response) ##

Y_cat <- as.factor(meta$Group)
lookup <- c("CTR"=0,"CRC"=1)
Y <- lookup[Y_cat]


## ----CRCknockoff, message=FALSE, warning=FALSE, echo=TRUE, results='hide'-----
fit_CRC <- fit.zinck(X,num_clusters = 8, method="ADVI",tuned=FALSE,seed=11, alpha_param = 0.1)

theta_CRC <- fit_CRC$theta

beta_CRC <- fit_CRC$beta

X_tilde_CRC <- generateKnockoff(X,theta_CRC,beta_CRC,seed=1)


## ----random-forest------------------------------------------------------------
selected_genera = zinck.filter(X,X_tilde_CRC,Y,model="Random Forest",fdr=0.1,offset=1,seed=11)

## ----feature-statistics,  fig.width=6, fig.height=4---------------------------
set.seed(1)
W <- stat.random_forest(X,X_tilde_CRC,Y)
T <- knockoff.threshold(W,fdr=0.1,offset=1)

names <- colnames(X[,which(W>=T)])
data.genus <- data.frame(impscores=sort(W[which(W>=T)],decreasing = FALSE), 
                         name=factor(names, levels=names), y=seq(length(names))*0.9)

plot.genus <- ggplot(data.genus)+geom_col(aes(impscores,name),
                                          fill="black",width=0.6)+theme_bw()+
  ylab("Genera")+xlab("Feature Statistic")

plot.genus


## ----optimal-k, warning=FALSE, message=FALSE, echo=TRUE, fig.width=6, fig.height=4----

load("count.genus.RData")
dcount <- count.genus[,order(decreasing=T,colSums(count.genus,na.rm=T), apply(count.genus,2L,paste,collapse=''))] 
## ordering the columns w/ decreasing abundance
X <- dcount
kmin =8 
kmax=17
K_values <- seq(kmin, kmax, by=1 )


## Run this command to find the optimal K value ##
#optimal_k(X,kmin=8,kmax=17,seed_list=list(1,11,1,1,1,1,1,1,1,1))

## On running optimal_k(), the Jensen-Shannon divergence values 
js_values <- c(0.266242705918838,0.373342570843659,0.36512930638952,0.375253322971532,0.430181098394843,0.435903467160891,0.444518304173072,0.451035137033785,0.515442875511868,0.442322087559003)

data <- data.frame(K_values, js_values)
p <- ggplot(data, aes(x = K_values, y = js_values)) +
    geom_line() +
    geom_point(size = 4) +
    labs(
      title = "Optimal Number of clusters based on Jensen-Shannon Divergence",
      x = "Number of clusters (K)",
      y = "Average JS Divergence"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "grey"),
      panel.grid.minor = element_line(size = 0.5, linetype = 'solid', color = "grey")
    )
print(p)


