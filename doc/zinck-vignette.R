## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE,message=FALSE)
#knitr::opts_knit$set(root.dir = "/Users/Patron/Documents/zinck research")
library(zinck)
library(knockoff)
library(dplyr)
library(ggplot2)
library(reshape2)
library(glmnet)
library(gridExtra)
library(randomForest)
library(cowplot)
library(kosel)

## ----toy-matrix,message=FALSE, warning=FALSE,results='hide'-------------------
load("/Users/Patron/Documents/zinck research/count.RData") 
load("/Users/Patron/Documents/zinck research/meta.RData") 

norm_count <- count/rowSums(count)
col_means <- colMeans(norm_count > 0)
indices <- which(col_means > 0.2)
sorted_indices <- indices[order(col_means[indices], decreasing=TRUE)]
dcount <- count[,sorted_indices][,1:400]

set.seed(123) 
selected_rows <- sample(1:nrow(dcount), 20)     ## Randomly select 20 subjects
selected_cols <- sample(1:ncol(dcount),30)      ## Randomly select 30 taxa
X <- dcount[selected_rows,selected_cols]      ## Resulting OTU matrix of dimensions 20*30


## ----knockoff-toy-matrix, message=FALSE, warning=FALSE, results='hide'--------
model_zinck <- fit.zinck(X, num_clusters=13, method="ADVI", seed=2, 
                         boundary_correction = TRUE,prior_ZIGD = TRUE)
Theta <- model_zinck$theta
Beta <- model_zinck$beta
X_zinck <- generateKnockoff(X,Theta,Beta,seed=2)


## ----heatmap-comparison, fig.width=6, fig.height=4----------------------------

heat1 <- draw_heatmap(X, "Original")

rownames(X_zinck) <- rownames(X)
heat2 <- draw_heatmap(X_zinck, "Knockoff")

plot_grid(heat1, heat2, ncol = 2, align="v")


## ----CRC-example--------------------------------------------------------------
generate_data <- function(p,seed){
  dcount <- count[,order(decreasing=T,colSums(count,na.rm=T),apply(count,2L,paste,collapse=''))]
  ## ordering the columns w/ decreasing abundance
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

X <- generate_data(p=ntaxa, seed=1)$X
Y <- generate_data(p=ntaxa, seed=1)$Y


## ----learn-parameters, message=FALSE, warning=FALSE, echo=TRUE, results='hide'----
fit <- fit.zinck(X,num_clusters = 15,method="ADVI",seed=12,elbo_samples = 100)
theta <- fit$theta
beta <- fit$beta


## ----knockoffgen--------------------------------------------------------------
X_tilde <- zinck::generateKnockoff(X,theta,beta,seed=1) ## getting the knockoff copy

## ----rf-model-----------------------------------------------------------------
cts_rf <- suppressWarnings(zinck.filter(X,X_tilde,Y,model="Random Forest",fdr=0.2,offset=0,mtry=200,
                                        seed=12,metric = "Accuracy", rftuning = TRUE))
index_est <- cts_rf[["selected"]]


## ----evaluation---------------------------------------------------------------
index <- 1:30
### Evaluation metrics ###
FN <- sum(index %in% index_est == FALSE) ## False Negatives
FP <- sum(index_est %in% index == FALSE) ## False Positives
TP <- sum(index_est %in% index == TRUE) ## True Positives

# ## False Discovery Proportion ##
estimated_FDR <- FP/(FP+TP)
print(estimated_FDR)

## Power ##
estimated_power <- TP/(TP+FN)
print(estimated_power)

## ----CRCspecies---------------------------------------------------------------
norm_count <- count/rowSums(count)
col_means <- colMeans(norm_count>0)
indices <- which(col_means > 0.2)
sorted_indices <- indices[order(col_means[indices],decreasing = TRUE)]
dcount <- count[,sorted_indices]

X <- dcount
Y <- ifelse(meta$Group=="CRC",1,0)

## ----CRCknockoff, message=FALSE, warning=FALSE, echo=TRUE, results='hide'-----
fitCRC <- fit.zinck(X,num_clusters=20,method="ADVI",seed=1,boundary_correction = TRUE,
                    elbo_samples = 100,importance_resampling = FALSE,
                    prior_ZIGD = TRUE)
theta <- fitCRC$theta
beta <- fitCRC$beta
X_tilde <- zinck::generateKnockoff(X,theta,beta,seed=1) ## Generating the knockoff copy


## ----random-forest------------------------------------------------------------
filter_zinck <- zinck.filter(as.matrix(X),as.matrix(X_tilde),as.factor(Y),
                model="Random Forest",offset = 1,seed=1,mtry=28,
                rftuning=TRUE,metric = "Gini")

selected_species <- filter_zinck$selected


## Importance scores ##
W <- filter_zinck$W

## Threshold ##
T <- filter_zinck$T

names <- colnames(X[,which(W>=T)])

### Creating the data frame with Feature Importance Statistics
data.species <- data.frame(             
  impscores = sort(W[which(W>=T)], decreasing=FALSE) , 
  name = factor(names, levels = names),
  y = seq(length(names)) * 0.9
)


## ----feature-statistics,  fig.width=6, fig.height=4---------------------------
data.species <- data.frame(impscores=sort(W[which(W>=T)],decreasing = FALSE), 
                         name=factor(names, levels=names), y=seq(length(names))*0.9)

plot.species <- ggplot(data.species)+geom_col(aes(impscores,name),
                                          fill="black",width=0.6)+theme_bw()+
  ylab("Species")+xlab("Feature Statistic")

plot.species


## ----optimal-k, warning=FALSE, message=FALSE, echo=TRUE, fig.width=6, fig.height=4----

load("/Users/Patron/Documents/zinck research/meta.RData") 
dcount <- count.genus[,order(decreasing=T,colSums(count.genus,na.rm=T), 
                             apply(count.genus,2L,paste,collapse=''))] 
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


