############################## Codes to generate heatmaps for Figure 1(a) #####################################
###############################################################################################################
library(zinLDA)
library(zinck)
library(knockoff)
library(rstan)
library(topicmodels)
library(ggplot2)

load("count.RData") ## Loading the CRC species level data from the zinck package
dcount <- count[,order(decreasing=T,colSums(count,na.rm=T),apply(count,2L,paste,collapse=''))][,1:300] ## Ordering the columns w/ decreasing abundance

set.seed(123) 
selected_rows <- sample(1:nrow(dcount), 20)     ## Randomly select 20 subjects
selected_cols <- sample(1:ncol(dcount),30)      ## Randomly select 30 taxa
OTU <- dcount[selected_rows,selected_cols]      ## Resulting OTU matrix of dimensions 20*30
X <- OTU

draw.heatmap <- function(X, title="") {         ## Function to generate heatmaps
  reshape2::melt(asinh(X)) %>%
    dplyr::rename(sample = Var1, taxa = Var2, asinh.abun = value) %>%
    ggplot2::ggplot(., aes (x = taxa, y = sample, fill = asinh.abun)) +
    ggplot2::geom_tile() + ggplot2::theme_bw() + ggplot2::ggtitle(title) +
    ggplot2::labs(fill = "arcsinh\nabundance") +
    ggplot2::theme(plot.title = element_text(hjust = 0.5),
                   axis.title.x=element_blank(), axis.title.y=element_blank(),
                   axis.text.x = element_text(size=3, angle=90), axis.text.y = element_text(size=4)) +
   viridis::scale_fill_viridis(discrete = FALSE, direction = -1, na.value = "grey") +
    theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),
          axis.title.y = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
          panel.background = element_blank(),panel.border = element_blank())+
    theme(legend.position="none") +
    ggplot2::coord_fixed(ratio = 1)  # Fixing the aspect ratio
}
################# Zinck ########################
################################################

zinck_fit <- zinLDA::zinLDA(X, K = 9, alpha = 0.1, pi = 0.75, a = .5, b = 10)  ## Fitting the zinck model
posteriorEsts = zinLDA::posterior(zinck_fit)
theta <- posteriorEsts$theta  
beta <- posteriorEsts$beta   
full_col_names <- colnames(OTU)
### Padding zero columns so that the knockoff copy has same number of columns as X ###
new_beta <- matrix(0, nrow=9, ncol=30)
colnames(new_beta) <- full_col_names
for (col in colnames(beta)) {
  new_beta[,col] <- beta[,col]
}

X_tilde.zinck <- zinck::generateKnockoff(X,theta,new_beta,seed=1) ### Knockoff copy of X 

############# Model-X Knockoffs ##################
##################################################

X[X == 0] <- 0.5  ## Replacing zero entries with 0.5 so that log does not explode!
W <- log(X)  
W_tilde.KF <- create.second_order(W,method="equi",shrink=T) ## Generating second-order Gaussian knockoffs
X_tilde.KF <- exp(W_tilde.KF)
X_tilde.KF[X_tilde.KF<0.5] = 0

############ Vanilla LDA knockoffs ###############
##################################################

df.LDA <- as(as.matrix(X),"dgCMatrix")

vanilla.LDA <- LDA(df.LDA,k=9,method="VEM") ## Training the vanilla LDA model
theta.LDA <- vanilla.LDA@gamma
beta.LDA <- vanilla.LDA@beta
beta.LDA <- t(apply(beta.LDA, 1, function(row) row/sum(row))) ## Normalizing the betas

X_tilde.LDA <- zinck::generateKnockoff(X,theta.LDA,beta.LDA,seed=1) ## Generating vanilla LDA knockoffs

# Calculate the sparsity of each column for the Original OTU matrix
sparsity <- apply(X, 2, function(col) 1 - mean(col > 0))
# Order the matrix by decreasing sparsity
X <- X[, order(sparsity, decreasing = FALSE)]
# Draw the heatmap for the original OTU matrix
draw.heatmap(X)

## Similarly do this for zinck, Model-X and vanilla LDA knockoffs ##
sparsity <- apply(X_tilde.zinck, 2, function(col) 1 - mean(col > 0))
X_tilde.zinck <- X_tilde.zinck[, order(sparsity, decreasing = FALSE)]
draw.heatmap(X_tilde.zinck)

sparsity <- apply(X_tilde.KF, 2, function(col) 1 - mean(col > 0))
X_tilde.KF <- X_tilde.KF[, order(sparsity, decreasing = FALSE)]
draw.heatmap(X_tilde.KF)

sparsity <- apply(X_tilde.LDA, 2, function(col) 1 - mean(col > 0))
X_tilde.LDA <- X_tilde.LDA[, order(sparsity, decreasing = FALSE)]
draw.heatmap(X_tilde.LDA)








