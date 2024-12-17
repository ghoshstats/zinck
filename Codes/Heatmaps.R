############################## Codes to generate heatmaps for Figure 1(a) #####################################
###############################################################################################################
library(zinck)
library(knockoff)
library(rstan)
library(topicmodels)
library(ggplot2)

load("count.RData") ## Loading the CRC species level data from the zinck package
norm_count <- count/rowSums(count)
col_means <- colMeans(norm_count > 0)
indices <- which(col_means > 0.2)
sorted_indices <- indices[order(col_means[indices], decreasing=TRUE)]
dcount <- count[,sorted_indices][,1:400]

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
dlt <- rep(0,ncol(X)) ## Initializing deltas with the individual column sparsities
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
  K = 13,
  V = ncol(X),
  D = nrow(X),
  n = X,
  delta = dlt
)

stan.model = stan_model(model_code = zinck_code)
set.seed(2) ## Set seeds carefully since vb is sensitive to starting points. If there is an error for iteration i switch to seed = 11 ##
fit <- vb(stan.model, data=zinLDA_stan_data, algorithm="meanfield", iter=10000)
theta <- fit@sim[["est"]][["theta"]]
beta <- fit@sim[["est"]][["beta"]]

X_tilde.zinck <- zinck::generateKnockoff(X,theta,beta,seed=2) ### Knockoff copy of X 

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








