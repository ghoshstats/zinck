#' fit.zinck
#'
#' Fit the Zinck model to the data using either ADVI or Gibbs sampling methods.
#'
#' @param X An OTU matrix with dimensions \eqn{D \times p}.
#' @param num_clusters An integer specifying the number of clusters.
#' @param method A character string, either "ADVI" or "Gibbs", specifying the method to fit the model.
#' @param seed An integer used to set the seed for reproducibility.
#' @param init_values A list of initial values for the ADVI algorithm. This parameter is optional and should be used only with the ADVI method. If NULL, the algorithm uses default initialization.
#' @param boundary_correction A logical value. If TRUE, it adds and subtracts small numbers so that the log-likelihood doesn't blow up at the boundary points 0 or 1, default = FALSE.
#' @param prior_ZIGD A logical value. If TRUE and method is "ADVI", the model will impose Gamma priors on the ZIGD hyperparameters.
#' @param alpha_param A positive real. The symmetric smoothed Dirichlet parameter for the cluster distributions, default=0.1
#' @param elbo_samples A positive integer. The number of samples for Monte Carlo estimate of ELBO (objective function), defaulting to 100. (ELBO stands for "the evidence lower bound".)name description
#' @param importance_resampling	 Logical scalar (defaulting to FALSE) indicating whether to do importance resampling to adjust the draws at the optimum to be more like draws from the posterior distribution.
#' 
#' @references
#' Kucukelbir, A., Tran, D., Ranganath, R., Gelman, A., and Blei, D. (2017).
#' Automatic Differentiation Variational Inference.
#' Journal of Machine Learning Research 18(14).
#'
#' Deek, R., and Li, H. (2021).
#' A Zero-Inflated Latent Dirichlet Allocation Model for Microbiome Studies.
#' Frontiers in Genetics 11.
#'
#' @return A list containing the posterior estimates of beta and theta.
#' @export
#'
#' @examples
#' # Generate a random OTU matrix
#' X <- matrix(rpois(20, lambda = 5), nrow = 5)
#'
#' # Fit the Zinck model using ADVI with default initialization
#' fit <- fit.zinck(X, num_clusters = 3, method = "ADVI", alpha_param=1)
#'
#' # Fit the Zinck model using ADVI with custom initial values
#' init_vals <- list(theta = array(runif(15), dim = c(5, 3)),
#'                   zeta = array(runif(12), dim = c(3, 4)))
#' fit_with_init <- fit.zinck(X, num_clusters = 3, method = "ADVI", init_values = init_vals, alpha_param=1)
#'


fit.zinck <- function(X, num_clusters, method = c("ADVI", "Gibbs"), seed=NULL, init_values=NULL, alpha_param = 0.1, boundary_correction = FALSE,importance_resampling=FALSE,elbo_samples=500, prior_ZIGD = FALSE) {

  # Check if X is a matrix
  if (!is.matrix(X)) {
    stop("Error: X should be a matrix.")
  }

  # Check if X is a count matrix: all elements should be non-negative integers
  if (!all(X >= 0 & X == floor(X))) {
    stop("Error: X should be a count matrix (all elements must be non-negative integers).")
  }

  # Check other input errors
  if (!is.numeric(num_clusters) || length(num_clusters) != 1 || num_clusters <= 0) {
    stop("Error: num_clusters should be a positive integer.")
  }

  method <- match.arg(method)


  ### Initialize Deltas for zinck ###
  if(boundary_correction == FALSE){
  dlt <- 1 - colMeans(X > 0)
  } else {
  dlt <- 1 - colMeans(X > 0)
  # Adding and subtracting small numbers so that the log-likelihood doesn't blow up at the boundary points [0,1]
  dlt[dlt == 0] <- dlt[dlt == 0] + 0.01
  dlt[dlt == 1] <- dlt[dlt == 1] - 0.01
  }
  if (!is.null(init_values)) {
    if (any(init_values[["theta"]] < 0)) {
      stop("Invalid initial values: Elements of 'theta' must be non-negative.")
    }
    init_values[["theta"]] <- sweep(init_values[["theta"]],1,rowSums(init_values[["theta"]]),FUN="/")
  }

  if (method == "ADVI" && prior_ZIGD == TRUE) {
  zinck_init_data <- list(
      K = num_clusters,
      V = ncol(X),
      D = nrow(X),
      n = X,
      delta = dlt
    )

    zinck_init_code <- "data {
    int<lower=1> K; // num clusters
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
    init_stan.model = stan_model(model_code = zinck_init_code)
    ####### Initial ADVI fit #########
    # Use init_values if provided
    if (!is.null(init_values)) {
      set.seed(seed)
      fit_initial <- suppressWarnings(vb(init_stan.model, data=zinck_init_data,init=init_values, algorithm="meanfield", iter=10000,importance_resampling=importance_resampling,elbo_samples=elbo_samples ))
    } else {
      set.seed(seed)
      fit_initial <- suppressWarnings(vb(init_stan.model, data=zinck_init_data, algorithm="meanfield", iter=10000,importance_resampling=importance_resampling,elbo_samples=elbo_samples))
      theta <- fit_initial@sim[["est"]][["theta"]]
      beta <- fit_initial@sim[["est"]][["beta"]]
      out = list(beta = beta,
                 theta = theta)
      return(out)
}} else if (method == "ADVI" && prior_ZIGD == FALSE) {
    zinck_default_data <- list(
      K = num_clusters,
      V = ncol(X),
      D = nrow(X),
      n = X,
      alpha = rep(alpha_param, num_clusters),
      gamma1 = rep(0.5, ncol(X)),
      gamma2 = rep(10,ncol(X)),
      delta = dlt
    )
    default_zinck <- "data {
  int<lower=1> K; // num clusters
  int<lower=1> V; // num taxa
  int<lower=0> D; // num subjects
  int<lower=0> n[D, V]; // OTU matrix

  // hyperparameters
  vector<lower=0>[K] alpha;
  vector<lower=0>[V] gamma1;
  vector<lower=0>[V] gamma2;
  vector<lower=0, upper=1>[V] delta;
}
parameters {
  simplex[K] theta[D]; // cluster mixtures
  vector<lower=0,upper=1>[V] zeta[K]; // zero-inflated betas
}


transformed parameters {
  vector<lower=0>[V] beta[K];
  for (k in 1:K) {
	beta[k,1] =  zeta[k,1];
  for (m in 2:V) {
    beta[k,m] = zeta[k,m]*prod(1 - zeta[k,1:(m - 1)]);  // stick breaking process
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
default_stan.model = stan_model(model_code = default_zinck)
if (!is.null(init_values)) {
  set.seed(seed)
  fit_default <- suppressWarnings(vb(default_stan.model, data=zinck_default_data,init=init_values, algorithm="meanfield", iter=10000,importance_resampling=importance_resampling,elbo_samples=elbo_samples))
} else {
  set.seed(seed)
  fit_default <- suppressWarnings(vb(default_stan.model, data=zinck_default_data, algorithm="meanfield", iter=10000,importance_resampling=importance_resampling,elbo_samples=elbo_samples))
}
###### Posterior Estimates of theta and beta ########

theta <- fit_default@sim[["est"]][["theta"]]
beta <- fit_default@sim[["est"]][["beta"]]
out = list(beta = beta,
           theta = theta)
return(out)
  } else if (method == "Gibbs") {
    if (tuned) warning("Tuning is not applicable for the Gibbs method. Argument 'tuned' will be ignored.")
    fit_Gibbs = zinLDA::zinLDA(X, K = num_clusters, alpha = 0.1, pi = (1-mean(X>0)), a = 0.5, b = 10)
    posteriorEsts = zinLDA::posterior(fit_Gibbs)
    theta <- posteriorEsts$theta
    beta <- posteriorEsts$beta

    out = list(beta = beta,
               theta = theta)
    return(out)
  }
}

