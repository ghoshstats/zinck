#' @title Optimal Number of Clusters based on Jensen-Shannon Divergence
#' @description This function identifies the optimal number of clusters for fitting a zinck model
#' using Jensen-Shannon Divergence. The function fits the model for various values of clusters and calculates
#' the average Jensen-Shannon Divergence for each, returning the number of clusters that minimizes this value.
#' @param X A OTU matrix with dimensions \eqn{D \times p}.
#' @param kmin Numeric; the minimum number of clusters to be considered.
#' @param kmax Numeric; the maximum number of clusters to be considered.
#' @param seed_list List of numeric values; seeds for reproducibility for each k value, default is NULL.
#' @return The optimal number of clusters, K, that minimizes the average Jensen-Shannon Divergence.
#' @examples
#' \dontrun{
#' data("combo")
#' X <- combo[["abund_list"]][["Genus"]]
#' result <- optimal_k(X, kmin=2, kmax=10, seed_list=list(1,1,2,4,123,6,123,123,12))
#' print(result)
#' }
#'
#' @export
#' @importFrom rstan stan_model vb
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal theme

optimal_k <- function(X, kmin, kmax, seed_list = NULL)
{
  # Check for input errors
  if (!is.matrix(X)) stop("X must be a matrix.")
  if (!is.numeric(kmin) || kmin <= 0) stop("kmin must be a positive numeric value.")
  if (!is.numeric(kmax) || kmax <= 0) stop("kmax must be a positive numeric value.")
  if (kmin > kmax) stop("kmin must be less than or equal to kmax.")
  if (!is.null(seed_list)) {
    if (!is.list(seed_list)) stop("seed_list must be a list of numeric values.")
    if (length(seed_list) != kmax - kmin + 1) stop("Length of seed_list must be equal to the range of k values.")
    if (any(sapply(seed_list, function(x) !is.numeric(x) || x <= 0))) stop("Each seed in seed_list must be a positive numeric value.")
  }
  # Function to calculate Jensen-Shannon Divergence
  js_divergence <- function(p, q) {
    m <- 0.5 * (p + q)
    0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
  }

  # Loop over various values of K
  K_values <- seq(kmin, kmax, by=1 )
  js_values <- numeric(length(K_values))
  dlt <- rep(0,ncol(X))
  for(t in (1:ncol(X)))
  {
    dlt[t] <- 1-mean(X[,t]>0)
  }

  for (i in seq_along(K_values)) {
    num_clusters <- K_values[i]
    zinck_stan_data <- list(
      K = num_clusters,
      V = ncol(X),
      D = nrow(X),
      n = X,
      delta = dlt
    )
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
    stan.model = stan_model(model_code = zinck_code)
    set.seed(seed_list[[i]])
    fit_zinck <- suppressWarnings(vb(stan.model, data=zinck_stan_data, algorithm="meanfield", importance_resampling=TRUE, iter=10000,tol_rel_obj=0.01,elbo_samples=500))

    # Extract cluster-term distributions
    beta <- fit_zinck@sim[["est"]][["beta"]]

    # Calculate average JS Divergence for this K
    js_sum <- 0
    count <- 0

    for (cluster1 in 1:(num_clusters - 1)) {
      for (cluster2 in (cluster1 + 1):num_clusters) {
        js_sum <- js_sum + js_divergence(beta[,cluster1], beta[,cluster2])
        count <- count + 1
      }
    }

    js_avg <- js_sum / count
    js_values[i] <- js_avg
  }

  data <- data.frame(K_values, js_values)

  ggplot(data, aes(x = K_values, y = js_values)) +
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
  return(K_values[which.min(js_values)])
}
