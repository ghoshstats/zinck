% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.zinck.R
\name{fit.zinck}
\alias{fit.zinck}
\title{fit.zinck}
\usage{
fit.zinck(
  X,
  num_clusters,
  method = c("ADVI", "Gibbs"),
  seed = NULL,
  init_values = NULL,
  alpha_param = 0.1,
  boundary_correction = FALSE,
  importance_resampling = FALSE,
  elbo_samples = 500,
  prior_ZIGD = FALSE
)
}
\arguments{
\item{X}{An OTU matrix with dimensions \eqn{D \times p}.}

\item{num_clusters}{An integer specifying the number of clusters.}

\item{method}{A character string, either "ADVI" or "Gibbs", specifying the method to fit the model.}

\item{seed}{An integer used to set the seed for reproducibility.}

\item{init_values}{A list of initial values for the ADVI algorithm. This parameter is optional and should be used only with the ADVI method. If NULL, the algorithm uses default initialization.}

\item{alpha_param}{A positive real. The symmetric smoothed Dirichlet parameter for the cluster distributions, default=0.1}

\item{boundary_correction}{A logical value. If TRUE, it adds and subtracts small numbers so that the log-likelihood doesn't blow up at the boundary points 0 or 1, default = FALSE.}

\item{importance_resampling}{Logical scalar (defaulting to FALSE) indicating whether to do importance resampling to adjust the draws at the optimum to be more like draws from the posterior distribution.}

\item{elbo_samples}{A positive integer. The number of samples for Monte Carlo estimate of ELBO (objective function), defaulting to 100. (ELBO stands for "the evidence lower bound".)name description}

\item{prior_ZIGD}{A logical value. If TRUE and method is "ADVI", the model will impose Gamma priors on the ZIGD hyperparameters.}
}
\value{
A list containing the posterior estimates of beta and theta.
}
\description{
Fit the Zinck model to the data using either ADVI or Gibbs sampling methods.
}
\examples{
# Generate a random OTU matrix
X <- matrix(rpois(20, lambda = 5), nrow = 5)

# Fit the Zinck model using ADVI with default initialization
fit <- fit.zinck(X, num_clusters = 3, method = "ADVI", alpha_param=1)

# Fit the Zinck model using ADVI with custom initial values
init_vals <- list(theta = array(runif(15), dim = c(5, 3)),
                  zeta = array(runif(12), dim = c(3, 4)))
fit_with_init <- fit.zinck(X, num_clusters = 3, method = "ADVI", init_values = init_vals, alpha_param=1)

}
\references{
Kucukelbir, A., Tran, D., Ranganath, R., Gelman, A., and Blei, D. (2017).
Automatic Differentiation Variational Inference.
Journal of Machine Learning Research 18(14).

Deek, R., and Li, H. (2021).
A Zero-Inflated Latent Dirichlet Allocation Model for Microbiome Studies.
Frontiers in Genetics 11.
}
