#' Generate Knockoff Copy of Microbial Sample Data
#'
#' This function generates a knockoff copy of microbial sample data given a matrix `X`
#' and matrices `Theta` and `Beta`. The function adjusts the column structure of `Beta`
#' to match `X`, generates samples based on `Theta` and `Beta`, and then compiles these
#' into a knockoff count matrix.
#'
#' @param X A numeric matrix representing original microbial sample data.
#' @param Theta A numeric matrix representing cluster mixing probabilities.
#' @param Beta A numeric matrix representing feature proportions for each cluster.
#' @param seed An optional integer seed for reproducibility of random generation.
#' @return A numeric matrix representing the knockoff copy of the microbial sample data.
#'
#' @examples
#' X <- matrix(runif(40), nrow = 10)
#' colnames(X) <- paste("Taxa", 1:ncol(X))
#' Theta <- matrix(runif(30), nrow = 10)
#' Beta <- matrix(runif(20), nrow = 5)
#' knockoff_data <- generateKnockoff(X, Theta, Beta)
#'
#' @export
#'
generateKnockoff <- function(X, Theta, Beta, seed = NULL) {
  D=nrow(Theta); V=ncol(Beta); N=rowSums(X)
  if(V == ncol(X)){
    colnames(Beta) <- colnames(X)
  }else if(V < ncol(X)){
    # Get the names of columns that are missing in Beta
    missing_cols <- setdiff(colnames(X), colnames(Beta))

    # For each missing column, add a column of zeros to Beta
    for (col in missing_cols) {
      Beta <- cbind(Beta, 0)
      colnames(Beta)[ncol(Beta)] <- col
    }

    # Reorder the columns of Beta to match the order in X
    Beta <- Beta[, colnames(X)]
  }
  # generate 1 sample
  generateSample <- function(N, theta, beta) {

    sample <- vector(length = N)
    z_d <- vector(length = N)
    for (n in 1:N) {
      z_n <- stats::rmultinom(1, 1, theta)
      w_n <- stats::rmultinom(1, 1, beta[which(z_n == 1),])
      sample[n] <- colnames(beta)[which(w_n == 1)]
      z_d[n] <- which(z_n == 1)
      names(z_d)[n] <- sample[n]
    }
    return(list(sample = sample, z = z_d))
  }

  # generate n samples
  cohort <- vector(mode = "list", length = D)
  z <- vector(mode = "list", length = D)
  for (d in 1:D) {
    sample.d <- generateSample(N[d], Theta[d, ], Beta)
    cohort[[d]] <- sample.d[["sample"]]
    z[[d]] <- sample.d[["z"]]
  }

  # collapse list of vectors into list of count tables
  sampleTaxaFreq <- lapply(cohort, table)

  # count matrix
  x_tilde <- matrix(data = 0, nrow = D, ncol = ncol(X))
  #rownames(x_tilde) <- rownames(Theta)
  colnames(x_tilde) <- colnames(Beta)

  for (d in 1:D) {
    x_tilde[d, names(sampleTaxaFreq[[d]])] <- sampleTaxaFreq[[d]]
  }

  return(x_tilde)
}

