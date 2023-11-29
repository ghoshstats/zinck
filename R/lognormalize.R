#' Log-Normalization for Microbiome Compositional Data
#'
#' This function performs log-normalization on a given matrix, typically used 
#' in microbiome data analysis. In microbiome studies, data are often compositional 
#' and contain many zero counts. Log-normalization, with the addition of a pseudo-count, 
#' is a standard approach to handle zeros and maintain the compositional nature of the data. 
#' This function takes a numeric matrix, adds a pseudo-count to zero values, 
#' and then applies a log transformation, preserving the relative proportions 
#' in the data.
#'
#' The function first adds a pseudo-count of 0.5 to zero entries to handle zeros 
#' in the data. It then divides each count by the total counts in its row (sample) 
#' to make the data compositional. Finally, it applies a natural logarithm 
#' transformation to the normalized data.
#'#'
#' @param X A numeric matrix to log-normalize.
#' @return A matrix of the same dimensions as \code{X}, with each element
#'   being the log-normalized value of the corresponding element in \code{X}.
#'   The base of the logarithm used for normalization is e (natural logarithm).
#'
#' @examples
#' # Create a sample matrix
#' mat <- matrix(1:4, nrow = 2)
#' 
#' # Perform log-normalization
#' log_normalized_mat <- log_normalize(mat)
#'
#' @export
#'
log_normalize <- function(X) {
  # Ensure input is a matrix
  if (!is.matrix(X)) {
    stop("Input must be a matrix.")
  }
  
  # Ensure all elements in the matrix are numeric
  if (!all(is.numeric(X))) {
    stop("Matrix contains non-numeric values.")
  }
  
  # Add a pseudo-count of 0.5 to zero values
  X[X == 0] <- 0.5
  
  # Make compositional by dividing each value by its row sum
  X <- sweep(X, 1, rowSums(X), FUN="/")
  
  # Apply log transformation
  Z <- log(X)
  
  return(Z)
}

