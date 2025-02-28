#' zinck: A Package for using zero-inflated compositional knockoffs for microbial variable selection
#'
#' \strong{Description:}
#' zinck is a novel knockoff-based framework tailored for microbiome data exploiting a flexible generative model.
#' It can properly capture the zero-inflation and complex dependence structure among microbes, enjoying the property of simultaneous variable selection and FDR control.
#'
#' \strong{Main Functions:}
#' \describe{
#'   \item{\code{\link{fit.zinck}}}{Fits the zero-inflated hierarchical model to a count matrix using ADVI or Gibbs sampling.}
#'   \item{\code{\link{generateKnockoff}}}{Generates a knockoff copy of the original matrix once the posterior estimates of the latent parameters are obtained.}
#'   \item{\code{\link{optimal_k}}}{Finds the optimal number of clusters minimizing the Jensen-Shannon Divergence.}
#'   \item{\code{\link{log_normalize}}}{Performs log-normalization of a given matrix.}
#'   \item{\code{\link{draw_heatmap}}}{Creates a heatmap for a microbial sample taxa matrix.}
#'   \item{\code{\link{zinck.filter}}}{Performs FDR-controlled variable selection by fitting a glmnet or a Random Forest model.}
#' }
#'
#' \strong{Datasets:}
#' \describe{
#'   \item{\code{\link{count.genus}}}{Genus level CRC data.}
#'   \item{\code{\link{count}}}{Species level CRC data.}
#' }
#'
#' For a complete list of functions and datasets, see the \code{INDEX} section.
#'
#' @docType package
#' @name zinck-package
#' @aliases zinck
#' @keywords package
NULL
