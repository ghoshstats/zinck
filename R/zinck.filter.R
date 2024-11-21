#' @title Zinck Filter for Variable Selection
#' @description Performs variable selection by fitting the augmented set of features to the response, using a glmnet or Random Forest model
#' @param X An OTU matrix with dimensions \eqn{D \times p}.
#' @param X_tilde A knockoff matrix corresponding to X with the same dimensions.
#' @param Y The response variable, either continuous or binary.
#' @param model A string; the model to be used ("glmnet" or "Random Forest").
#' @param fdr Numeric; the false discovery rate. Default is 0.1
#' @param offset Numeric; either 0 or 1. Default is 1
#' @param seed Numeric; the seed for reproducibility.
#' @param ntrees Numeric; the number of trees for Random Forest. Default is 1000.
#' @param tune_mtry Logical; whether to tune mtry in Random Forest. Default is FALSE.
#' @param mtry Numeric; the number of variables randomly sampled as candidates at each split in Random Forest. Default is NULL (auto-tuned).
#' @return A vector of selected variables at a target false discovery rate
#' @importFrom knockoff knockoff.threshold
#' @importFrom randomForest tuneRF randomForest importance
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats coef
#' @export
zinck.filter <- function(X, X_tilde, Y, model, fdr = 0.1, offset = 1, seed = NULL, ntrees = 1000, tune_mtry = FALSE, mtry = NULL) {
  # Check for input errors
  if (!is.matrix(X) || !is.matrix(X_tilde)) stop("X and X_tilde must be matrices.")
  if (ncol(X) != ncol(X_tilde)) stop("X and X_tilde must have the same number of columns.")
  if (length(Y) != nrow(X)) stop("Length of Y must match the number of rows in X.")
  if (!(model %in% c("glmnet", "Random Forest"))) stop("Invalid model. Choose 'glmnet' or 'Random Forest'.")
  if (!is.numeric(fdr) || fdr <= 0 || fdr >= 1) stop("fdr must be a numeric value between 0 and 1.")
  if (!is.numeric(offset) || !(offset %in% c(0, 1))) stop("offset must be 0 or 1.")
  if (!is.null(seed) && (!is.numeric(seed) || seed <= 0)) stop("seed must be a positive numeric value.")
  if (model == "Random Forest" && (!is.numeric(ntrees) || ntrees <= 0)) stop("ntrees must be a positive numeric value.")
  if (model == "Random Forest" && !is.logical(tune_mtry)) stop("tune_mtry must be TRUE or FALSE.")
  if (!is.null(mtry) && (!is.numeric(mtry) || mtry <= 0)) stop("mtry must be a positive numeric value.")
  
  X_aug <- cbind(X, X_tilde)
  
  if (model == "glmnet") {
    if (is.factor(Y) || length(unique(Y)) == 2) { # Binary response
      W <- stat.lasso_coefdiff_bin(X, X_tilde, Y)
    } else { # Continuous response
      W <- stat.lasso_coefdiff(X, X_tilde, Y)
    }
    T <- knockoff.threshold(W, fdr = fdr, offset = offset)
    selected_glm <- sort(which(W >= T))
    cat('\nSelected variables:\n')
    print(selected_glm)
    return(selected_glm)
  } else if (model == "Random Forest") {
    set.seed(seed)
    if (is.null(mtry)) { # Tune mtry if not provided
      if (tune_mtry) {
        bestmtry <- tuneRF(X_aug, Y, stepFactor = 1.5, improve = 1e-5, ntree = ntrees, trace = FALSE)
        mtry <- bestmtry[as.numeric(which.min(bestmtry[, "OOBError"])), 1]
      } else {
        mtry <- floor(sqrt(ncol(X_aug))) # Default mtry
      }
    }
    set.seed(seed)
    model_rf <- randomForest(X_aug, Y, ntree = ntrees, mtry = mtry, importance = TRUE)
    cf <- importance(model_rf)[, 1]
    W <- abs(cf[1:ncol(X)]) - abs(cf[(ncol(X) + 1):ncol(X_aug)])
    T <- knockoff.threshold(W, fdr = fdr, offset = offset)
    selected_rf <- sort(which(W >= T))
    cat('\nSelected variables:\n')
    print(selected_rf)
    return(selected_rf)
  }
}
