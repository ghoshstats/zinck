# zinck

## Overview
`zinck` is a novel knockoff-based framework specifically designed for microbiome data analysis. Microbiome datasets are often high-dimensional, compositional, and zero-inflated, posing unique challenges for statistical analysis and interpretation. `zinck` addresses these challenges by employing a flexible generative model that effectively captures the zero-inflation and complex dependence structure characteristic of microbial communities. This approach allows for simultaneous variable selection and false discovery rate (FDR) control, making `zinck` particularly suited for taxonomic variable selection in microbiome studies. 


# Installation Guide for `zinck`

## System Requirements

Before installing `zinck`, ensure your system meets these requirements:
- R version 4.0.0 or higher
- Rtools (for Windows users)
- A C++ compiler compatible with R

## Dependencies

`zinck` relies on several R packages. Most dependencies are automatically installed when you install `zinck`. However, special attention is needed for `rstan` and `StanHeaders` due to version compatibility.

### Core Dependencies

- `dplyr`
- `reshape2`
- `knockoff`
- `glmnet`
- `randomForest`
- `caret`
- `rstan` (version 2.21.8)
- `StanHeaders` (version 2.26.25)
- `stats`
- `fitdistrplus`
- `ggplot2`
- `MLmetrics`
- `phyloseq`
- `GUniFrac`
- `kosel`
- `gridExtra`
- `zinLDA`

### Suggested Packages

- `knitr`
- `rmarkdown`

## Installing `rstan` and `StanHeaders`

To use `zinck`, it's critical to install specific versions of `rstan` and `StanHeaders`. Follow these steps:

1. **Remove Older Versions**: If you have different versions of `rstan` or `StanHeaders` installed, remove them using:
   
    ```r
    remove.packages(c("rstan", "StanHeaders"))
    ```

2. **Install `StanHeaders` 2.26.25**:
   
    ```r
    packageurl <- "https://cran.r-project.org/src/contrib/Archive/StanHeaders/StanHeaders_2.26.25.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
    ```

3. **Install `rstan` 2.21.8**:

    ```r
    packageurl <- "https://cran.r-project.org/src/contrib/Archive/rstan/rstan_2.21.8.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
    ```

    Note: Installation from source requires Rtools (Windows) or the appropriate development tools for other operating systems.

4. **Verify Installation**: Ensure that the correct versions of `rstan` and `StanHeaders` are installed by running:

    ```r
    packageVersion("rstan") # Should return 2.21.8
    packageVersion("StanHeaders") # Should return 2.26.25
    ```

## Installing `zinck`

After ensuring that the correct versions of `rstan` and `StanHeaders` are installed, you can proceed with installing `zinck`.

1. **Install `devtools` Package**: If not already installed, you need `devtools` to install `zinck` from GitHub.

    ```r
    if (!require("devtools")) install.packages("devtools")
    ```

2. **Install `zinck`**:

    ```r
    devtools::install_github("ghoshstats/zinck", build_vignettes = TRUE)
    ```

## Post-Installation

Once `zinck` is successfully installed, load it in R to begin your analysis:

```r
library(zinck)
```

## Examples and vignettes

To get started with zinck, explore the package documentation:

```r
# View package documentation
?zinck
```

## Getting Started

After installation, you can verify and get a feel for the package's capabilities through a practical example. This guide walks you through loading a dataset included with zinck, selecting samples, generating knockoffs, and performing variable selection. This example utilizes the `count.Rdata` dataset included in the zinck package. Ensure it's loaded correctly into your R session.

```r
# Load the dataset
data("count", package = "zinck")

# Ordering the columns by decreasing abundance
dcount <- count[, order(decreasing = TRUE, colSums(count, na.rm = TRUE), apply(count, 2L, paste, collapse = ''))][, 1:300]

# Random Subject Selection
set.seed(1)
sel_index <- rbinom(nrow(dcount), size = 1, prob = 0.5)
selected_samples <- which(sel_index == 1)
X <- dcount[selected_samples, ]

# Prepare signals so that each block sums up to zero
Five1 <- c(-3,3,2.5,-1,-1.5)
Five2 <- c(3,3,-2,-2,-2)
Five3 <- c(1,-1,3,-2,-1)
Five4 <- c(-1,1,2,-1,-1)
Five5 <- c(3,3,-3,-2,-1)
Five6 <- c(-1,1,-2,1,1)
Five_all <- c(Five1,Five2,Five3,Five4,Five5,Five6)
randBeta <- rep(0,300)

set.seed(1)
rand_indices <- sample(1:200,size=30,replace=FALSE) # Randomly Injecting the signals among the first 200 most abundant species

set.seed(1)
randBeta[rand_indices] <- sample(Five_all, size=30, replace=FALSE) # Randomly assigning each signal amplitude to these indices



# Generate response
n = nrow(X)
p = ncol(X)
W <- log_normalize(X)
set.seed(1)
eps = rnorm(n, mean = 0, sd = 1)
Y <- W %*% randBeta + eps # Generate Y using a log-contrast model

# Fit the zinck model
species_fit <- fit.zinck(X, num_clusters = 6, method = "ADVI", seed = 123, alpha_param = 0.1)

# Extract model parameters
theta <- species_fit$theta
beta <- species_fit$beta

# Generate knockoff features
X_tilde <- generateKnockoff(X, theta, beta, seed = 1)
W_tilde <- log_normalize(X_tilde)

# Perform variable selection
selected_species <- zinck.filter(W, W_tilde, Y, model = "glmnet", fdr = 0.1, offset = 1)


index <- rand_indices
index_est <- selected_species
neg_index <- (1:p)[-c(index)]
if(length(index_est) == 0) {
  neg_index_est <- 1:p
} else {
  neg_index_est <- (1:p)[-c(index_est)]
}

# Evaluating model performance
TP <- sum(selected_species %in% rand_indices) # True Positives
FP <- sum(!selected_species %in% rand_indices) # False Positives
TN <- sum(!neg_index_est %in% rand_indices) # True Negatives
FN <- length(rand_indices) - TP # False Negatives

estimated_FDR <- FP / (FP + TP) # Evaluating the empirical False Discovery Rate
estimated_power <- TP / (TP + FN) # Evaluating the empirical Power or TPR

print(paste("Estimated FDR:", estimated_FDR))
print(paste("Estimated Power:", estimated_power))
```
Explore the package vignettes for detailed tutorials and examples:

```r
# View available vignettes
browseVignettes(package = "zinck")
```

For any issues or further assistance, consult the zinck pdf manual or visit the GitHub repository's Issues section.







