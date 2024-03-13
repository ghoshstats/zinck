# zinck

## Overview
`zinck` is a novel knockoff-based framework specifically designed for microbiome data analysis. Microbiome datasets are often high-dimensional, compositional, and zero-inflated, posing unique challenges for statistical analysis and interpretation. `zinck` addresses these challenges by employing a flexible generative model that effectively captures the zero-inflation and complex dependence structure characteristic of microbial communities. This approach allows for simultaneous variable selection and false discovery rate (FDR) control, making `zinck` particularly suited for taxonomic variable selection in microbiome studies. 


## Installation Guide for `zinck` 

### System Requirements

Before installing `zinck`, ensure your system meets these requirements:
- R version 4.0.0 or higher
- Rtools (for Windows users)
- A C++ compiler compatible with R

## Dependencies

`zinck` relies on several R packages. Most dependencies are automatically installed when you install `zinck`. However, special attention is needed for `rstan` due to version compatibility.

### Core Dependencies

- `dplyr`
- `reshape2`
- `knockoff`
- `glmnet`
- `randomForest`
- `caret`
- `rstan` (version 2.21.8)
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

## Installing `rstan` (Version 2.21.8)

To use `zinck`, it's critical to install `rstan` version 2.21.8. Follow these steps:

1. **Remove Older Versions**: If you have a different version of `rstan` installed, remove it using `remove.packages("rstan")`.

2. **Install `rstan` 2.21.8**: You can install the exact version from the RStan archives. Use the following command in R:

    ```r
    packageurl <- "https://cran.r-project.org/src/contrib/Archive/rstan/rstan_2.21.8.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
    ```

    Note: Installation from source requires Rtools (Windows) or the appropriate development tools for other operating systems.

3. **Verify Installation**: Ensure that the correct version of `rstan` is installed by running `packageVersion("rstan")` in R. You should see `2.21.8` as the output.

## Installing `zinck`

After ensuring that the correct version of `rstan` is installed, you can proceed with installing `zinck`.

1. **Install `devtools` Package**: If not already installed, you need `devtools` to install `zinck` from GitHub.

    ```r
    if (!require("devtools")) install.packages("devtools")
    ```

2. **Install `zinck`**:

    ```r
    devtools::install_github("ghoshstats/zinck", build_vignettes = TRUE)
    ```



## Usage

Once installed, you can load `zinck` in R:

```r
library(zinck)
```

## Examples and vignettes

To get started with zinck, explore the package vignettes and documentation:

```r
# View available vignettes
browseVignettes(package = "zinck")

# View package documentation
?zinck
```


