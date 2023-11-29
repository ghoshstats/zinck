# zinck

## Overview
`zinck` is a novel knockoff-based framework specifically designed for microbiome data analysis. Microbiome datasets are often high-dimensional, compositional, and zero-inflated, posing unique challenges for statistical analysis and interpretation. `zinck` addresses these challenges by employing a flexible generative model that effectively captures the zero-inflation and complex dependence structure characteristic of microbial communities. This approach allows for simultaneous variable selection and false discovery rate (FDR) control, making `zinck` particularly suited for taxonomic variable selection in microbiome studies. 

## Installation

You can install the latest development version of `zinck` from GitHub:

```r
# Install the development version from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ghoshstats/zinck",build_vignettes = TRUE)
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


