
<!-- README.md is generated from README.Rmd. Please edit that file -->

# icWGCNA

## Overview

Iterative Correcting Weighted Gene Co-expression Network Analysis
function for constructing a gene network from a gene expression matrix.
The algorithm:

    1. Constructs a signed wgcna network
    2. Drops correlated modules based on kurtosis
    3. Regresses out the largest community from the expression data
    4. Repeats steps 1-3 until a maximum number of communities or iterations is reached

Some differences from standard [WGNCA
(Horvath/Langfelder)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559):

-   Makes heavy use of
    [Rfast](https://cran.r-project.org/web/packages/Rfast/) to compute
    adjacencies and [topological overlap measure
    (TOM)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559),
    which enables iterative network creation on \> 20K features.
-   Uses [signed
    adjacency](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-10-327)
    in order to avoid possible distortions of community signatures
    (eigengenes).
-   Iteratively regresses out strongest community in order to facilitate
    discovery of communities possibly obscured larger module(s).
-   Clustering does not focus on merging communities but dropping to
    identify strongest module(s), keeping communities with higher
    membership kurtosis.
-   In addition to Pearson correlation, there is an option for Spearman
    correlation as the base constructing adjacency matrix measure in
    order to enable robust application in RNA-seq and micro-array data.
    Future updates may include mutual information and other measures of
    similarity.

## Installation

Install the released version of icWGCNA from BMS RStudio Package Manager
(BRAN):

``` r
install.packages("icWGCNA", 
                 repos = "http://pm.rdcloud.bms.com:4242/bms-cg-biogit-bran/latest")
```

Or install the development version from BMS BioGit with:

``` r
remotes::install_github(repo = "Systems-Immunology/icWGCNA", 
                        host = "https://biogit.pri.bms.comapi/api/v3")
#or use install_git
remotes::install_git('https://biogit.pri.bms.com/Systems-Immunology/icWGCNA.git')
```

## Lifecycle

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

icWGCNA is still in development.

## Example

EXAMPLE CODE TBD

``` r
library(icWGCNA)
## basic example code
```

## Code of Conduct

Please note that the icWGCNA project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
