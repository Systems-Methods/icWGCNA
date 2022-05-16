
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

Some differences from standard WGNCA (Horvath/Langfelder):

-   Makes heavy use of [Rfast
    package](https://cran.r-project.org/web/packages/Rfast/) to compute
    adjacencies and TOM to enable iterative network creation on \> 20K
    features.
-   Uses signed adjacency in order to avoid possible distortions of
    community signatures (eigengenes).
-   Iteratively regresses out strongest community in order to facilitate
    discovery of communities possibly obscured larger module(s).
-   Clustering does not focus on merging communities but dropping to
    identify strongest module(s).
-   Enables Spearman correlation for constructing adjacency matrix
    instead of Pearson to enable robust application in RNA-seq and
    micro-array data. Future updates may include mutual information

## Installation

You can currently install the development version of icWGCNA from the
BMS BioGit repository:

``` r
remotes::install_github("masonm3/icWGCNA", 
                        host = "https://biogit.pri.bms.com")
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
