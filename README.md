
<!-- README.md is generated from README.Rmd. Please edit that file -->

# icWGCNA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/Systems-Methods/icWGCNA/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Systems-Methods/icWGCNA?branch=main)
<!-- badges: end -->

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
    discovery of communities possibly obscured by larger module(s).
-   Clustering does not focus on merging communities but dropping to
    identify strongest module(s), keeping communities with higher
    membership kurtosis.
-   In addition to Pearson correlation, there is an option for Spearman
    correlation as the base constructing adjacency matrix measure in
    order to enable robust application in RNA-seq and micro-array data.
    Future updates may include mutual information and other measures of
    similarity.

## Installation

Install the released version of icWGCNA from Github:

``` r
remotes::install_github(repo = "Systems-Methods/icWGCNA")
```

Or install the development version from BMS BioGit with:

``` r
remotes::install_github(repo = "Systems-Methods/icWGCNA", 
                        ref = "develop")
```

Or:

``` r
remotes::install_git(
  'https://github.com/Systems-Methods/icWGCNA.git'
)
```

## Example

### Example Data

First we can use the
[{UCSCXenaTools}](https://cran.r-project.org/web/packages/UCSCXenaTools/vignettes/USCSXenaTools.html)
package to download an example mRNASeq dataset.

``` r
luad <- UCSCXenaTools::getTCGAdata(project = "LUAD", 
                                   mRNASeq = TRUE, 
                                   mRNASeqType = "normalized",
                                   clinical = FALSE, 
                                   download = TRUE)

ex <- as.matrix(data.table::fread(luad$destfiles), rownames = 1)
```

### Main Results

Next we can run the main `icwgcna()` function

``` r
library(icWGCNA)

results <- icwgcna(ex,   
                   expo = 6,
                   Method = "pearson",
                   q = 0.3,
                   maxIt = 10,
                   maxComm = 100,
                   corCut = 0.8,
                   covCut = 0.33,
                   mat_mult_method = "Rfast")
```

### Downstream Analysis

Finally, downstream analysis can be run on the Iterative Correcting
Weighted Gene Co-expression Network Analysis results.

``` r

compute_eigengene_matrix(ex, 
                         membership_matrix = results$community_membership, 
                         cutoff = 5,
                         pc_flag = TRUE)



pangDB <- data.table::fread(pangDB_link)
compute_panglaoDB_enrichment(results$community_membership,
                             K = 100,
                             memb_cut = 0.65,
                             pangDB = pangDB,
                             prolif = prolif_names)
```

## Code of Conduct

Please note that the icWGCNA project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.