
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
    2. Drops correlated modules based on membership kurtosis
    3. Regresses out the strongest community from the expression data
    4. Repeats steps 1-3 until a maximum number of communities or iterations is reached

Some differences from standard [WGNCA
(Horvath/Langfelder)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559):

- Makes heavy use of
  [Rfast](https://cran.r-project.org/web/packages/Rfast/) to compute
  adjacencies and [topological overlap measure
  (TOM)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559),
  which enables iterative network creation on \> 20K features.
- Always uses [signed
  adjacency](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-10-327)
  in order to avoid possible distortions of community signatures
  (eigengenes).
- Iteratively regresses out strongest community in order to facilitate
  discovery of communities possibly obscured by larger module(s).
- Clustering does not focus on merging communities but dropping to
  identify strongest module(s), keeping communities with higher
  membership kurtosis.
- In addition to Pearson correlation, there is an option for Spearman
  correlation as the base constructing adjacency matrix measure in order
  to enable robust application in RNA-seq and micro-array data. Future
  updates may include mutual information and other measures of
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

#### Compute community signatures (eigengenes)

``` r

eigengene_mat <- compute_eigengene_matrix(
  ex, 
  membership_matrix = results$community_membership, 
  cutoff = 5,
  pc_flag = TRUE
)
```

#### Compute panglaoDB collection enrichments for each community

``` r
pangDB <- data.table::fread(pangDB_link)
panglaoDB_enrichment <- compute_panglaoDB_enrichment(
  results$community_membership,
  K = 100,
  memb_cut = 0.65,
  pangDB = pangDB,
  prolif = prolif_names, 
  p_cut = 0.001
)
```

#### Compute MSigDB collection enrichments for each community

``` r
MSigDB_enrichment <- compute_MSigDB_enrichment(
  results$community_membership,
  K = 100,
  memb_cut = .65,
  cats = c("H", "C3", "C6", "C7", "C8"), 
  p_cut = 0.001
)
```

#### Compute xCell collection enrichments for each community

``` r
xCell_enrichment <- compute_xCell_enrichment(
  results$community_membership,
  K = 100,
  memb_cut = .65, 
  p_cut = 0.001
)
```

#### Display UMAP of Community Membership with text overlays

``` r
network_umap <- make_network_umap(
  results$community_membership,
  community_memb_cut_main = 0.7,
  community_n_main = 20,
  community_memb_cut_secondary = 0.8,
  community_n_secondary = 5,
  gene_memb_cut_main = 0.75,
  gene_memb_cut_secondary = 0.65,
  community_labels = NULL,
  umap_specs = umap::umap.defaults
)

network_umap$umap_w_annotation
```

#### Identify Top Genes, with Values, of all Communities

``` r
top_genes <- display_top_genes(
  results$community_membership,
  K = 10,
  output = "both"
)
```

#### Identify Top Gene of Communities that are unique (only belong to one community)

``` r
unique_top_genes <- find_unique_top_genes(
  results$community_membership,
  K = 10,
  maxIt = 10
)
```

#### Map icWGCNA Eigengenes on a Seurat Object

``` r
unique_top_genes <- map_eigengenes_on_seurat(
  SeuratObject::pbmc_small,
  results$community_membership,
  cutoff_method = "top_gene",
  top_genes_cutoff = 20
)
```

## Code of Conduct

Please note that the icWGCNA project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
