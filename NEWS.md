# icWGCNA 0.3.0

### New Features

- `map_eigengenes_on_seurat()` to map eigengenes on a Seurat Object

# icWGCNA 0.2.5

### New Features

- `find_unique_top_genes()` to identify top gene of communities that are unique 
(only belong to one community)

### Other Changes

- Removed uncorrected_community_signature output from `icwgcna()`
    - Now users should run `compute_eigengene_matrix()` after 
    `icwgcna()`if uncorrected community signatures are needed.

# icWGCNA 0.2.4

- Fixed bug in `compute_eigengene_matrix()` when all genes are above cutoff

# icWGCNA 0.2.3

### New Features

`expression_compression()` function for converting between gene types (i.e. ENTREZ to
Hugo), which many options for compression of duplicate rows.

### Other Changes

- Minor bug fixes
- Switching to  `fastcluster::hclust()` from stats pacakge (#16)

# icWGCNA 0.2.2

- Bug fixes and internal function improvements
- Add p value output for `compute_MSigDB_enrichment()`
- Add full_metaGenes and full_eigenGenes output for `icwgcna()`

# icWGCNA 0.2.1

Bug fixes and checking for PC1 > 35% in input expression

`icWGCNA` returns `uncorrected_community_signature` instead of 
`full_community_membership` and `full_community_signature`

# icWGCNA 0.2.0

### New Features

#### Enrichments Functions
- `compute_MSigDB_enrichment()` (with parallel and distributed processing)
- `compute_xCell_enrichment()`

#### UMAP Plotting
- `make_network_umap()`

# icWGCNA 0.1.1

Released on public Github. 

# icWGCNA 0.1.0

Initial release to be shared publicly

#### Bug Fixes

Better panglaoDB variable name matches and error catching (#18)

Downstream functions improvement and error catching

# icWGCNA 0.0.0.9002

#### Bug Fixes

Fixed error when genes have 0 standard deviation (#16)


# icWGCNA 0.0.0.9001

First Version

