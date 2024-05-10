# code to make testing data used
luad <- suppressMessages(
  UCSCXenaTools::getTCGAdata(
    project = "LUAD",
    mRNASeq = TRUE,
    mRNASeqType = "normalized",
    clinical = FALSE,
    download = TRUE
  )
)
ex <- data.table::fread(luad$destfiles, data.table = FALSE)
rownames(ex) <- ex[, 1]
mu <- apply(as.matrix(ex[, -1]), 1, mean)
SD <- apply(as.matrix(ex[, -1]), 1, sd)
n_picked <- 50
testing_data <- withr::with_seed(
  seed = 148654315,
  code = ex[mu > 2.5 & SD > 1.5, sample.int(ncol(ex), n_picked)]
)
# Adding rows with 0 SD, witch should get removed
testing_data <- rbind(
  rep(1, ncol(testing_data)),
  testing_data,
  rep(1, ncol(testing_data))
)
saveRDS(testing_data, file = testthat::test_path("fixtures", "testing_data.rds"))

# saving results file
testing_results <- icwgcna(testing_data, maxIt = 3, covCut = .66, mat_mult_method = "RcppEigen")
saveRDS(testing_results, file = testthat::test_path("fixtures", "testing_results.rds"))

# saving compute eigengene matrix results
testing_eigengene_matrix <- compute_eigengene_matrix(testing_data, testing_results$community_membership)
saveRDS(testing_eigengene_matrix, file = testthat::test_path("fixtures", "testing_eigengene_matrix.rds"))


# saving pangDB results
testing_pangDB <- data.table::fread(pangDB_link, showProgress = FALSE)
saveRDS(testing_pangDB,
        file = testthat::test_path("fixtures",
                                   "testing_pangDB.rds"))

testing_panglaoDB_enrichment <- compute_panglaoDB_enrichment(testing_results$community_membership,
  pangDB = testing_pangDB
)
saveRDS(testing_panglaoDB_enrichment,
        file = testthat::test_path("fixtures",
                                   "testing_panglaoDB_enrichment.rds"))

# saving MSigDB results
testing_MSigDB_enrichment <- withr::with_collate(
  "C",
  compute_MSigDB_enrichment(testing_results$community_membership)
)
saveRDS(testing_MSigDB_enrichment,
        file = testthat::test_path("fixtures",
                                   "testing_MSigDB_enrichment.rds"))


# saving xCell results
testing_xCell_enrichment <- compute_xCell_enrichment(
  testing_results$community_membership
)
saveRDS(testing_xCell_enrichment,
        file = testthat::test_path("fixtures",
                                   "testing_xCell_enrichment.rds"))

# saving UMAP ggplots results
custom_umap_specs <- umap::umap.defaults
custom_umap_specs$random_state <- 94124456
testing_UMAP_results <- make_network_umap(
    testing_results$community_membership,
    umap_specs = custom_umap_specs,
    community_labels = data.frame(community = 'mA1', lab = 'Extra')
  )

saveRDS(testing_UMAP_results,
        file = testthat::test_path("fixtures",
                                   "testing_UMAP_results.rds"))
saveRDS(list(layout = testing_UMAP_results$layout[,1:2]),
        file = testthat::test_path("fixtures",
                                   "UMAP_testing_layout.rds"))

