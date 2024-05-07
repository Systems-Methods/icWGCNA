# grabbing large file for testing
testing_data <- readRDS(test_path("fixtures", "testing_data.rds"))
testing_results <- readRDS(test_path("fixtures", "testing_results.rds"))

testing_pangDB <-
  readRDS(test_path("fixtures", "testing_pangDB.rds"))
testing_panglaoDB_enrichment <-
  readRDS(test_path("fixtures", "testing_panglaoDB_enrichment.rds"))

testing_MSigDB_enrichment <-
  readRDS(test_path("fixtures", "testing_MSigDB_enrichment.rds"))

testing_xCell_enrichment <-
  readRDS(test_path("fixtures", "testing_xCell_enrichment.rds"))

testing_eigengene_matrix <-
  readRDS(test_path("fixtures", "testing_eigengene_matrix.rds"))

testing_UMAP_results <-
  readRDS(test_path("fixtures", "testing_UMAP_results.rds"))
UMAP_testing_layout <-
  readRDS(test_path("fixtures", "UMAP_testing_layout.rds"))

testing_Seurat <-
  readRDS(test_path("fixtures", "testing_Seurat.rds"))

