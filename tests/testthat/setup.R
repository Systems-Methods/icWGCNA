# grabbing large file for testing
testing_data <- readRDS(test_path("fixtures", "testing_data.rds"))
testing_results <- readRDS(test_path("fixtures", "testing_results.rds"))
testing_enrichment <- readRDS(test_path("fixtures", "testing_enrichment.rds"))
testing_eigengene_matrix <- readRDS(
  test_path("fixtures", "testing_eigengene_matrix.rds")
)

testing_pangDB <- data.table::fread(pangDB_link, showProgress = FALSE)
