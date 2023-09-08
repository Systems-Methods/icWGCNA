test_that("gene_mapping success", {

  tmp_expr_data <- data.frame(
    ID_1 = 1:10,
    ID_2 = 10:1,
    ID_3 = c(2,3,6,7,10,1000,1,1,1,1)
  )
  rownames(tmp_expr_data) <- paste0("Probe_", 1:10)
  tmp_mapping_file <- data.frame(
    probe = paste0("Probe_", c(1,1,2,2,2,3,4,5,6,7)),
    gene_symbol = paste0("Gene_", c(1,2,1,3,4,5,3,1,6,3))
  )


  default_results <- purrr::quietly(
    ~ gene_mapping(tmp_expr_data, tmp_mapping_file)
  )()

  test_results <- data.frame(
    ID_1 = c(5, 1, 4, 2, 3, 6),
    ID_2 = c(6, 10, 7, 9, 8, 5),
    ID_3 = c(10, 2, 7, 3, 6, 1000)
  )
  rownames(test_results) = paste0("Gene_", 1:6)

  expect_equal(
    default_results$result,
    test_results
  )
  expect_equal(
    default_results$messages,
    "3 exprs_data rows could not be linked using the mapping file, resulting in 7 rows. \nThese rows link to 6 distinct gene symbols. \nWill compress duplicate rows using the highest_mean method with no transformation.\n"
  )

  test_results <- data.frame(
    ID_1 = c(2.15443469003188, 1, 3.82586236554478,
             2, 3, 6),
    ID_2 = c(8.14325284978472, 10, 6.31635959765638, 9,
             8, 5),
    ID_3 = c(3.91486764116886, 2, 2.75892417638112, 3, 6,
             1000)
  )
  rownames(test_results) = paste0("Gene_", 1:6)

  expect_equal(
    gene_mapping(tmp_expr_data, tmp_mapping_file, "mean", "log",
                 verbose = FALSE),
    test_results
  )

})


test_that("gene_mapping no dups", {

  tmp_expr_data <- data.frame(
    ID_1 = 1:10,
    ID_2 = 10:1,
    ID_3 = c(2,3,6,7,10,1000,1,1,1,1)
  )
  rownames(tmp_expr_data) <- paste0("Probe_", 0:9)
  tmp_mapping_file <- data.frame(
    probe = paste0("Probe_", 0:9),
    gene_symbol = paste0("Gene_", 0:9)
  )

  default_results <- purrr::quietly(
    ~ gene_mapping(tmp_expr_data, tmp_mapping_file)
  )()

  test_results <- tmp_expr_data
  rownames(test_results) = paste0("Gene_", 0:9)

  expect_equal(
    default_results$result,
    test_results
  )
  expect_equal(
    default_results$messages,
    "0 exprs_data rows could not be linked using the mapping file, resulting in 10 rows. \nThese rows link to 10 distinct gene symbols. \nWill compress duplicate rows using the highest_mean method with no transformation.\n"
  )
})

test_that("gene_mapping errors", {
  tmp_expr_data <- data.frame(
    ID_1 = 1:10,
    ID_2 = 10:1,
    ID_3 = c(0,3,6,7,10,1000,1,1,1,1)
  )
  rownames(tmp_expr_data) <- paste0("Probe_", 0:9)
  tmp_mapping_file <- data.frame(
    probe = paste0("Probe_", 0:9),
    gene_symbol = paste0("Gene_", 0:9)
  )

  expect_error(
    gene_mapping(tmp_expr_data, tmp_mapping_file, "mean", "log",
                 verbose = FALSE),
    "Can't do log transformation with exprs_data values <= 0"
  )
  tmp_mapping_file <- data.frame(
    probe = paste0("Probe_", 10:19),
    gene_symbol = paste0("Gene_", 0:9)
  )
  expect_error(
    gene_mapping(tmp_expr_data, tmp_mapping_file, verbose = FALSE),
    "Could not link rownames\\(exprs_data\\) to mapping_file\\[,1\\]\\!"
  )

})
