test_that("expression_compression success defaults", {

  tmp_expr_data <- data.frame(
    ID_1 = 1:10,
    ID_2 = 10:1,
    ID_3 = c(2,3,6,7,10,1000,1,1,1,1)
  )
  rownames(tmp_expr_data) <- paste0("Probe_", 1:10)
  tmp_mapping_df <- data.frame(
    probe = paste0("Probe_", c(1,1,2,2,2,3,4,5,6,7)),
    gene_symbol = paste0("Gene_", c(1,2,1,3,4,5,3,1,6,3))
  )

  expect_snapshot(expression_compression(tmp_expr_data, tmp_mapping_df))

})


test_that("expression_compression success log mean", {

  tmp_expr_data <- data.frame(
    ID_1 = 1:10,
    ID_2 = 10:1,
    ID_3 = c(2,3,6,7,10,1000,1,1,1,1)
  )
  rownames(tmp_expr_data) <- paste0("Probe_", 1:10)
  tmp_mapping_df <- data.frame(
    probe = paste0("Probe_", c(1,1,2,2,2,3,4,5,6,7)),
    gene_symbol = paste0("Gene_", c(1,2,1,3,4,5,3,1,6,3))
  )

  expect_snapshot(expression_compression(tmp_expr_data, tmp_mapping_df,
                               "highest_mean", "log",
                               verbose = FALSE))

})


test_that("expression_compression pca", {

  tmp_expr_data <- data.frame(
    ID_1 = 1:10,
    ID_2 = c(2,3,6,7,10,1000,1,1,1,1),
    ID_3 = rev(c(2,3,6,7,10,1000,1,1,1,1))
  )
  rownames(tmp_expr_data) <- paste0("Probe_", 1:10)
  tmp_mapping_df <- data.frame(
    probe = paste0("Probe_", c(1,1,2,2,2,3,4,5,6,7)),
    gene_symbol = paste0("Gene_", c(1,2,1,3,4,5,3,1,6,3))
  )

  expect_snapshot(expression_compression(tmp_expr_data, tmp_mapping_df,
                               compress_fun = "pc1"))
})

test_that("expression_compression no dups", {

  tmp_expr_data <- data.frame(
    ID_1 = 1:10,
    ID_2 = 10:1,
    ID_3 = c(2,3,6,7,10,1000,1,1,1,1)
  )
  rownames(tmp_expr_data) <- paste0("Probe_", 0:9)
  tmp_mapping_df <- data.frame(
    probe = paste0("Probe_", 0:9),
    gene_symbol = paste0("Gene_", 0:9)
  )

  expect_snapshot(expression_compression(tmp_expr_data, tmp_mapping_df))
})



test_that("expression_compression errors", {
  tmp_expr_data <- data.frame(
    ID_1 = 1:10,
    ID_2 = 10:1,
    ID_3 = c(0,3,6,7,10,1000,1,1,1,1)
  )
  rownames(tmp_expr_data) <- paste0("Probe_", 0:9)
  tmp_mapping_df <- data.frame(
    probe = paste0("Probe_", 0:9),
    gene_symbol = paste0("Gene_", 0:9)
  )

  expect_error(
    expression_compression(tmp_expr_data, tmp_mapping_df, "mean", "log",
                 verbose = FALSE),
    "Can't do log transformation with exprs_data values <= 0"
  )
  tmp_mapping_df <- data.frame(
    probe = paste0("Probe_", 10:19),
    gene_symbol = paste0("Gene_", 0:9)
  )
  expect_error(
    expression_compression(tmp_expr_data, tmp_mapping_df, verbose = FALSE),
    "Could not link rownames\\(exprs_data\\) to mapping_df\\[,1\\]\\!"
  )

})
