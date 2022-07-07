test_that("static test data", {

  testing_data <- readRDS(test_path("fixtures", "testing_data.rds"))
  testing_results <- readRDS(test_path("fixtures", "testing_results.rds"))

  results_plus <- purrr::quietly(
    ~ icwgcna(testing_data, maxIt = 3,covCut = .66, mat_mult_method = 'RcppEigen')
  )()

  expect_equal(results_plus$result, testing_results)
  expect_equal(
    results_plus$messages,
    c("Computing 1263 x 1263 TOM distance for subset of genes with higher variance\n",
      "number of modules found is 11\n",
      "eigegenes trimmed to 11 due to correlation > 0.8 max eigenCor = 0.77\n",
      "0.201680.084970.06978\n",
      "Done with iteration: 1 : current number of gene communities is 11 \n\n\n",
      "Computing 913 x 913 TOM distance for subset of genes with higher variance\n",
      "number of modules found is 4\n",
      "eigegenes trimmed to 4 due to correlation > 0.8 max eigenCor = 0.15\n",
      "0.115470.102250.06965\n",
      "eigegenes trimmed to 13 due to correlation > 0.8 max eigenCor = 0.77\n",
      "Done with iteration: 2 : current number of gene communities is 13 \n\n\n",
      "Computing 913 x 913 TOM distance for subset of genes with higher variance\n",
      "number of modules found is 10\n",
      "eigegenes trimmed to 10 due to correlation > 0.8 max eigenCor = 0.7\n",
      "0.121860.077380.07039\n",
      "eigegenes trimmed to 18 due to correlation > 0.8 max eigenCor = 0.73\n",
      "Done with iteration: 3 : current number of gene communities is 18 \n\n\n",
      "Reached maximimum number of iterations\n"
    )
  )
  expect_equal(results_plus$output, '')
  expect_equal(results_plus$warnings, character(0))

})


test_that("input checking", {
  # ex must be numeric
  expect_error(
    icwgcna(cbind(paste('a', testing_data[,1]), testing_data)),
    "all 'ex' columns must be numeric"
  )


})
