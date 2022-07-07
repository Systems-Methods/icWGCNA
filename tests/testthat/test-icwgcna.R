

test_that("input checking", {
  # ex must be numeric
  expect_error(
    icwgcna(cbind(paste('a', testing_data[,1]), testing_data)),
    "all 'ex' columns must be numeric"
  )

  expect_error(
    icwgcna(testing_data, maxIt = 26),
    "maxIt must be between 1 and 25"
  )
  expect_error(
    icwgcna(testing_data, maxIt = 0),
    "maxIt must be between 1 and 25"
  )

  expect_error(icwgcna(testing_data, expo = 0), "expo must be >0 and <=10, or NULL")
  expect_error(icwgcna(testing_data, expo = 11), "expo must be >0 and <=10, or NULL")

  expect_error(icwgcna(testing_data, q = 0), "q must be >0 and <1")
  expect_error(icwgcna(testing_data, q = 1), "q must be >0 and <1")

  expect_error(icwgcna(testing_data, corCut = 0), "corCut must be >0 and <1")
  expect_error(icwgcna(testing_data, corCut = 1), "corCut must be >0 and <1")

  expect_error(icwgcna(testing_data, covCut = 0), "covCut must be >0 and <1")
  expect_error(icwgcna(testing_data, covCut = 1), "covCut must be >0 and <1")



})


test_that("static test data", {

  results_plus <- purrr::quietly(
    ~ icwgcna(testing_data, maxIt = 3,covCut = .66, mat_mult_method = 'RcppEigen')
  )()

  expect_equal(results_plus$result, testing_results)
  expect_equal(
    results_plus$messages,
    c("Removing 2 genes with a 0 standard deviation\n",
      "Computing 1263 x 1263 TOM distance for subset of genes with higher variance\n",
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

  #checking rfast mult method
  results_plus_Rfast <- purrr::quietly(
    ~ icwgcna(testing_data, maxIt = 3,covCut = .66, mat_mult_method = 'Rfast')
  )()
  expect_equal(results_plus_Rfast, results_plus)


})



test_that("expo = NULL (angular dist)", {

  # testing angular distance results (expo = NULL)
  results_plus_expoNULL <- purrr::quietly(
    ~ icwgcna(testing_data, maxIt = 2, covCut = .66,
              mat_mult_method = 'RcppEigen', expo = NULL)
  )()
  expect_equal(
    results_plus_expoNULL$messages,
    c(
      "Removing 2 genes with a 0 standard deviation\n",
      "Computing 1263 x 1263 TOM distance for subset of genes with higher variance\n",
      "number of modules found is 10\n",
      "eigegenes trimmed to 10 due to correlation > 0.8 max eigenCor = 0.77\n",
      "0.201680.084970.06978\n",
      "Done with iteration: 1 : current number of gene communities is 10 \n\n\n",
      "Computing 913 x 913 TOM distance for subset of genes with higher variance\n",
      "number of modules found is 3\n",
      "eigegenes trimmed to 3 due to correlation > 0.8 max eigenCor = 0.27\n",
      "0.11760.102240.0697\n",
      "eigegenes trimmed to 12 due to correlation > 0.8 max eigenCor = 0.79\n",
      "Done with iteration: 2 : current number of gene communities is 12 \n\n\n",
      "Reached maximimum number of iterations\n"
    )
  )

})

test_that("spearman test", {

  results_plus_spearman <- purrr::quietly(
    ~ icwgcna(testing_data, maxIt = 2, covCut = .66,
              mat_mult_method = 'RcppEigen', Method = 'spearman')
  )()
  expect_equal(
    results_plus_spearman$messages,
    c(
      "Removing 2 genes with a 0 standard deviation\n",
      "Computing 1263 x 1263 TOM distance for subset of genes with higher variance\n",
      "number of modules found is 10\n",
      "eigegenes trimmed to 10 due to correlation > 0.8 max eigenCor = 0.75\n",
      "0.201680.084970.06978\n",
      "Done with iteration: 1 : current number of gene communities is 10 \n\n\n",
      "Computing 913 x 913 TOM distance for subset of genes with higher variance\n",
      "number of modules found is 4\n",
      "eigegenes trimmed to 3 due to correlation > 0.8 max eigenCor = 0\n",
      "0.115320.102510.06966\n",
      "eigegenes trimmed to 12 due to correlation > 0.8 max eigenCor = 0.75\n",
      "Done with iteration: 2 : current number of gene communities is 12 \n\n\n",
      "Reached maximimum number of iterations\n")
  )

})


