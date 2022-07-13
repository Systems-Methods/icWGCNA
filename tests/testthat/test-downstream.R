test_that("compute eigengene matrix input checking", {

  # ex must be numeric
  expect_error(
    compute_eigengene_matrix(cbind(paste('a', testing_data[,1]), testing_data),
                             testing_results$community_membership),
    "all 'ex' columns must be numeric"
  )
  expect_error(
    compute_eigengene_matrix(cbind(testing_data[,1] - 1, testing_data),
                             testing_results$community_membership),
    "all values of ex must be >=0"
  )
  expect_warning(
    suppressMessages(
      compute_eigengene_matrix(
        cbind(testing_data[,1] + c(rep(100, 10),
                                   rep(0, nrow(testing_data) - 10)),
              testing_data),
      testing_results$community_membership)),
    "some values of ex are >100, strongly indicating ex is not in log space"
  )

  bad_data <- testing_data
  rownames(bad_data) <- 1:nrow(bad_data)
  expect_error(
    compute_eigengene_matrix(bad_data,
                             testing_results$community_membership),
    "No matching rownames in ex and membership_matrix"
  )

})

test_that("compute eigengene matrix status results", {
  results_plus <- purrr::quietly(
    ~ compute_eigengene_matrix(testing_data,testing_results$community_membership)
  )()

  expect_equal(results_plus$result,testing_eigengene_matrix)
  expect_equal(
    results_plus$messages,
    "Removing 2 genes with a 0 standard deviation\n"
  )
  expect_equal(results_plus$output, '')
  expect_equal(results_plus$warnings, character(0))

})






test_that("panglaoDB enrichment input checking", {

  expect_error(
    compute_panglaoDB_enrichment(testing_results, pangDB = testing_pangDB),
    "membership_matrix must be a martix or data.frame"
  )
  expect_error(
    compute_panglaoDB_enrichment(testing_results$community_signature,
                                 pangDB = testing_pangDB),
    "membership_matrix values can't be <-1 or >1"
  )

  expect_error(
    compute_panglaoDB_enrichment(testing_results$community_membership,
                                 pangDB = data.frame(aaa = 1:10,bbb = 11:20)),
    'expecting pangDB variables "cell.type" and "official.gene.symbol" \\(after making syntactically valid names using make.names\\(\\) function\\)'
  )

})

test_that("panglaoDB enrichment status results", {
  expect_equal(
    compute_panglaoDB_enrichment(testing_results$community_membership,
                                 pangDB = testing_pangDB),
    testing_enrichment
  )
})
