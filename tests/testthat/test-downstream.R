test_that("panglaoDB enrichment input checking", {

  expect_error(
    compute_panglaoDB_enrichment(testing_results, pangDB = testing_pangDB),
    "t_memb must be a martix or data.frame"
  )
  expect_error(
    compute_panglaoDB_enrichment(testing_results$community_signature,
                                 pangDB = testing_pangDB),
    "t_memb values can't be <-1 or >1"
  )
  bad_data <- testing_results$community_membership
  rownames(bad_data) <- 1:nrow(bad_data)
  expect_error(
    compute_panglaoDB_enrichment(bad_data,
                                 pangDB = testing_pangDB),
    'No rownames of "t_memb" are in the provided "prolif" list'
  )

  expect_error(
    compute_panglaoDB_enrichment(testing_results$community_membership,
                                 pangDB = data.frame(aaa = 1:10,bbb = 11:20)),
    'expecting pangDB variables "cell.type" and "official.gene.symbol" \\(after making syntactically valid names using make.names\\(\\) function\\)'
  )

})
