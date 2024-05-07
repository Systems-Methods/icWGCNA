test_that("compute eigengene matrix input checking", {

  # ex must be numeric
  expect_error(
    compute_eigengene_matrix(
      cbind(paste("a", testing_data[, 1]), testing_data),
      testing_results$community_membership
    ),
    "all 'ex' columns must be numeric"
  )
  expect_error(
    compute_eigengene_matrix(
      cbind(testing_data[, 1] - 1, testing_data),
      testing_results$community_membership
    ),
    "all values of ex must be >=0"
  )
  expect_warning(
    suppressMessages(
      compute_eigengene_matrix(
        cbind(
          testing_data[, 1] + c(
            rep(100, 10),
            rep(0, nrow(testing_data) - 10)
          ),
          testing_data
        ),
        testing_results$community_membership
      )
    ),
    "some values of ex are >100, strongly indicating ex is not in log space"
  )

  bad_data <- testing_data
  rownames(bad_data) <- 1:nrow(bad_data)
  expect_error(
    compute_eigengene_matrix(
      bad_data,
      testing_results$community_membership
    ),
    "No matching rownames in ex and membership_matrix"
  )
})

test_that("compute eigengene matrix status results", {
  results_plus <- purrr::quietly(
    ~ compute_eigengene_matrix(testing_data, testing_results$community_membership)
  )()

  expect_equal(results_plus$result, testing_eigengene_matrix)
  expect_equal(
    results_plus$messages,
    "Removing 2 genes with a 0 standard deviation\n"
  )
  expect_equal(results_plus$output, "")
  expect_equal(results_plus$warnings, character(0))
})


test_that("compute eigengene matrix more options", {
  expect_snapshot(
    compute_eigengene_matrix(
      testing_data,
      testing_results$community_membership,
      pc_flag = FALSE
    )
  )
})

test_that("compute eigengene matrix even more options", {
  expect_snapshot(
    compute_eigengene_matrix(
      testing_data,
      testing_results$community_membership,
      cutoff = 0
    )
  )
})






test_that("panglaoDB enrichment input checking", {
  expect_error(
    compute_panglaoDB_enrichment(testing_results, pangDB = testing_pangDB),
    "membership_matrix must be a martix or data.frame"
  )
  expect_error(
    compute_panglaoDB_enrichment(testing_results$community_signature,
      pangDB = testing_pangDB
    ),
    "membership_matrix values can't be <-1 or >1"
  )

  expect_error(
    compute_panglaoDB_enrichment(testing_results$community_membership,
      pangDB = data.frame(aaa = 1:10, bbb = 11:20)
    ),
    'expecting pangDB variables "cell.type" and "official.gene.symbol" \\(after making syntactically valid names using make.names\\(\\) function\\)'
  )
})

test_that("panglaoDB enrichment status results", {
  expect_equal(
    compute_panglaoDB_enrichment(testing_results$community_membership,
      pangDB = testing_pangDB
    ),
    testing_panglaoDB_enrichment
  )
})






test_that("MSigDB enrichment input checking", {
  expect_error(
    compute_MSigDB_enrichment(testing_results),
    "membership_matrix must be a martix or data.frame"
  )
  expect_error(
    compute_MSigDB_enrichment(testing_results$community_signature),
    "membership_matrix values can't be <-1 or >1"
  )
  expect_error(
    compute_MSigDB_enrichment(testing_results$community_membership, cats = 'M'),
    'No "cats" found in MSigDB. Must use at least one of: C1, C2, C3, C4, C5, C6, C7, C8, H'
  )
})

test_that('requireNamespace stubbing (MSigDB)', {
  mockery::stub(compute_MSigDB_enrichment, 'requireNamespace', FALSE)

  expect_error(
    compute_MSigDB_enrichment(testing_results$community_membership),
    "Must have the following R packages installed for this function: msigdbr, foreach, tidyr")

})


test_that('parallel process issue', {
  mockery::stub(compute_MSigDB_enrichment, 'table',
                table(
                  rep(paste0('mA',1:10), 10),
                  rep(c("C3", "C6", "C7", "C8", "H"), 20))
                )

  expect_error(
    compute_MSigDB_enrichment(testing_results$community_membership, cats = 'H'),
    "problem computing enrichments, try a different parallel computing setup")

})

test_that("MSigDB enrichment status results", {
  results_plus <- purrr::quietly(
    ~ compute_MSigDB_enrichment(testing_results$community_membership,
                                cats = c("H", "AAAAAAA","C3", "C6", "C7", "C8",
                                         "BAD_CAT"))
  )()

  expect_equal(results_plus$result$top_enr, testing_MSigDB_enrichment$top_enr)
  expect_equal(results_plus$result$full_enr, testing_MSigDB_enrichment$full_enr)
  expect_equal(
    results_plus$messages,
    c("No parallel processing has been detected\n",
      "working on H\n",
      "working on C3\n",
      "working on C6\n",
      "working on C7\n",
      "working on C8\n"
    ))
  expect_equal(results_plus$output, "")
  expect_equal(results_plus$warnings,
               "The following \"cats\" are not in MSigDB: AAAAAAA, BAD_CAT")
})


test_that("MSigDB enrichment parallel", {
  cl <- parallel::makePSOCKcluster(2)
  doParallel::registerDoParallel(cl)

  results_plus <- purrr::quietly(
    ~ compute_MSigDB_enrichment(testing_results$community_membership)
  )()


  expect_equal(results_plus$result$top_enr, testing_MSigDB_enrichment$top_enr)
  expect_equal(results_plus$result$full_enr, testing_MSigDB_enrichment$full_enr)
  expect_equal(
    results_plus$messages,
    c("Using doParallelSNOW with 2 workers\n",
      "working on H\n",
      "working on C3\n",
      "working on C6\n",
      "working on C7\n",
      "working on C8\n"
    ))
  expect_equal(results_plus$output, "")
  expect_equal(results_plus$warnings, character(0))

  on.exit({
    try({
      doParallel::stopImplicitCluster()
      parallel::stopCluster(cl)
    })
  })
})




test_that("xCell enrichment input checking", {
  expect_error(
    compute_xCell_enrichment(testing_results),
    "membership_matrix must be a martix or data.frame"
  )
  expect_error(
    compute_xCell_enrichment(testing_results$community_signature),
    "membership_matrix values can't be <-1 or >1"
  )
})


test_that('requireNamespace stubbing (xCell)', {
  mockery::stub(compute_xCell_enrichment, 'requireNamespace', FALSE)

  expect_error(
    compute_xCell_enrichment(testing_results$community_membership),
    "Must have the following R packages installed for this function: xCell")

})


test_that("xCell enrichment status results", {
  expect_equal(
    compute_xCell_enrichment(testing_results$community_membership),
    testing_xCell_enrichment
  )
})






test_that("UMAP plotting input checking", {
  expect_error(
    make_network_umap(testing_results),
    "membership_matrix must be a martix or data.frame"
  )
  expect_error(
    make_network_umap(testing_results$community_signature),
    "membership_matrix values can't be <-1 or >1"
  )

  expect_error(
    make_network_umap(testing_results$community_membership,
                      community_memb_cut_main = 1,
                      community_memb_cut_secondary = 1),
    "Must have at least 2 communities after filtering. Try less restrictive cutoffs."
  )

  expect_error(
    make_network_umap(testing_results$community_membership,
                      gene_memb_cut_main = 1,
                      gene_memb_cut_secondary = 1),
    "Must have at least 2 genes after filtering. Try less restrictive cutoffs."
  )

  expect_error(
    make_network_umap(testing_results$community_membership, community_labels = 1),
    "community_labels must be a data.frame with 2 columns and a column named \"community\""
  )


})


test_that('requireNamespace stubbing (UMAP)', {
  mockery::stub(make_network_umap, 'requireNamespace', FALSE)

  expect_error(
    make_network_umap(testing_results$community_membership),
    "Must have the following R packages installed for this function: ggplot2, rlang, umap")

})

test_that("UMAP Success", {
  mockery::stub(make_network_umap, 'umap::umap', UMAP_testing_layout)

  custom_umap_specs <- umap::umap.defaults
  custom_umap_specs$random_state <- 94124456
  results_plus <- purrr::quietly(
    ~ make_network_umap(testing_results$community_membership,
                        umap_specs = custom_umap_specs,
                        community_labels = data.frame(community = 'mA1',
                                                      lab = 'Extra'))
  )()

  expect_equal(results_plus$result$layout,
               testing_UMAP_results$layout)
  expect_equal(results_plus$result$umap_w_legend[c("data",
                                                    "scales",
                                                    "coordinates",
                                                    "facet",
                                                    "labels")],
               testing_UMAP_results$umap_w_legend[c("data",
                                                    "scales",
                                                    "coordinates",
                                                    "facet",
                                                    "labels")],
               ignore_attr = TRUE, ignore_function_env = TRUE)
  expect_equal(results_plus$result$umap_w_annotation[c("data",
                                                       "scales",
                                                       "coordinates",
                                                       "facet",
                                                       "labels")],
               testing_UMAP_results$umap_w_annotation[c("data",
                                                        "scales",
                                                        "coordinates",
                                                        "facet",
                                                        "labels")],
               ignore_attr = TRUE, ignore_function_env = TRUE)
  expect_equal(
    results_plus$messages,
    c("Filtering from 18 communites to 15 communities for plotting.\n",
      "Then filtering from 2685 genes to 513 genes for plotting.\n"
    ))
  expect_equal(results_plus$output, "")
  expect_equal(results_plus$warnings, character(0))
})

test_that("find_unique_top_genes Success", {
  expect_snapshot(
    find_unique_top_genes(
      testing_results$community_membership
    )
  )
})

test_that("find_unique_top_genes high k", {
  expect_snapshot(
    find_unique_top_genes(
      testing_results$community_membership,
      K = 100
    ),
    error = TRUE
  )
})

test_that("find_unique_top_genes input checking", {
  expect_error(
    find_unique_top_genes(testing_results),
    "membership_matrix must be a martix or data.frame"
  )
  expect_error(
    find_unique_top_genes(testing_results$community_signature),
    "membership_matrix values can't be <-1 or >1"
  )
})


test_that('requireNamespace stubbing (UMAP)', {
  mockery::stub(map_eigengenes_on_seurat, 'requireNamespace', FALSE)

  expect_error(
    map_eigengenes_on_seurat(
      testing_Seurat,
      testing_results$community_membership
    ),
    "Must have the following R packages installed for this function: Seurat, UCell")

})


test_that("map_eigengenes_on_seurat run", {

  # Only meta.data is getting updated
  expect_snapshot(
    map_eigengenes_on_seurat(
      testing_Seurat,
      testing_results$community_membership
    )@meta.data
  )

})

test_that("map_eigengenes_on_seurat prefix and both method", {
  expect_snapshot(
    map_eigengenes_on_seurat(
      testing_Seurat,
      testing_results$community_membership,
      prefix = "Test",
      cutoff_method = 'both'
    )@meta.data
  )
})




