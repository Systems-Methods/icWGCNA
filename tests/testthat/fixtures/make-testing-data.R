# code to make testing data used
luad <- suppressMessages(
  UCSCXenaTools::getTCGAdata(project = "LUAD",
                             mRNASeq = TRUE,
                             mRNASeqType = "normalized",
                             clinical = FALSE,
                             download = TRUE)
  )
ex <- data.table::fread(luad$destfiles, data.table = FALSE)
mu <- apply(as.matrix(ex[,-1]),1,mean)
SD <- apply(as.matrix(ex[,-1]),1,sd)
n_picked <- 50
testing_data <- withr::with_seed(
  seed = 148654315,
  code = ex[mu > 2.5 & SD > 1.5, sample.int(ncol(ex), n_picked)]
)
# Adding rows with 0 SD, witch should get removed
testing_data <- rbind(rep(1,  ncol(testing_data)),
                      testing_data,
                      rep(1,  ncol(testing_data)))
saveRDS(testing_data, file = testthat::test_path('fixtures','testing_data.rds'))

# saving results file
results <- icwgcna(testing_data, maxIt = 3,covCut = .66, mat_mult_method = 'RcppEigen')
saveRDS(results, file = testthat::test_path('fixtures','testing_results.rds'))

