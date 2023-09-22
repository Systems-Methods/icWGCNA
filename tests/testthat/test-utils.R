test_that("calcEigenGene one row", {
  expect_equal(
    calcEigenGene(matrix(1:100, nrow = 1)),
    matrix(1:100, nrow = 1)
  )
})

test_that("dropModuels logFlag and Kurts", {
  expect_snapshot(
    dropModuels(testing_data[2:10,1:10],
                Kurts = 1:9,
                corCut = .01,
                logFlag = TRUE)
  )
})

test_that("dropModuels large corCut", {
  expect_snapshot(
    dropModuels(testing_data[2:10,1:10],
                corCut = .01,
                logFlag = TRUE)
  )
})

