test_that("EnsembleClustering is an exported function", {
  expect_true(exists("EnsembleClustering"))
  expect_true(is.function(EnsembleClustering))
})

test_that("EHC is an exported function", {
  expect_true(exists("EHC"))
  expect_true(is.function(EHC))
})

test_that("EnsembleClustering errors when no executable is configured", {
  data(mosaic_toy)
  L <- mosaic_toy$List
  # executable = FALSE (default): no consensus is possible in pure R, so the
  # function must stop with a clear, informative message rather than crash.
  expect_error(
    EnsembleClustering(List = L, type = "data",
                       distmeasure = c("euclidean", "euclidean"),
                       normalize = c(FALSE, FALSE), method = c(NULL, NULL),
                       clust = "agnes", linkage = c("flexible", "flexible"),
                       nrclusters = c(7, 7), ensembleMethod = "CSPA",
                       executable = FALSE),
    regexp = "executable"
  )
})

test_that("EHC (METIS) errors when no executable is configured", {
  data(mosaic_toy)
  L <- mosaic_toy$List
  # graphPartitioning = "METIS" with executable = FALSE must stop informatively.
  expect_error(
    EHC(List = L, type = "data",
        distmeasure = c("euclidean", "euclidean"),
        normalize = c(FALSE, FALSE), method = c(NULL, NULL),
        clust = "agnes", linkage = c("flexible", "flexible"),
        graphPartitioning = "METIS", optimalk = 7, executable = FALSE),
    regexp = "executable"
  )
})

test_that("EnsembleClustering errors when a configured executable is missing", {
  data(mosaic_toy)
  L <- mosaic_toy$List
  expect_error(
    EnsembleClustering(List = L, type = "data",
                       distmeasure = c("euclidean", "euclidean"),
                       normalize = c(FALSE, FALSE), method = c(NULL, NULL),
                       clust = "agnes", linkage = c("flexible", "flexible"),
                       nrclusters = c(7, 7), ensembleMethod = "CSPA",
                       executable = "this_executable_does_not_exist_12345"),
    regexp = "locate|executable"
  )
})
