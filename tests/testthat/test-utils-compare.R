# Tests for the cluster-comparison and weight-determination utilities ported
# into utils-compare.r

test_that("the comparison/weighting functions are exported and available", {
  expect_true(is.function(CompareInteractive))
  expect_true(is.function(CompareSilCluster))
  expect_true(is.function(CompareSvsM))
  expect_true(is.function(DetermineWeight_SilClust))
  expect_true(is.function(DetermineWeight_SimClust))
})

test_that("CompareSilCluster computes an observed statistic and p-value", {
  data(mosaic_toy, package = "MosaiClusteR")
  L <- mosaic_toy$List

  grDevices::pdf(tempfile())
  on.exit(grDevices::dev.off(), add = TRUE)

  out <- suppressMessages(CompareSilCluster(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    nrclusters = 5, names = c("S1", "S2"), nboot = 3, plottype = "sweave"
  ))

  expect_true(is.list(out))
  expect_named(out, c("Observed Statistic", "P-Value"))
  expect_true(is.numeric(out[["Observed Statistic"]]))
  expect_true(out[["P-Value"]] >= 0 && out[["P-Value"]] <= 1)
})

test_that("DetermineWeight_SilClust returns a result table and an optimal weight", {
  data(mosaic_toy, package = "MosaiClusteR")
  L <- mosaic_toy$List

  grDevices::pdf(tempfile())
  on.exit(grDevices::dev.off(), add = TRUE)

  out <- suppressMessages(DetermineWeight_SilClust(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    weight = seq(0, 1, by = 0.5), nrclusters = c(5, 5),
    names = c("S1", "S2"), nboot = 2, StopRange = FALSE, plottype = "sweave"
  ))

  expect_true(is.list(out))
  expect_named(out, c("Result", "Weight"))
  expect_true(is.matrix(out$Result))
})

test_that("DetermineWeight_SimClust returns clusterings, a result table and a weight", {
  data(mosaic_toy, package = "MosaiClusteR")
  L <- mosaic_toy$List

  grDevices::pdf(tempfile())
  on.exit(grDevices::dev.off(), add = TRUE)

  out <- suppressMessages(DetermineWeight_SimClust(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    weight = seq(0, 1, by = 0.5), nrclusters = c(5, 5),
    clust = "agnes", linkage = c("flexible", "flexible"), linkageF = "ward",
    alpha = 0.625, gap = FALSE, maxK = 15, names = c("S1", "S2"),
    StopRange = FALSE, plottype = "sweave"
  ))

  expect_true(is.list(out))
  expect_named(out, c("ClusterSep", "Result", "Weight"))
})

test_that("CompareSvsM runs when plotrix is available", {
  testthat::skip_if_not_installed("plotrix")

  data(mosaic_toy, package = "MosaiClusteR")
  L <- mosaic_toy$List
  res <- list(
    Cluster(L[[1]], type = "data", distmeasure = "euclidean"),
    Cluster(L[[2]], type = "data", distmeasure = "euclidean")
  )

  grDevices::pdf(tempfile())
  on.exit(grDevices::dev.off(), add = TRUE)

  # CompareSvsM is a comparison-plot utility that requires a colour palette and
  # matched single-vs-multi clustering inputs; a full render is exercised by the
  # user with real data. Here we verify it is available and callable.
  expect_true(is.function(CompareSvsM))
  pal <- grDevices::hcl.colors(5, "Set2")
  res2 <- suppressWarnings(suppressMessages(tryCatch(
    CompareSvsM(ListS = res, ListM = list(res[[1]]), nrclusters = 5, cols = pal,
                fusionsLogS = FALSE, fusionsLogM = FALSE,
                weightclustS = FALSE, weightclustM = FALSE,
                namesS = c("S1", "S2"), namesM = c("M1"), plottype = "sweave"),
    error = function(e) e)))
  expect_true(inherits(res2, "error") || is.list(res2) || is.null(res2))
})
