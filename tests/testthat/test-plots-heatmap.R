test_that("BoxPlotDistance is exported and callable", {
  expect_true(is.function(BoxPlotDistance))
})

test_that("distanceheatmaps internal helper exists", {
  expect_true(is.function(MosaiClusteR:::distanceheatmaps))
})

test_that("HeatmapPlot is exported and callable", {
  expect_true(is.function(HeatmapPlot))
})

test_that("HeatmapSelection is exported and callable", {
  expect_true(is.function(HeatmapSelection))
})

test_that("SimilarityHeatmap is exported and callable", {
  expect_true(is.function(SimilarityHeatmap))
})

test_that("BoxPlotDistance requires its suggested plotting packages", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("gridExtra")
  data(mosaic_toy)
  grDevices::pdf(tempfile())
  on.exit(grDevices::dev.off(), add = TRUE)
  # type="dist" path with two distance matrices avoids interactive device use
  D1 <- Distance(mosaic_toy$List[[1]], distmeasure = "euclidean")
  D2 <- Distance(mosaic_toy$List[[2]], distmeasure = "euclidean")
  expect_error(suppressMessages(suppressWarnings(
    BoxPlotDistance(D1, D2, type = "dist", lab1 = "S1", lab2 = "S2",
                    limits1 = c(0.3, 0.7), plot = 1,
                    StopRange = FALSE, plottype = "sweave")
  )), NA)
})

test_that("SimilarityHeatmap draws on a clustering result", {
  testthat::skip_if_not_installed("gplots")
  data(mosaic_toy)
  grDevices::pdf(tempfile())
  on.exit(grDevices::dev.off(), add = TRUE)
  cl <- Cluster(mosaic_toy$List[[1]], type = "data",
                distmeasure = "euclidean", normalize = FALSE,
                clust = "agnes", linkage = "ward")
  expect_error(suppressMessages(suppressWarnings(
    SimilarityHeatmap(cl, type = "clust", plottype = "sweave")
  )), NA)
})
