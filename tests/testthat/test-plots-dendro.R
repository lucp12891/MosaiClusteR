# Tests for the dendrogram / comparison plotting utilities ported from the
# original MoSaiClusteR project: ComparePlot, Cyclogram, ColorsNames,
# ClusterCols, ClusterPlot, ColorPalette, LabelCols and LabelPlot.

test_that("plotting helpers are exported functions", {
  expect_true(is.function(ComparePlot))
  expect_true(is.function(Cyclogram))
  expect_true(is.function(ColorsNames))
  expect_true(is.function(ClusterCols))
  expect_true(is.function(ClusterPlot))
  expect_true(is.function(ColorPalette))
  expect_true(is.function(LabelCols))
  expect_true(is.function(LabelPlot))
})

test_that("ColorPalette returns the requested number of hex colours", {
  cols <- ColorPalette(c("red", "green", "blue"), ncols = 6)
  expect_length(cols, 6)
  expect_true(all(grepl("^#", cols)))
})

test_that("ColorsNames maps a reorder matrix to colour names", {
  M <- matrix(c(1, 1, 2, 2,
                1, 2, 1, 2), nrow = 2, byrow = TRUE)
  pal <- ColorPalette(c("red", "green"), ncols = 2)
  Names <- ColorsNames(M, pal)
  expect_equal(dim(Names), dim(M))
  expect_false(any(is.na(Names)))
})

test_that("ClusterPlot and LabelPlot render without error on the toy data", {
  data(mosaic_toy)
  cl <- Cluster(mosaic_toy$List[[1]], type = "data",
                distmeasure = "euclidean", linkage = "ward")

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  cols <- ColorPalette(c("red", "green", "blue"), ncols = 3)
  expect_error(
    ClusterPlot(cl, nrclusters = 3, cols = cols, plottype = "sweave"),
    NA)

  sel <- rownames(mosaic_toy$List[[1]])[1:5]
  expect_error(
    LabelPlot(cl, sel1 = sel, col1 = "darkorchid"),
    NA)
})

test_that("ComparePlot produces a rectangular comparison on the toy data", {
  testthat::skip_if_not_installed("plotrix")
  data(mosaic_toy)
  c1 <- Cluster(mosaic_toy$List[[1]], type = "data",
                distmeasure = "euclidean", linkage = "ward")
  c2 <- Cluster(mosaic_toy$List[[2]], type = "data",
                distmeasure = "euclidean", linkage = "ward")
  L <- list(c1, c2)
  cols <- ColorPalette(c("red", "green", "blue", "orange"), ncols = 8)

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(
    ComparePlot(L, nrclusters = 3, cols = cols, names = c("S1", "S2"),
                plottype = "sweave"),
    NA)
})

test_that("Cyclogram requires the circlize package", {
  testthat::skip_if_not_installed("circlize")
  succeed()
})
