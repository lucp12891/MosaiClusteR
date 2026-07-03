test_that("plot/feature functions are exported functions", {
  expect_true(is.function(BinFeaturesPlot_SingleData))
  expect_true(is.function(BinFeaturesPlot_MultipleData))
  expect_true(is.function(ContFeaturesPlot))
  expect_true(is.function(PlotPathways))
  expect_true(is.function(ProfilePlot))
})

test_that("BinFeaturesPlot_MultipleData needs plotrix", {
  testthat::skip_if_not_installed("plotrix")
  expect_true(is.function(BinFeaturesPlot_MultipleData))
})

test_that("PlotPathways needs MLP/biomaRt/org.Hs.eg.db", {
  testthat::skip_if_not_installed("MLP")
  testthat::skip_if_not_installed("biomaRt")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  expect_true(is.function(PlotPathways))
})

test_that("ProfilePlot needs Biobase for ExpressionSet input", {
  testthat::skip_if_not_installed("Biobase")
  expect_true(is.function(ProfilePlot))
})
