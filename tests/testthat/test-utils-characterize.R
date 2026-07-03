test_that("characterization and tracking utilities are exported functions", {
  expect_true(is.function(CharacteristicFeatures))
  expect_true(is.function(ChooseCluster))
  expect_true(is.function(FeatSelection))
  expect_true(is.function(FeaturesOfCluster))
  expect_true(is.function(TrackCluster))
})

test_that("ChooseCluster gene-expression path needs suggested packages", {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("a4Core")
  testthat::skip_if_not_installed("a4Base")
  expect_true(is.function(ChooseCluster))
})
