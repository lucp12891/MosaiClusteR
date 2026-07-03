test_that("Distance is symmetric with zero diagonal", {
  set.seed(1)
  m <- matrix(rnorm(8 * 5), nrow = 8)
  D <- Distance(m, "euclidean")
  expect_equal(dim(D), c(8, 8))
  expect_true(isSymmetric(unname(round(D, 10))))
  expect_equal(unname(diag(D)), rep(0, 8))

  b <- matrix(rbinom(8 * 6, 1, 0.5), nrow = 8)
  Dt <- Distance(b, "tanimoto")
  expect_true(all(Dt >= 0 & Dt <= 1))
})

test_that("Normalization range maps to [0,1]", {
  m <- matrix(rnorm(50), 10, 5)
  r <- Normalization(m, "Range")
  expect_gte(min(r), 0); expect_lte(max(r), 1)
  s <- Normalization(m, "S")
  expect_equal(unname(round(colMeans(s), 8)), rep(0, 5))
})

test_that("cluster_agreement is 1 for identical partitions up to relabelling", {
  a <- c(1, 1, 2, 2, 3, 3)
  b <- c(3, 3, 1, 1, 2, 2)
  ag <- cluster_agreement(a, b)
  expect_equal(unname(ag["ARI"]), 1)
  expect_equal(unname(ag["Jaccard"]), 1)
})

test_that("mosaic_labels extracts a partition of the right size", {
  data(mosaic_toy)
  fit <- M_ABCpp(mosaic_toy$List, numsim = 30, NC = 3)
  lab <- mosaic_labels(fit, k = 3)
  expect_length(lab, 60)
  expect_equal(length(unique(lab)), 3)
})

test_that("ABCdist.SingleInMultiple yields a scaled [0,1] dissimilarity", {
  data(mosaic_toy)
  out <- ABCdist.SingleInMultiple(mosaic_toy$List[[1]], numsim = 30)
  expect_true(all(out$D >= 0 & out$D <= 1))
  expect_equal(unname(diag(out$D)), rep(0, 60))
})
