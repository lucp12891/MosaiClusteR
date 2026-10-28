# Regression tests for Cluster() normalisation forwarding and input orientation.

test_that("Cluster forwards the normalisation method to Distance", {
  set.seed(1)
  # a feature on a very different scale, so normalisation changes the distances
  X <- cbind(matrix(rnorm(20 * 4), 20, 4), big = rnorm(20, sd = 1000))

  d_none <- Cluster(X, distmeasure = "euclidean", normalize = FALSE)$DistM
  d_std  <- Cluster(X, distmeasure = "euclidean", normalize = "Standardize")$DistM

  # forwarding a method must actually change the distances
  expect_false(isTRUE(all.equal(d_none, d_std)))
  # and must match calling Distance() with that method directly
  d_direct <- as.matrix(Distance(X, distmeasure = "euclidean", normalize = "Standardize"))
  expect_equal(unname(d_std), unname(d_direct))
})

test_that("Cluster(normalize = TRUE) selects a method instead of erroring", {
  set.seed(2)
  X <- matrix(rnorm(15 * 3), 15, 3)
  # Previously normalize = TRUE reached Normalization(method = TRUE) and failed.
  expect_silent(fit <- Cluster(X, distmeasure = "euclidean", normalize = TRUE))
  # default method is "Standardize"
  expect_equal(
    unname(fit$DistM),
    unname(Cluster(X, normalize = "Standardize")$DistM)
  )
  # `method` is honoured when normalize = TRUE
  expect_equal(
    unname(Cluster(X, normalize = TRUE, method = "Range")$DistM),
    unname(Cluster(X, normalize = "Range")$DistM)
  )
})

test_that("Cluster treats objects as rows (no transpose needed)", {
  data(mosaic_toy)
  X <- mosaic_toy$List[[1]]          # objects in rows (n x p)
  fit <- Cluster(X, type = "data", distmeasure = "euclidean")
  expect_equal(nrow(fit$DistM), nrow(X))
  lab <- mosaic_labels(fit, k = 3)
  expect_length(lab, nrow(X))
  expect_equal(length(unique(lab)), 3L)
})

test_that("mosaic_labels validates k with a clear message", {
  data(mosaic_toy)
  fit <- Cluster(mosaic_toy$List[[1]], distmeasure = "euclidean")
  expect_error(mosaic_labels(fit, NA), "single positive integer")
  expect_error(mosaic_labels(fit, 2.7), "single positive integer")
  expect_error(mosaic_labels(fit, c(2, 3)), "single positive integer")
  expect_error(mosaic_labels(fit, 0), "single positive integer")
  expect_length(mosaic_labels(fit, 3), nrow(mosaic_toy$List[[1]]))
})
