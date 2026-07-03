# Data-nugget weighting on the ABC engine (the .abc_zf internal + ABCpp).

test_that("internal .abc_zf returns a probability vector for every scheme", {
  zf <- getFromNamespace(".abc_zf", "MosaiClusteR")
  set.seed(1)
  X <- matrix(rnorm(40 * 15), nrow = 40)        # features x objects (internal)
  for (w in c("var", "cv", "nugget", "equal")) {
    p <- zf(X, w, nugget_args = list(max_nuggets = 8, seed = 1))
    expect_length(p, nrow(X))
    expect_equal(sum(p), 1, tolerance = 1e-8)
    expect_true(all(p >= 0))
  }
})

test_that("nugget weighting concentrates probability on informative features", {
  zf <- getFromNamespace(".abc_zf", "MosaiClusteR")
  set.seed(3)
  X <- matrix(rnorm(30 * 20), nrow = 30)        # 30 features x 20 objects
  X[1:5, 1:10] <- X[1:5, 1:10] + 5              # 5 informative features
  p <- zf(X, "nugget", nugget_args = list(max_nuggets = 6, seed = 3))
  expect_gt(mean(p[1:5]), mean(p[6:30]))
})

test_that("data-nugget weighting recovers structure through M_ABCpp", {
  data(mosaic_toy)
  fit <- M_ABCpp(mosaic_toy$List, numsim = 60, NC = 3,
                 weighting = c("nugget", "nugget"),
                 nugget_args = list(max_nuggets = 25, seed = 1))
  expect_gt(cluster_agreement(mosaic_labels(fit, 3), mosaic_toy$truth)[["ARI"]], 0.6)
})
