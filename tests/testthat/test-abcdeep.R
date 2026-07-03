test_that("ABCdeep.SingleInMultiple returns the ABCpp label contract", {
  data(mosaic_toy)
  lab <- ABCdeep.SingleInMultiple(mosaic_toy$List[[1]], numsim = 30, NC = 3,
                                  ae_epochs = 30, seed = 1)
  expect_s3_class(lab, "data.frame")
  expect_equal(nrow(lab), 30)
  expect_equal(ncol(lab), nrow(mosaic_toy$List[[1]]))
  expect_true(all(as.matrix(lab) >= 0))                 # 0 = not selected
  expect_true(max(as.matrix(lab)) <= 3)                 # NC base clusters
})

test_that("ABCdeep recovers strong simulated clusters", {
  data(mosaic_toy)
  set.seed(2)
  fit <- M_ABCdeep(mosaic_toy$List, numsim = 60, NC = 3, ae_epochs = 40, seed = 2)
  cut <- stats::cutree(stats::as.hclust(fit$Clust), 3)
  expect_gt(cluster_agreement(cut, mosaic_toy$truth)["ARI"], 0.6)
})

test_that("M_ABCdeep returns an M_ABC object usable by the helpers", {
  data(mosaic_toy)
  fit <- M_ABCdeep(mosaic_toy$List, numsim = 40, NC = 3, ae_epochs = 30, seed = 3)
  expect_s3_class(fit, "M_ABC")
  expect_true(all(c("DistM", "Clust") %in% names(fit)))
  expect_equal(dim(fit$DistM), c(nrow(mosaic_toy$List[[1]]), nrow(mosaic_toy$List[[1]])))
  expect_identical(attr(fit, "method"), "M-ABCdeep")
})

test_that("kmeans latent clustering also works and there is no fallback path", {
  data(mosaic_toy)
  lab <- ABCdeep.SingleInMultiple(mosaic_toy$List[[1]], numsim = 20, NC = 3,
                                  ae_epochs = 20, cluster_in_latent = "kmeans",
                                  seed = 4)
  expect_equal(nrow(lab), 20)
  expect_error(
    ABCdeep.SingleInMultiple(mosaic_toy$List[[1]], cluster_in_latent = "pca"),
    "ward|kmeans"
  )
})
