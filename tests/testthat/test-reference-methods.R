# Tests for the reference-derived methods: spectral, NEMO, data-nugget
# clustering, intNMF and the LUCID wrapper.

test_that("spectral_clustering separates two well-separated blobs", {
  set.seed(1)
  X <- rbind(matrix(rnorm(25 * 4, 0), 25), matrix(rnorm(25 * 4, 8), 25))
  A <- exp(-as.matrix(dist(X))^2 / 4)
  sc <- spectral_clustering(A, k = 2)
  expect_length(sc$cluster, 50)
  expect_equal(cluster_agreement(sc$cluster, rep(1:2, each = 25))[["ARI"]], 1)
})

test_that("NEMO returns a fused affinity and recovers planted structure", {
  data(mosaic_toy)
  fit <- NEMO(mosaic_toy$List, k = 3, NN = 15)
  expect_true(all(c("FusedM", "DistM", "cluster", "Clust") %in% names(fit)))
  expect_equal(dim(fit$FusedM), c(60, 60))
  expect_gt(cluster_agreement(fit$cluster, mosaic_toy$truth)[["ARI"]], 0.6)
})

test_that("NEMO tolerates partial data (a source missing some objects)", {
  data(mosaic_toy)
  L <- mosaic_toy$List
  L[[2]] <- L[[2]][1:50, ]                 # source 2 covers only 50/60 objects
  fit <- NEMO(L, k = 3, NN = 12)
  expect_equal(nrow(fit$FusedM), 60)       # union of objects
})

test_that("Wkmeans minimises WWCSS and recovers blobs", {
  set.seed(2)
  x <- rbind(matrix(rnorm(30 * 3, 0), 30), matrix(rnorm(30 * 3, 6), 30))
  w <- runif(60, 1, 4)
  km <- Wkmeans(x, k = 2, weights = w)
  expect_equal(cluster_agreement(km$cluster, rep(1:2, each = 30))[["ARI"]], 1)
  expect_true(is.finite(km$wwcss))
})

test_that("Whclust returns a valid hclust over weighted points", {
  set.seed(3)
  x <- rbind(matrix(rnorm(8 * 3, 0), 8), matrix(rnorm(8 * 3, 7), 8))
  hc <- Whclust(x, weights = rep(3, 16))
  expect_s3_class(hc, "hclust")
  expect_equal(cluster_agreement(stats::cutree(hc, 2), rep(1:2, each = 8))[["ARI"]], 1)
})

test_that("nugget_cluster compresses then clusters back to all objects", {
  data(mosaic_toy)
  X <- do.call(cbind, mosaic_toy$List)     # fuse sources, objects in rows
  fit <- nugget_cluster(X, k = 3, max_nuggets = 30)
  expect_length(fit$cluster, 60)
  expect_gt(cluster_agreement(fit$cluster, mosaic_toy$truth)[["ARI"]], 0.6)
})

test_that("intNMF integrates sources and recovers structure", {
  data(mosaic_toy)
  fit <- intNMF(mosaic_toy$List, k = 3, nstart = 5, seed = 1)
  expect_equal(ncol(fit$W), 3)
  expect_equal(dim(fit$DistM), c(60, 60))
  expect_gt(cluster_agreement(fit$cluster, mosaic_toy$truth)[["ARI"]], 0.5)
})

test_that("LUCID is exported and guards the missing LUCIDus package", {
  expect_true(is.function(LUCID))
  skip_if(requireNamespace("LUCIDus", quietly = TRUE),
          "LUCIDus installed; guard not exercised")
  data(mosaic_toy)
  expect_error(
    LUCID(G = mosaic_toy$List[[2]][, 1:3], Z = mosaic_toy$List[[1]],
          Y = as.numeric(mosaic_toy$truth), K = 3),
    "LUCIDus")
})
