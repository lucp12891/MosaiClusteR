test_that("M_ABCpp returns a standard MosaiClusteR result", {
  data(mosaic_toy)
  fit <- M_ABCpp(mosaic_toy$List, numsim = 40, NC = 3,
                 distmeasure = c("euclidean", "euclidean"))
  expect_true(all(c("DistM", "Clust") %in% names(fit)))
  expect_true(inherits(fit$Clust, c("agnes", "twins")))
  expect_equal(dim(fit$DistM), c(60, 60))
  expect_true(isSymmetric(unname(fit$DistM)))
  expect_equal(unname(diag(fit$DistM)), rep(0, 60))
})

test_that("M_ABCpp recovers structure with var and data-nugget weighting", {
  data(mosaic_toy)
  truth <- mosaic_toy$truth
  set.seed(11)
  fit_v <- M_ABCpp(mosaic_toy$List, numsim = 60, NC = 3,
                   weighting = c("var", "var"),
                   distmeasure = c("euclidean", "euclidean"))
  fit_n <- M_ABCpp(mosaic_toy$List, numsim = 60, NC = 3,
                   weighting = c("nugget", "nugget"),
                   distmeasure = c("euclidean", "euclidean"),
                   nugget_args = list(max_nuggets = 25, seed = 1))
  expect_gt(cluster_agreement(mosaic_labels(fit_v, 3), truth)[["ARI"]], 0.6)
  expect_gt(cluster_agreement(mosaic_labels(fit_n, 3), truth)[["ARI"]], 0.6)
})

test_that("ABCpp.SingleInMultiple runs every weighting scheme", {
  data(mosaic_toy)
  for (w in c("var", "cv", "nugget", "equal")) {
    lab <- ABCpp.SingleInMultiple(mosaic_toy$List[[1]], numsim = 20, NC = 3,
                                  weighting = w,
                                  nugget_args = list(max_nuggets = 20, seed = 1))
    expect_equal(dim(lab), c(20, 60))
  }
})

test_that("M_ABCdist fuses sources into a valid dissimilarity", {
  data(mosaic_toy)
  fit <- M_ABCdist(mosaic_toy$List, numsim = 40,
                   weighting = c("var", "nugget"))
  expect_equal(dim(fit$DistM), c(60, 60))
  expect_true(all(fit$DistM >= 0))
  expect_length(fit$source_D, 2)
})

test_that("M_ABCdist.WC fuses via WeightedClust", {
  data(mosaic_toy)
  fit <- suppressWarnings(M_ABCdist.WC(mosaic_toy$List, numsim = 30,
                      wc_weight = seq(1, 0, -0.25), wc_weightclust = 0.5))
  expect_equal(dim(fit$DistM), c(60, 60))
  expect_true(inherits(fit$Clust, c("agnes", "twins")))
})

test_that("compiled C++ kernels produce the same co-clustering counts as pure R", {
  cc <- getFromNamespace("count_coclustering_cpp", "MosaiClusteR")
  set.seed(7)
  res <- matrix(sample(0:3, 8 * 6, replace = TRUE), nrow = 8)  # R x N labels
  out <- cc(res)
  # pure-R reference
  N <- ncol(res); co <- matrix(0, N, N); nc <- matrix(0, N, N)
  for (i in 1:N) for (j in i:N) {
    both <- res[, i] != 0 & res[, j] != 0
    co[i, j] <- co[j, i] <- sum(both)
    nc[i, j] <- nc[j, i] <- sum(both & res[, i] != res[, j])
  }
  expect_equal(out$co_sel, co)
  expect_equal(out$not_co_clust, nc)
})
