test_that("ConsensusClustering (IVC) returns a consensus partition", {
  data(mosaic_toy)
  set.seed(1)
  fit <- ConsensusClustering(
    List = mosaic_toy$List, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    clust = "agnes", linkage = c("ward", "ward"),
    nrclusters = c(3, 3), gap = FALSE,
    votingMethod = "IVC", optimalk = 3
  )
  expect_equal(attr(fit, "method"), "Ensemble")
  cl <- fit$Clust$Clusters
  expect_length(cl, nrow(mosaic_toy$List[[1]]))
  expect_true(all(cl %in% seq_len(3)))
})

test_that("ConsensusClustering (IPVC) returns a consensus partition", {
  data(mosaic_toy)
  set.seed(2)
  fit <- ConsensusClustering(
    List = mosaic_toy$List, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    clust = "agnes", linkage = c("ward", "ward"),
    nrclusters = c(3, 3), gap = FALSE,
    votingMethod = "IPVC", optimalk = 3
  )
  cl <- fit$Clust$Clusters
  expect_length(cl, nrow(mosaic_toy$List[[1]]))
  expect_named(cl)
})

test_that("WonM returns a consensus dissimilarity matrix and clustering", {
  data(mosaic_toy)
  set.seed(3)
  n <- nrow(mosaic_toy$List[[1]])
  fit <- WonM(
    List = mosaic_toy$List, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    nrclusters = list(seq(3, 6), seq(3, 6)),
    clust = "agnes", linkage = c("ward", "ward")
  )
  expect_equal(attr(fit, "method"), "WonM")
  expect_equal(dim(fit$DistM), c(n, n))
  expect_true(all(fit$DistM >= 0 & fit$DistM <= 1))
  expect_s3_class(fit$Clust, "agnes")
  lab <- stats::cutree(fit$Clust, k = 3)
  expect_length(lab, n)
})
