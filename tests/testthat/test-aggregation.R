test_that("HBGF (spectral) returns a partition over all objects", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  out <- HBGF(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    clust = "agnes", linkage = c("ward", "ward"),
    nrclusters = c(7, 7), graphPartitioning = "Spec", optimalk = 7
  )

  expect_null(out$DistM)
  expect_equal(attr(out, "method"), "Ensemble")
  cl <- out$Clust$Clusters
  expect_length(cl, nrow(L[[1]]))
  expect_equal(names(cl), rownames(L[[1]]))
  expect_true(length(unique(cl)) >= 1)
})

test_that("ClusteringAggregation Aggl yields a named partition", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  out <- ClusteringAggregation(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    clust = "agnes", linkage = c("ward", "ward"), alpha = 0.625,
    nrclusters = c(7, 7), agglMethod = "Aggl",
    improve = TRUE, distThresh_B = 0.5, distThresh_A = 0.8
  )

  expect_null(out$DistM)
  expect_equal(attr(out, "method"), "Ensemble")
  cl <- out$Clust$Clusters
  expect_length(cl, nrow(L[[1]]))
  expect_setequal(names(cl), rownames(L[[1]]))
})

test_that("ClusteringAggregation Balls yields a named partition", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  out <- ClusteringAggregation(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    clust = "agnes", linkage = c("ward", "ward"), alpha = 0.625,
    nrclusters = c(7, 7), agglMethod = "Balls",
    improve = FALSE, distThresh_B = 0.5, distThresh_A = 0.8
  )

  cl <- out$Clust$Clusters
  expect_length(cl, nrow(L[[1]]))
  expect_setequal(names(cl), rownames(L[[1]]))
})

test_that("ClusteringAggregation Furthest yields a named partition", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  out <- ClusteringAggregation(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE), method = c(NULL, NULL),
    clust = "agnes", linkage = c("ward", "ward"), alpha = 0.625,
    nrclusters = c(7, 7), agglMethod = "Furthest",
    improve = FALSE, distThresh_B = 0.5, distThresh_A = 0.8
  )

  cl <- out$Clust$Clusters
  expect_length(cl, nrow(L[[1]]))
  expect_setequal(names(cl), rownames(L[[1]]))
})
