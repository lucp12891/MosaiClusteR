test_that("CVAA runs on mosaic_toy (CVAA, optimalk == nrclustersR)", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  # Build a reference partition from the first data source.
  Ref <- Cluster(L[[1]], type = "data", distmeasure = "euclidean",
                 normalize = FALSE, method = NULL, clust = "agnes",
                 linkage = "ward", alpha = 0.625)
  attr(Ref, "method") <- "Single"

  res <- CVAA(Reference = Ref, nrclustersR = 5, List = L, typeL = "data",
              distmeasure = c("euclidean", "euclidean"),
              normalize = c(FALSE, FALSE), method = c(NULL, NULL),
              clust = "agnes", linkage = c("ward", "ward"), alpha = 0.625,
              nrclusters = c(5, 5), gap = FALSE, maxK = 15,
              votingMethod = "CVAA", optimalk = 5)

  expect_type(res, "list")
  expect_equal(attr(res, "method"), "Ensemble")
  expect_true(!is.null(res$Clust$Clusters))
  expect_equal(length(res$Clust$Clusters), nrow(L[[1]]))
})

test_that("CVAA runs on mosaic_toy (W-CVAA)", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  Ref <- Cluster(L[[1]], type = "data", distmeasure = "euclidean",
                 linkage = "ward")
  attr(Ref, "method") <- "Single"

  res <- CVAA(Reference = Ref, nrclustersR = 5, List = L, typeL = "data",
              distmeasure = c("euclidean", "euclidean"),
              linkage = c("ward", "ward"), nrclusters = c(5, 5),
              votingMethod = "W-CVAA", optimalk = 5)

  expect_equal(attr(res, "method"), "Ensemble")
  expect_equal(length(res$Clust$Clusters), nrow(L[[1]]))
})

test_that("EvidenceAccumulation runs on mosaic_toy (SL_agnes, pure R)", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  res <- EvidenceAccumulation(List = L, type = "data",
                              distmeasure = c("euclidean", "euclidean"),
                              normalize = c(FALSE, FALSE), method = c(NULL, NULL),
                              clust = "agnes", linkage = c("ward", "ward"),
                              alpha = 0.625, nrclusters = c(5, 5),
                              gap = FALSE, maxK = 15,
                              graphPartitioning = "SL_agnes")

  expect_type(res, "list")
  expect_equal(attr(res, "method"), "Single Clustering")
  expect_equal(length(res$Clust$Clusters), nrow(L[[1]]))
})

test_that("EvidenceAccumulation runs on mosaic_toy (SL threshold path)", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  res <- EvidenceAccumulation(List = L, type = "data",
                              distmeasure = c("euclidean", "euclidean"),
                              nrclusters = c(5, 5),
                              graphPartitioning = "SL", t = 0.5)

  expect_equal(attr(res, "method"), "Ensemble")
  expect_equal(length(res$Clust$Clusters), nrow(L[[1]]))
})

test_that("EvidenceAccumulation MST path requires igraph", {
  skip_if_not_installed("igraph")
  data(mosaic_toy)
  L <- mosaic_toy$List

  res <- EvidenceAccumulation(List = L, type = "data",
                              distmeasure = c("euclidean", "euclidean"),
                              nrclusters = c(5, 5),
                              graphPartitioning = "MST")
  expect_equal(attr(res, "method"), "Ensemble")
})
