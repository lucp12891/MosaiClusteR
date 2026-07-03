test_that("LinkBasedClustering CTS (pure R) runs on mosaic_toy", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  res <- LinkBasedClustering(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE),
    linkBasedMethod = "cts", decayfactor = 0.8, niter = 5,
    linkBasedLinkage = "ward", nrclusters = c(3, 3)
  )

  expect_type(res, "list")
  expect_equal(attr(res, "method"), "LinkBased")
  expect_true(all(c("DistM", "Clust") %in% names(res)))
  expect_equal(dim(res$DistM), c(nrow(L[[1]]), nrow(L[[1]])))
  expect_s3_class(res$Clust, "agnes")

  labs <- stats::cutree(res$Clust, k = 3)
  expect_length(labs, nrow(L[[1]]))
  expect_equal(length(unique(labs)), 3)
})

test_that("LinkBasedClustering SRS (pure R) runs on mosaic_toy", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  res <- LinkBasedClustering(
    List = L, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE),
    linkBasedMethod = "srs", decayfactor = 0.8, niter = 5,
    linkBasedLinkage = "ward", nrclusters = c(3, 3)
  )

  expect_equal(attr(res, "method"), "LinkBased")
  expect_equal(dim(res$DistM), c(nrow(L[[1]]), nrow(L[[1]])))
  expect_s3_class(res$Clust, "agnes")
})

test_that("LinkBasedClustering ASRS requires an executable", {
  data(mosaic_toy)
  L <- mosaic_toy$List

  expect_error(
    LinkBasedClustering(
      List = L, type = "data",
      distmeasure = c("euclidean", "euclidean"),
      linkBasedMethod = "asrs", nrclusters = c(3, 3),
      executable = FALSE
    ),
    "asrs"
  )
})
