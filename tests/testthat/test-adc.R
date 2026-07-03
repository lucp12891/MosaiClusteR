test_that("ADC returns expected structure", {
  data(mosaic_toy)
  L <- mosaic_toy$List
  res <- ADC(List = L, distmeasure = "euclidean", normalize = FALSE,
             method = NULL, clust = "agnes", linkage = "flexible", alpha = 0.625)

  expect_type(res, "list")
  expect_named(res, c("AllData", "DistM", "Clust"))
  expect_equal(attr(res, "method"), "ADC")

  # Fused data has objects in rows and all features stacked in columns
  expect_equal(nrow(res$AllData), nrow(L[[1]]))
  expect_equal(ncol(res$AllData), ncol(L[[1]]) + ncol(L[[2]]))

  # Distance matrix is square over the objects
  expect_equal(dim(res$DistM), c(nrow(L[[1]]), nrow(L[[1]])))

  # Clustering result is an agnes object
  expect_s3_class(res$Clust, "agnes")
})

test_that("ADC errors on non-list input", {
  expect_error(ADC(List = matrix(1:4, 2, 2)), "list")
})

test_that("ADEC with random sampling and fixed nrclusters runs", {
  data(mosaic_toy)
  L <- mosaic_toy$List
  res <- ADEC(List = L, distmeasure = "euclidean", normalize = FALSE,
              method = NULL, t = 3, r = 50, nrclusters = 5,
              clust = "agnes", linkage = "flexible", alpha = 0.625)

  expect_type(res, "list")
  expect_named(res, c("AllData", "DistM", "Clust"))
  expect_equal(attr(res, "method"), "ADEC")

  # Incidence/co-association matrix is square over the objects
  expect_equal(dim(res$DistM), c(nrow(L[[1]]), nrow(L[[1]])))
  expect_s3_class(res$Clust, "agnes")
})

test_that("ADEC with sequence of nrclusters runs", {
  data(mosaic_toy)
  L <- mosaic_toy$List
  res <- ADEC(List = L, distmeasure = "euclidean", normalize = FALSE,
              method = NULL, t = 5, r = NULL, nrclusters = seq(2, 5, 1),
              clust = "agnes", linkage = "flexible", alpha = 0.625)

  expect_equal(attr(res, "method"), "ADEC")
  expect_equal(dim(res$DistM), c(nrow(L[[1]]), nrow(L[[1]])))
  expect_s3_class(res$Clust, "agnes")
})

test_that("ADEC errors when nrclusters is NULL", {
  data(mosaic_toy)
  expect_error(ADEC(List = mosaic_toy$List, nrclusters = NULL),
               "number of cluters")
})
