test_that("SelectnrClusters returns silhouette selection info without plotting", {
  data(mosaic_toy)
  sel <- SelectnrClusters(
    mosaic_toy$List, type = "data",
    distmeasure = c("euclidean", "euclidean"),
    normalize = c(FALSE, FALSE),
    nrclusters = seq(2, 5), plottype = "none"
  )

  expect_type(sel, "list")
  expect_named(sel, c("Silhoutte_Widths", "Optimal_Nr_of_CLusters"))

  # silhouette table: one row per candidate k, one column per source + average
  sw <- sel$Silhoutte_Widths
  expect_true(is.matrix(sw))
  expect_equal(nrow(sw), length(seq(2, 5)))
  expect_equal(ncol(sw), length(mosaic_toy$List) + 1)
  expect_true(all(is.finite(sw)))

  # chosen number of clusters is reported and lies in the requested range
  opt <- sel$Optimal_Nr_of_CLusters
  expect_s3_class(opt, "data.frame")
  expect_true(all(as.numeric(opt[1, ]) %in% seq(2, 5)))
})

test_that("SimilarityMeasure returns numeric similarities for a label matrix", {
  M <- rbind(
    c(1, 1, 2, 2, 3, 3),
    c(1, 1, 2, 3, 3, 3),
    c(2, 2, 1, 1, 3, 3)
  )
  sim <- SimilarityMeasure(M)

  expect_true(is.numeric(sim))
  expect_length(sim, nrow(M))
  expect_true(all(sim >= 0 & sim <= 1))
  # reference row compared to itself is a perfect match
  expect_equal(sim[1], 1)
})

test_that("SimilarityMeasure is exported and callable", {
  expect_true(is.function(SimilarityMeasure))
})
