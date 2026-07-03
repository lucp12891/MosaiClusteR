# Tests for the navigation / comparison utilities ported into utils-navigation.r

# Build a realistic input: two single-source clustering results on mosaic_toy,
# each a list with $Clust (agnes) and $DistM and attr "method" == "Cluster".
make_results <- function(nrclusters = 5) {
  data(mosaic_toy, package = "MosaiClusteR")
  L <- mosaic_toy$List
  res <- list(
    Cluster(L[[1]], type = "data", distmeasure = "euclidean"),
    Cluster(L[[2]], type = "data", distmeasure = "euclidean")
  )
  res
}

test_that("the navigation functions are exported and available", {
  expect_true(is.function(ReorderToReference))
  expect_true(is.function(FindCluster))
  expect_true(is.function(FindElement))
  expect_true(is.function(SharedComps))
})

test_that("ReorderToReference relabels clusterings to the reference", {
  res <- make_results()
  M <- suppressMessages(ReorderToReference(
    res, nrclusters = 5, fusionsLog = FALSE,
    weightclust = FALSE, names = c("S1", "S2")
  ))

  expect_true(is.matrix(M))
  # one row per method, one column per object
  expect_equal(nrow(M), 2L)
  expect_equal(ncol(M), nrow(res[[1]]$DistM))
  expect_equal(rownames(M), c("S1", "S2"))
  # columns are labelled by the reference ordering
  expect_setequal(colnames(M), rownames(res[[1]]$DistM))
  # the reference row is itself relabelled to 1..k in dendrogram order
  expect_true(all(M[1, ] %in% seq_len(5)))
})

test_that("FindCluster returns the members of a chosen cluster", {
  res <- make_results()
  comps <- suppressMessages(FindCluster(
    List = res, nrclusters = 5, select = c(1, 2),
    fusionsLog = FALSE, weightclust = FALSE
  ))

  expect_true(is.character(comps))
  # all returned objects are valid object names
  expect_true(all(comps %in% rownames(res[[1]]$DistM)))

  # cross-check against ReorderToReference: method 1, cluster 2
  M <- suppressMessages(ReorderToReference(
    res, nrclusters = 5, fusionsLog = FALSE, weightclust = FALSE, names = NULL
  ))
  expect_setequal(comps, names(which(M[1, ] == 2)))
})

test_that("FindElement recursively extracts a named element", {
  obj <- list(a = 1, b = list(TopDE = 42, c = list(TopDE = 7)))
  found <- FindElement(what = "TopDE", object = obj)

  expect_true(is.list(found))
  expect_true(length(found) >= 1)
  expect_true(all(grepl("^TopDE_", names(found))))
  expect_true(42 %in% unlist(found))

  # absent name yields an empty list
  expect_length(FindElement(what = "DOES_NOT_EXIST", object = obj), 0L)
})

test_that("SharedComps returns shared objects per cluster", {
  res <- make_results()
  sc <- suppressMessages(SharedComps(
    List = res, nrclusters = 5, fusionsLog = FALSE,
    weightclust = FALSE, names = c("S1", "S2")
  ))

  expect_true(is.list(sc))
  expect_length(sc, 5L)
  expect_true(all(grepl("^Cluster ", names(sc))))
  # every shared object is a genuine object name
  allnames <- rownames(res[[1]]$DistM)
  expect_true(all(unlist(sc) %in% allnames))
})

test_that("functions error informatively on bad input", {
  # an empty list cannot be reordered
  expect_error(suppressMessages(ReorderToReference(list(), nrclusters = 3)))
})
