# Regression tests for cluster_agreement() edge cases.

test_that("identical partitions agree perfectly, including degenerate cases", {
  # normal identical-up-to-relabelling
  expect_equal(unname(cluster_agreement(c(1, 1, 2, 2), c(2, 2, 1, 1))),
               c(1, 1, 1))
  # single-cluster, identical (ARI/NMI formulae are otherwise 0/0)
  expect_equal(unname(cluster_agreement(rep(1, 6), rep(1, 6))), c(1, 1, 1))
  # all singletons, identical
  expect_equal(unname(cluster_agreement(1:6, 1:6)), c(1, 1, 1))
})

test_that("genuinely different partitions do not report perfect agreement", {
  ag <- cluster_agreement(c(1, 1, 2, 2, 3, 3), c(1, 2, 3, 1, 2, 3))
  expect_lt(ag[["ARI"]], 1)
  # a single cluster vs a real partition is not perfect agreement
  expect_lt(cluster_agreement(rep(1, 6), c(1, 1, 1, 2, 2, 2))[["ARI"]], 1)
})

test_that("missing labels and length mismatch raise clear errors", {
  expect_error(cluster_agreement(c(1, 1, NA, 2), c(1, 1, 2, 2)),
               "must not contain NA")
  expect_error(cluster_agreement(c(1, 1, 2), c(1, 2)),
               "equal length")
})

test_that("agreement is invariant to label permutation", {
  set.seed(1)
  a <- sample(1:3, 30, replace = TRUE)
  perm <- c(`1` = 2L, `2` = 3L, `3` = 1L)
  b <- unname(perm[as.character(a)])
  expect_equal(unname(cluster_agreement(a, b)), c(1, 1, 1))
})
