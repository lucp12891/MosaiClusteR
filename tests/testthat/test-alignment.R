# Regression tests for cross-source object alignment in the direct methods.

mk <- function(n = 20, p = 3, named = TRUE, seed = 1) {
  set.seed(seed)
  m <- matrix(stats::rnorm(n * p), n, p)
  if (named) rownames(m) <- sprintf("o%02d", seq_len(n))
  m
}

test_that(".align_sources aligns reordered rows and errors clearly otherwise", {
  a <- mk(); b <- mk(seed = 2)
  # reordered second source is realigned to the first
  br <- b[sample(nrow(b)), ]
  out <- MosaiClusteR:::.align_sources(list(a, br))
  expect_identical(rownames(out[[2]]), rownames(a))
  expect_equal(out[[2]], b[rownames(a), ])

  # no row names: positional, returned unchanged
  an <- mk(named = FALSE); bn <- mk(named = FALSE, seed = 2)
  expect_silent(MosaiClusteR:::.align_sources(list(an, bn)))

  # different object sets, different nrow, duplicated ids -> clear errors
  bmix <- mk(seed = 2); rownames(bmix)[1] <- "zzz"
  expect_error(MosaiClusteR:::.align_sources(list(a, bmix)), "same objects")
  expect_error(MosaiClusteR:::.align_sources(list(a, mk(n = 18))),
               "same number of objects")
  bdup <- mk(seed = 2); rownames(bdup)[2] <- rownames(bdup)[1]
  expect_error(MosaiClusteR:::.align_sources(list(bdup, a)), "duplicated")
})

test_that("ADC is invariant to row order of a source and handles unnamed input", {
  a <- mk(); b <- mk(seed = 2)
  lab_aligned <- mosaic_labels(ADC(list(a, b), distmeasure = "euclidean",
                                   normalize = FALSE), 3)
  br <- b[sample(nrow(b)), ]
  lab_reorder <- mosaic_labels(ADC(list(a, br), distmeasure = "euclidean",
                                   normalize = FALSE), 3)
  expect_equal(cluster_agreement(lab_aligned, lab_reorder)[["ARI"]], 1)

  # unnamed sources now work positionally instead of failing cryptically
  expect_silent(ADC(list(mk(named = FALSE), mk(named = FALSE, seed = 2)),
                    distmeasure = "euclidean", normalize = FALSE))
})
