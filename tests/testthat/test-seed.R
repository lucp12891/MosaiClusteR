# Reproducibility: seeded stochastic methods give identical results.

test_that("M_ABC is reproducible with seed=", {
  data(mosaic_toy); L <- mosaic_toy$List
  d1 <- M_ABC(L, numsim = 40, NC = 3, seed = 42)$DistM
  d2 <- M_ABC(L, numsim = 40, NC = 3, seed = 42)$DistM
  expect_equal(d1, d2)
  d3 <- M_ABC(L, numsim = 40, NC = 3, seed = 7)$DistM
  expect_false(isTRUE(all.equal(d1, d3)))
})

test_that("ADEC is reproducible with seed= when sampling features", {
  data(mosaic_toy); L <- mosaic_toy$List
  a1 <- ADEC(L, distmeasure = "euclidean", normalize = FALSE,
             t = 4, r = 50, nrclusters = 3, seed = 1)$DistM
  a2 <- ADEC(L, distmeasure = "euclidean", normalize = FALSE,
             t = 4, r = 50, nrclusters = 3, seed = 1)$DistM
  expect_equal(a1, a2)
})
