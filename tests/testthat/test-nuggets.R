test_that("create_data_nuggets compresses and preserves total weight", {
  set.seed(1)
  x <- rbind(matrix(rnorm(100 * 4, 0), ncol = 4),
             matrix(rnorm(100 * 4, 5), ncol = 4))
  dn <- create_data_nuggets(x, max_nuggets = 30, seed = 1)
  expect_s3_class(dn, "data_nugget")
  expect_lte(length(dn$weights), 30)
  expect_equal(sum(dn$weights), nrow(x))          # weights partition the data
  expect_equal(ncol(dn$centers), ncol(x))
  expect_true(all(dn$weights >= 1))
})

test_that("small samples fall back to one-nugget-per-observation", {
  x <- matrix(rnorm(20 * 3), ncol = 3)
  dn <- create_data_nuggets(x, seed = 1)        # n = 20 <= 50
  expect_equal(length(dn$weights), 20)
})

test_that("nugget feature weights rank signal above noise", {
  set.seed(2)
  signal <- c(rnorm(60, 0), rnorm(60, 6))       # separates two groups
  noise  <- rnorm(120)
  x <- cbind(signal = signal, noise = noise)
  dn <- create_data_nuggets(x, max_nuggets = 20, seed = 2)
  fw <- nugget_feature_weights(dn, "between")
  expect_gt(fw["signal"], fw["noise"])
  expect_true(all(fw >= 0))
})
