# The public block-count rule: m = floor(n_eff / (5 p)), the largest
# block count leaving five rows per released coefficient per block, with the
# two-part positive channel scaling n by a published zero rate carried as a
# dp.rate provenance object.  A realised rate must not reach the rule.

test_that("dp.subagg.blocks applies the five-rows-per-coefficient cap", {
  expect_identical(dp.subagg.blocks(2000, 6), 66L)
  expect_identical(dp.subagg.blocks(2000, 6, dp.rate.public(0.5, "test")), 33L)
  expect_error(dp.subagg.blocks(2000, 6, rate = 0.5), "dp.rate.public")
  expect_error(dp.subagg.blocks(50, 6))
  expect_error(dp.rate.public(1.0, "test"))
  expect_error(dp.rate.public(0.5, ""))
})

test_that("the engine self-defaults by the rule and refuses it under block", {
  set.seed(1)
  d <- data.frame(x1 = rnorm(300), x2 = rnorm(300))
  d$y <- 1 + d$x1 - d$x2 + rnorm(300)
  l <- synthexpln:::.subsample.aggregate.beta.loglink(
    d, y ~ x1 + x2, stats::gaussian(), "gaussian", "identity",
    eps = 5, delta = 1e-6, m = NULL, clip.lo = -5, clip.hi = 5, seed = 1)
  expect_equal(l$m, 20L)                       # floor(300 / (5 * 3))
  blk <- sample(rep_len(1:10, 300))
  expect_error(
    synthexpln:::.subsample.aggregate.beta.loglink(
      d, y ~ x1 + x2, stats::gaussian(), "gaussian", "identity",
      eps = 5, delta = 1e-6, m = NULL, clip.lo = -5, clip.hi = 5, seed = 1,
      block = blk),
    "m is required")
})

test_that("the ZI generator computes the per-channel pair from a published rate", {
  set.seed(2)
  n <- 600L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  pos <- rbinom(n, 1L, 0.5) == 1L
  d$y <- ifelse(pos, rgamma(n, 2, 1), 0)
  box <- c(`(Intercept)` = 3, x1 = 3, x2 = 3)
  expect_error(
    gen.syn.dp.projected.zi(d, y ~ x1 + x2, epsilon = 5,
      clip.lo = -box, clip.hi = box, clip.z.lo = -box, clip.z.hi = box,
      dispersion = dp.dispersion.public(1, "test"), seed = 1),
    "pos.rate")
  res <- gen.syn.dp.projected.zi(d, y ~ x1 + x2, epsilon = 5,
    pos.rate = dp.rate.public(0.5, "test rate"),
    clip.lo = -box, clip.hi = box, clip.z.lo = -box, clip.z.hi = box,
    dispersion = dp.dispersion.public(1, "test"), seed = 1)
  expect_identical(res$m, c(z = 40L, p = 20L))
  # floor(600 / 15) = 40 (zero-model), floor(600 * 0.5 / 15) = 20 (positive).
})
