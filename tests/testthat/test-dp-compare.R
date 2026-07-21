# dp.compare(): the one-call comparison report.  It binds the existing
# fidelity metrics with a dependence section and the privacy accounting the
# release carries, and its docs are explicit that a dataset comparison
# quantifies utility while DP is a property of the mechanism.

mk.cmp.d <- function(n = 300L, seed = 9L) {
  set.seed(seed)
  g  <- factor(sample(c("a", "b"), n, replace = TRUE), levels = c("a", "b"))
  x1 <- rnorm(n, c(a = 30, b = 60)[as.character(g)], 10)
  x2 <- 0.5 * x1 + rnorm(n, 20, 8)
  y  <- 1 + 0.05 * x1 - 0.02 * x2 + 2 * (g == "b") + rnorm(n)
  data.frame(y = y, x1 = x1, x2 = x2, g = g)
}
cmp.spec <- list(x1 = list(type = "continuous", bounds = c(-10, 110)),
                 x2 = list(type = "continuous", bounds = c(-10, 110)),
                 g  = list(type = "categorical", levels = c("a", "b")))

test_that("dp.compare reports all four sections for a subagg release", {
  d   <- mk.cmp.d()
  syn <- gen.syn.dp.full(d, y ~ x1 + x2 + g, x.spec = cmp.spec,
                          epsilon = 10, delta = 1e-6,
                          clip.lo = -10, clip.hi = 10,
                          disp.lo = 0, disp.hi = 9, seed = 2)
  cmp <- dp.compare(syn, d, y ~ x1 + x2 + g)

  expect_s3_class(cmp, "dp.comparison")
  # inferential: delegates to fidelity.inferential
  expect_true(is.finite(cmp$inferential$max.abs.pct))
  # marginals: one entry per shared column (y, x1, x2, g)
  expect_setequal(names(cmp$marginals), c("y", "x1", "x2", "g"))
  # dependence: continuous block error and the eta table (g x {y, x1, x2})
  expect_true(is.finite(cmp$dependence$cor.frob))
  expect_equal(nrow(cmp$dependence$eta), 3L)
  expect_true(all(c("eta.orig", "eta.syn") %in% names(cmp$dependence$eta)))
  # privacy: the mu accounting echoed, with the MI ceiling Phi(mu / sqrt(2))
  expect_equal(cmp$privacy$mu.total, syn$mu.total)
  expect_equal(cmp$privacy$mi.ceiling, pnorm(syn$mu.total / sqrt(2)))
  expect_match(cmp$privacy$guarantee, "mu-GDP")
  # print method runs and returns invisibly
  expect_output(ret <- withVisible(print(cmp)), "Privacy accounting")
  expect_false(ret$visible)
})

test_that("dp.compare reports NO guarantee for the retired local path", {
  d <- mk.cmp.d()
  suppressWarnings(
    syn <- gen.syn.dp.full(d, y ~ x1 + x2 + g, x.spec = cmp.spec,
                            epsilon = 10, sens.method = "local",
                            B.y = 20, seed = 3))
  cmp <- dp.compare(syn, d, y ~ x1 + x2 + g, calibrate = FALSE)
  expect_match(cmp$privacy$guarantee, "NONE")
  expect_true(is.na(cmp$privacy$mi.ceiling))
})

test_that("dp.compare handles a non-DP projection release", {
  d   <- mk.cmp.d()
  syn <- projection(d, y ~ x1 + x2 + g)
  cmp <- suppressMessages(dp.compare(syn, d, y ~ x1 + x2 + g))
  expect_match(cmp$privacy$guarantee, "Not a DP release")
  expect_true(is.na(cmp$privacy$mi.ceiling))
  # the projection reproduces the coefficients, so deviation is ~0
  expect_lt(cmp$inferential$max.abs.pct, 1e-4)
})

test_that("dp.compare validates its inputs", {
  d <- mk.cmp.d()
  expect_error(dp.compare(list(), d, y ~ x1), "synthexpln release")
  syn <- projection(d, y ~ x1 + x2 + g)
  expect_error(dp.compare(syn, d), "formula")
})
