test_that("fidelity.inferential reports zero deviation for exact projection", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  syn <- projection(d, y ~ x)
  m <- fidelity.inferential(syn, orig = d, formula = y ~ x)
  expect_equal(m$max.abs.pct, 0, tolerance = 1e-6)
  expect_equal(m$sign.agreement, 1)
  expect_equal(m$sig.agreement, 1)
})

test_that("fidelity.inferential computes calibrated SE", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  # V1 warns that it carries no formal DP (asserted in test-dp-analysis.R).
  # This test is about fidelity.inferential(), so the warning is noise here.
  syn <- suppressWarnings(
    gen.syn.dp.projected(d, y ~ x, epsilon = 10, delta = 1e-6))
  m <- fidelity.inferential(syn, orig = d, formula = y ~ x, calibrate = TRUE)
  expect_true(all(m$se.corr >= m$se.syn))
})

test_that("fidelity.distributional returns KS statistic per column", {
  set.seed(1)
  orig <- data.frame(a = rnorm(100), b = rnorm(100))
  syn  <- data.frame(a = rnorm(100), b = rnorm(100))
  m <- fidelity.distributional(syn, orig, type = "ks")
  expect_named(m, c("a", "b"))
  expect_true(all(m >= 0 & m <= 1))
})
