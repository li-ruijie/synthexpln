# V1 is the retired analysis-scope route.  It calibrates to the DFBETA local
# sensitivity and carries no formal DP, so it warns on every call.  The tests
# below exercise its MECHANICS, which are still what they were; the warning is
# asserted once, here, and suppressed at the other call sites.

test_that("gen.syn.dp.projected WARNS that it carries no formal DP", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  expect_warning(gen.syn.dp.projected(d, y ~ x, epsilon = 10, delta = 1e-6),
                 "NO formal \\(eps, delta\\)-DP guarantee")
  # and it points at the routes that do carry one
  w <- tryCatch({ gen.syn.dp.projected(d, y ~ x, epsilon = 10); NULL },
                warning = conditionMessage)
  expect_match(w, "gen\\.syn\\.dp\\.projected\\.subagg")
  expect_match(w, "gen\\.syn\\.dp\\.ols\\.public")
})

test_that("gen.syn.dp.projected recovers beta_DP exactly from synthetic data", {
  set.seed(6489)
  n <- 500
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  y <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n, sd = 0.5)
  d <- data.frame(y = y, x1 = x1, x2 = factor(x2))

  syn <- suppressWarnings(
    gen.syn.dp.projected(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6, seed = 1))
  beta.syn <- coef(lm(y ~ x1 + x2, syn$syn))
  expect_equal(unname(beta.syn), unname(syn$beta.dp), tolerance = 1e-8)
  expect_equal(syn$variant, "V1")
  expect_equal(syn$scope, "analysis-level")
})

test_that("gen.syn.dp.projected respects epsilon via Gaussian noise scale", {
  set.seed(1)
  n <- 500
  d <- data.frame(y = rnorm(n), x = rnorm(n))
  sigs <- suppressWarnings(sapply(c(1, 10), function(eps) {
    replicate(50, {
      s <- gen.syn.dp.projected(d, y ~ x, epsilon = eps, delta = 1e-6)
      s$beta.dp[2]
    })
  }))
  # variance ratio ~ (m(1)/m(10))^2 = 61 under the calibrated multiplier
  # (was (10/1)^2 = 100 under the classical 1/eps scaling); same 50% slack
  expect_gt(var(sigs[, 1]) / var(sigs[, 2]), 30)
})
