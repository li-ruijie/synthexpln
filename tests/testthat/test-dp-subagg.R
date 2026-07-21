test_that("gen.syn.dp.projected.subagg gives worst-case DP and recovers beta_DP (gaussian)", {
  set.seed(4321)
  n  <- 2000
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  y  <- 1 + 0.8 * x1 - 0.5 * x2 + 0.3 * x3 + rnorm(n)
  d  <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

  # sigma^2 = 1 by the DGP on the line above.  A design constant, declared.
  disp <- dp.dispersion.public(1, source = "unit residual variance, by the DGP")
  syn <- gen.syn.dp.projected.subagg(d, y ~ x1 + x2 + x3, epsilon = 5, m = 20,
           family = gaussian(), clip.lo = -5, clip.hi = 5,
           dispersion = disp, seed = 6489)

  expect_s3_class(syn, "synthexpln")
  expect_equal(syn$variant, "D2")

  # block aggregate close to the direct OLS estimate
  b.ols <- coef(lm(y ~ x1 + x2 + x3, d))
  expect_lt(max(abs(syn$beta.agg - b.ols)), 0.1)

  # exact inferential fidelity survives the fix: the refit still recovers
  # beta_DP, and now the declared dispersion too.
  expect_equal(unname(coef(lm(y ~ x1 + x2 + x3, syn$syn))), unname(syn$beta.dp),
               tolerance = 1e-8)
  expect_equal(summary(lm(y ~ x1 + x2 + x3, syn$syn))$sigma^2, 1,
               tolerance = 1e-6)

  # worst-case Gaussian sigma: s_j * sqrt(q) / mu, calibrated by the exact dual
  # (it used to use .dp.gm.multiplier, which falls back to the classical bound
  # below eps = 1 and so over-noised there)
  p      <- length(b.ols)
  sj     <- (5 - (-5)) / 20
  sigmaQ <- synthexpln:::.gdp.sigma(sqrt(p), 5, 1e-6)
  expect_equal(unname(syn$sigma[1]), sj * sigmaQ, tolerance = 1e-10)
})

test_that("subsample-aggregate noise scales with the exact dual", {
  set.seed(11)
  n  <- 1500
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- 1 + 0.5 * x1 - 0.3 * x2 + rnorm(n)
  d  <- data.frame(y = y, x1 = x1, x2 = x2)
  disp <- dp.dispersion.public(1, source = "unit residual variance, by the DGP")
  gen <- function(eps) gen.syn.dp.projected.subagg(
    d, y ~ x1 + x2, epsilon = eps, m = 20, family = gaussian(),
    clip.lo = -5, clip.hi = 5, dispersion = disp, seed = 2)
  s.lo <- gen(1); s.hi <- gen(10)

  # sigma is proportional to 1 / mu*(eps, delta), the exact dual
  r.expect <- synthexpln:::.gdp.mu.from.eps.delta(10, 1e-6) /
    synthexpln:::.gdp.mu.from.eps.delta(1, 1e-6)
  expect_equal(unname(s.lo$sigma[1] / s.hi$sigma[1]), r.expect, tolerance = 1e-8)
  expect_gt(r.expect, 1)   # more budget, less noise
})

test_that("the dispersion must be DECLARED, and its mode is honoured", {
  set.seed(51)
  n <- 800
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n, 0, 2)

  # A GLM has two parameters and DP must cover both.  Silence is not an option.
  expect_error(
    gen.syn.dp.projected.subagg(d, y ~ x1 + x2, epsilon = 10,
                                family = gaussian(), clip.lo = -5, clip.hi = 5),
    "dispersion is required")

  # public: free.  The value is reused as-is and beta_DP is p-dimensional.
  pub <- gen.syn.dp.projected.subagg(
    d, y ~ x1 + x2, epsilon = 10, family = gaussian(),
    clip.lo = -5, clip.hi = 5, seed = 3,
    dispersion = dp.dispersion.public(4, source = "known instrument variance"))
  expect_equal(pub$sigma2.dp, 4)
  expect_equal(summary(lm(y ~ x1 + x2, pub$syn))$sigma^2, 4, tolerance = 1e-6)

  # private: paid for.  sigma^2 rides along as a (p+1)-th coordinate, so the L2
  # sensitivity rises from sqrt(p) to sqrt(p+1) and beta_DP gets NOISIER.
  prv <- gen.syn.dp.projected.subagg(
    d, y ~ x1 + x2, epsilon = 10, family = gaussian(),
    clip.lo = -5, clip.hi = 5, seed = 3,
    dispersion = dp.dispersion.private(0, 16, source = "instrument caps sd at 4"))
  expect_true(prv$sigma2.dp >= 0 && prv$sigma2.dp <= 16 + 1e-6)
  expect_gt(prv$sigma[[1]], pub$sigma[[1]])
  expect_equal(unname(prv$sigma[[1]] / pub$sigma[[1]]), sqrt(4 / 3),
               tolerance = 1e-8)   # sqrt((p+1)/p), p = 3

  # a family with no free dispersion must not be handed one
  db <- data.frame(x1 = rnorm(200))
  db$y <- rbinom(200, 1, plogis(db$x1))
  expect_error(
    gen.syn.dp.projected.subagg(db, y ~ x1, epsilon = 10, family = binomial(),
                                clip.lo = -3, clip.hi = 3,
                                dispersion = dp.dispersion.public(1, source = "x")),
    "no free dispersion")
})

test_that("the release does not echo the empirical residual distribution", {
  # The old construction was y_virt = y + mu(beta_DP) - mu(beta_hat), which
  # retains the real y, so the release reproduced the real residuals in scale
  # and in shape.  A release that carries the empirical residual law does not
  # depend on the data "only through beta_DP", which is the hypothesis
  # the post-processing clause needs.
  #
  # Feed the generator a deliberately heavy-tailed, skewed residual law.  If the
  # old construction is ever reinstated, that shape flows into the release and
  # this test fails.
  set.seed(21)
  n <- 1000
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + (rt(n, df = 3) * 2 + rexp(n, 1))

  kurt <- function(z) mean((z - mean(z))^4) / stats::sd(z)^4
  skew <- function(z) mean((z - mean(z))^3) / stats::sd(z)^3

  r.real <- residuals(lm(y ~ x1 + x2, d))
  expect_gt(kurt(r.real), 5)          # the input is heavy-tailed

  syn <- gen.syn.dp.projected.subagg(
    d, y ~ x1 + x2, epsilon = 10, family = gaussian(),
    clip.lo = -5, clip.hi = 5, seed = 2,
    dispersion = dp.dispersion.public(16, source = "declared, not measured"))
  r.syn <- residuals(lm(y ~ x1 + x2, syn$syn))

  # the release is Gaussian at the declared dispersion and carries no trace of
  # the t3 + exponential input: kurtosis ~ 3, skew ~ 0
  expect_lt(kurt(r.syn), 3.6)
  expect_lt(abs(skew(r.syn)), 0.35)
  # and its scale is the declared one, not the sample's
  expect_equal(stats::sd(r.syn), 4, tolerance = 0.05)
  expect_false(isTRUE(all.equal(stats::sd(r.syn), stats::sd(r.real),
                                tolerance = 0.05)))
})

test_that(".smooth.sens.beta.loglink (D1) is finite, well-conditioned, and smooth >= local", {
  set.seed(7)
  n  <- 1500
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- rgamma(n, 2, scale = exp(0.6 + 0.4 * x1 - 0.3 * x2) / 2)
  mo <- glm(y ~ x1 + x2, data.frame(y = y, x1 = x1, x2 = x2),
            family = Gamma(link = "log"))
  s <- .smooth.sens.beta.loglink(mo, eps = 5, delta = 1e-6)
  expect_true(is.finite(s$S.smooth))
  expect_gt(s$lambda.min, 0)
  expect_gte(s$S.smooth, s$S.local)   # smooth sensitivity inflates the local one
})

test_that(".smooth.sens.beta.loglink.std (D1 viable) well-conditions a large-scale design", {
  set.seed(13)
  n  <- 1500
  x1 <- rnorm(n)
  x2 <- round(exp(rnorm(n, 7.5, 0.5)))     # income-like: large native scale
  d  <- data.frame(y = rgamma(n, 2, scale = exp(0.6 + 0.4 * x1 - 3e-4 * x2) / 2),
                   x1 = x1, x2 = x2)
  mo <- glm(y ~ x1 + x2, d, family = Gamma(link = "log"))

  raw <- .smooth.sens.beta.loglink(mo, eps = 5, delta = 1e-6)
  std <- .smooth.sens.beta.loglink.std(mo, eps = 5, delta = 1e-6,
           x.scale = c(`(Intercept)` = 1, x1 = 1, x2 = 2000))

  # the raw envelope diverges on the large-scale x2, whereas the public
  # standardisation well-conditions the Fisher information and reduces the
  # smooth to local ratio.
  expect_gt(raw$ratio, 1e3)
  expect_lt(std$ratio, raw$ratio)
  expect_true(all(is.finite(std$sigma.native)) && all(std$sigma.native > 0))
  # noise on the large-scale coefficient is shrunk by its public scale.
  expect_lt(std$sigma.native[["x2"]], std$sigma.native[["x1"]])
})
