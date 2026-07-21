test_that("projection() reproduces beta_hat exactly on Gaussian identity", {
  set.seed(6489)
  n <- 500
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4); x3 <- runif(n)
  y <- 1 + 0.5 * x1 - 0.3 * x2 + 0.2 * x3 + rnorm(n, sd = 0.5)
  d <- data.frame(y = y, x1 = x1, x2 = factor(x2), x3 = x3)
  beta.orig <- coef(lm(y ~ x1 + x2 + x3, d))

  syn <- projection(d, formula = y ~ x1 + x2 + x3, family = gaussian())
  beta.syn <- coef(lm(y ~ x1 + x2 + x3, syn$syn))

  expect_s3_class(syn, "synthexpln")
  expect_equal(unname(beta.syn), unname(beta.orig), tolerance = 1e-10)
})

test_that("projection() reproduces beta_hat on Gamma log within numerical tolerance", {
  set.seed(1)
  n <- 500
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  eta <- 1 + 0.4 * x1 - 0.2 * x2
  mu <- exp(eta)
  y <- rgamma(n, shape = 2, rate = 2 / mu)
  d <- data.frame(y = y, x1 = x1, x2 = factor(x2))
  beta.orig <- coef(glm(y ~ x1 + x2, family = Gamma("log"), d))

  syn <- projection(d, formula = y ~ x1 + x2, family = Gamma("log"))
  beta.syn <- coef(glm(y ~ x1 + x2, family = Gamma("log"), syn$syn))

  expect_equal(unname(beta.syn), unname(beta.orig), tolerance = 1e-6)
})

test_that("projection() reproduces beta_hat exactly on binomial logit with a valid [0, 1] response", {
  set.seed(6489)
  n <- 800
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(0.5 + 1.2 * x1 - 0.8 * x2))
  d <- data.frame(y = y, x1 = x1, x2 = factor(x2))
  beta.orig <- coef(glm(y ~ x1 + x2, family = binomial(), d))

  syn <- projection(d, formula = y ~ x1 + x2, family = binomial())

  # the release is a valid fractional binomial response in [0, 1], not the
  # unconstrained real-valued output of the old PATH-1 misroute
  expect_true(all(syn$syn$y >= 0 & syn$syn$y <= 1))
  # refit recovers beta_hat exactly (0.0% deviation, the continuous-family bar)
  beta.syn <- coef(suppressWarnings(glm(y ~ x1 + x2, family = binomial(), syn$syn)))
  expect_equal(unname(beta.syn), unname(beta.orig), tolerance = 1e-6)
})

test_that(".project.glm write-back stays aligned when the synthetic response has a zero the original lacks", {
  set.seed(6505)
  n <- 300
  x1 <- rnorm(n)
  mu <- exp(1 + 0.4 * x1)
  y <- rgamma(n, shape = 2, rate = 2 / mu)
  d <- data.frame(y = y, x1 = x1)
  beta.orig <- coef(glm(y ~ x1, family = Gamma("log"), d))

  # a zero in the synthetic response only: the fit subset is chosen by
  # zeros in the original data, the write-back must target the same rows
  d.syn <- d[sample(n, n, replace = TRUE), , drop = FALSE]
  d.syn <- `rownames<-`(d.syn, NULL)
  d.syn$y[42] <- 0

  expect_silent(
    v.proj <- synthexpln:::.project.glm(d.syn, d, "y", "x1", "gamma", "log"))
  expect_length(v.proj, n)
  expect_true(all(v.proj > 0))

  # near-exact, not exact: the injected zero is outside the Gamma
  # support, so the positivity clip after Pearson rescaling binds on
  # that row (documented residual-error source). The pre-fix failure
  # modes (recycling warning, surviving zero, shifted assignment) sit
  # orders of magnitude outside this tolerance.
  d.syn$y <- v.proj
  beta.syn <- coef(glm(y ~ x1, family = Gamma("log"), d.syn))
  expect_equal(unname(beta.syn), unname(beta.orig), tolerance = 1e-3)
})

test_that("projection() returns a synthexpln S3 object with expected slots", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  syn <- projection(d, y ~ x, gaussian())
  expect_named(syn, c("syn", "beta", "variant", "scope", "design.exact"), ignore.order = TRUE)
  expect_equal(syn$variant, "projection")
  expect_equal(syn$scope, "none")
})
