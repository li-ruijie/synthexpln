test_that("projection.twopart() preserves the zero pattern and reproduces coefficients and SEs", {
  set.seed(6489)
  n  <- 500
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.4)
  mu <- exp(0.7 + 0.5 * x1 - 0.3 * x2)
  y  <- ifelse(rbinom(n, 1, 0.6) == 1, rgamma(n, shape = 2, rate = 2 / mu), 0)
  d  <- data.frame(y = y, x1 = x1, x2 = factor(x2))

  fit.pos   <- glm(y ~ x1 + x2, family = Gamma("log"), data = d[d$y > 0, ])
  beta.orig <- coef(fit.pos)
  se.orig   <- summary(fit.pos)$coefficients[, "Std. Error"]

  syn <- projection.twopart(d, y ~ x1 + x2, family = Gamma("log"), seed = 1)
  expect_s3_class(syn, "synthexpln")

  # the zero/positive indicator is released unchanged, row for row
  expect_true(all((syn$syn$y > 0) == (d$y > 0)))

  fit.syn <- glm(y ~ x1 + x2, family = Gamma("log"), data = syn$syn[syn$syn$y > 0, ])
  # coefficients reproduced: the projection matches the score at beta_hat.
  # Observed max abs deviation ~2e-6 (the Newton score tolerance of 1e-8 plus
  # the refit's own convergence); 1e-5 leaves headroom without slackening it.
  expect_equal(unname(coef(fit.syn)), unname(beta.orig), tolerance = 1e-5)

  # standard errors reproduced too: the positive design is released unchanged
  # (X* = X), so phi (X'WX)^{-1} is the same matrix.  Observed max relative
  # deviation ~1e-6, far smaller than the coefficient bar.
  se.syn <- summary(fit.syn)$coefficients[, "Std. Error"]
  expect_equal(unname(se.syn), unname(se.orig), tolerance = 1e-5)
})

test_that("projection.twopart() rejects a non-positive-support link", {
  d <- data.frame(y = c(0, 1, 2, 0, 3, 1.5), x = rnorm(6))
  expect_error(projection.twopart(d, y ~ x, family = gaussian()),
               "positive-support link")
})

test_that("projection.twopart() rejects a response with no zero block", {
  set.seed(1)
  d <- data.frame(y = rgamma(50, shape = 2, rate = 1), x = rnorm(50))   # all positive
  expect_error(projection.twopart(d, y ~ x, family = Gamma("log")),
               "both zeros and positive")
})

test_that("projection.twopart() returns a synthexpln object with the twopart tags", {
  set.seed(1)
  n  <- 200
  x  <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x)
  y  <- ifelse(rbinom(n, 1, 0.6) == 1, rgamma(n, shape = 2, rate = 2 / mu), 0)
  d  <- data.frame(y = y, x = x)

  syn <- projection.twopart(d, y ~ x, family = Gamma("log"))
  expect_named(syn, c("syn", "beta", "variant", "scope", "design.exact"), ignore.order = TRUE)
  expect_equal(syn$variant, "twopart")
  expect_equal(syn$scope, "none")
})
