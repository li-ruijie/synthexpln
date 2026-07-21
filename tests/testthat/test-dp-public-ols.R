test_that("gen.syn.dp.ols.public yields OLS estimate equal to beta_DP", {
  set.seed(1)
  n <- 500
  X <- cbind(1, rnorm(n), runif(n))
  beta <- c(0.5, 1.0, -0.3)
  y <- as.numeric(X %*% beta) + rnorm(n, sd = 0.5)
  d <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3])

  # B.y from the DGP above, not from y: |y| <= |0.5| + 1.0 * |x1| + 0.3 * |x2|
  # + |e| with x1 ~ N(0, 1) (|x1| <= 6), x2 ~ U(0, 1) (|x2| <= 1), and
  # e ~ N(0, 0.5) (|e| <= 3), giving 0.5 + 6 + 0.3 + 3 = 9.8. Round up to 10.
  syn <- gen.syn.dp.ols.public(d, y ~ x1 + x2,
                                epsilon = 10, delta = 1e-6, B.y = 10, seed = 1)
  beta.syn <- coef(lm(y ~ x1 + x2, syn$syn))
  expect_equal(unname(beta.syn), unname(syn$beta.dp), tolerance = 1e-8)
  expect_equal(syn$variant, "V2")
  expect_equal(syn$scope, "joint-release")

  # the refit also recovers the DP dispersion exactly, so the exact-fidelity
  # property covers both parameters
  expect_equal(summary(lm(y ~ x1 + x2, syn$syn))$sigma^2, syn$sigma2.dp,
               tolerance = 1e-8)
})

test_that("beta is calibrated to the GLOBAL sensitivity, not to DFBETA", {
  set.seed(11)
  n <- 300
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  B.y <- 10
  syn <- gen.syn.dp.ols.public(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6,
                                B.y = B.y, seed = 7)

  m.X <- model.matrix(y ~ x1 + x2, d)
  expect_equal(syn$Delta.beta, synthexpln:::.gs.beta.ols.public(m.X, B.y))
  expect_equal(syn$sigma.beta, syn$Delta.beta / syn$mu.beta)

  # the noise is isotropic on the raw coefficient scale: one sd for every
  # coordinate.  Column standardisation was a prop for the DFBETA noise and is
  # gone, so there is no $s.j left to divide by.
  expect_length(syn$sigma.beta, 1L)
  expect_null(syn$s.j)

  # regression guard.  The retired calibration was the DFBETA local sensitivity
  # on the standardised scale.  If anyone reinstates it, the noise scale changes
  # and this fails.
  v.s      <- ifelse(apply(m.X, 2, sd) > 1e-8, apply(m.X, 2, sd), 1)
  m.dfb    <- dfbeta(lm(y ~ x1 + x2, d))
  LS.tilde <- 2 * max(sqrt(rowSums(sweep(m.dfb, 2, v.s, "*")^2)))
  retired  <- LS.tilde * synthexpln:::.dp.gm.multiplier(8, 1e-6)
  expect_false(isTRUE(all.equal(syn$sigma.beta, retired)))
})

test_that("the two channels compose in one mu quadrature", {
  set.seed(1)
  d <- data.frame(y = rnorm(200), x = rnorm(200))
  eps <- 10; delta <- 1e-6
  # y ~ N(0, 1) by construction, so six sd bounds it: B.y = 6. A design
  # constant of the DGP on the line above, not a statistic of d$y.
  syn <- gen.syn.dp.ols.public(d, y ~ x,
                                epsilon = eps, delta = delta,
                                eps.split = c(beta = 0.8, sigma = 0.2),
                                B.y = 6)

  # eps.split holds fractions of the mu^2 budget, not of epsilon
  expect_equal(syn$mu.beta, syn$mu.total * sqrt(0.8))
  expect_equal(syn$mu.v,    syn$mu.total * sqrt(0.2))
  expect_equal(syn$mu.check, syn$mu.total)
  expect_equal(sqrt(syn$mu.beta^2 + syn$mu.v^2), syn$mu.total)

  # and mu.total is the exact dual of (eps, delta)
  expect_equal(synthexpln:::.gdp.delta(eps, syn$mu.total), delta,
               tolerance = 1e-9)

  # the absolute-epsilon split is gone: a Laplace sigma^2 channel could not
  # enter a mu quadrature, so it was Gaussianised
  expect_null(syn$epsilon.split)
})

test_that("the joint-release generators refuse a NULL B.y", {
  # Guard on the refusal itself. B.y sets the noise scale, so it must be a
  # public constant. The removed default read it off the sample
  # (1.5 * max(abs(y))), which makes the scale a function of the private data
  # and leaves no sensitivity bound to instantiate the theorem with, at any
  # epsilon. If someone reinstates that default, this test fails.
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  x.spec <- list(x = list(type = "continuous", bounds = c(-5, 5)))
  rx <- "B\\.y is required"

  expect_error(
    gen.syn.dp.ols.public(d, y ~ x, epsilon = 10, delta = 1e-6, B.y = NULL),
    regexp = rx)
  expect_error(
    gen.syn.dp.ridge.public(d, y ~ x, lambda = 1, epsilon = 10, delta = 1e-6,
                            B.y = NULL),
    regexp = rx)
  # V4's sound path needs no B.y at all (subsample-and-aggregate takes its
  # sensitivity from the public coefficient box, so it holds however unbounded
  # y is).  Only its retired legacy path still demands one.
  expect_error(
    suppressWarnings(
      gen.syn.dp.full(d, y ~ x, x.spec = x.spec, epsilon = 10, delta = 1e-6,
                      sens.method = "local", B.y = NULL)),
    regexp = rx)

  # B.y is not defaulted anywhere in the call chain either: dp.project()
  # forwards it, so omitting it still reaches the refusal the user sees.
  expect_error(
    dp.project(d, y ~ x, epsilon = 10, delta = 1e-6, x.public = TRUE),
    regexp = rx)
})

test_that("sigma2.dp draws at the K_X-corrected GAUSSIAN scale", {
  set.seed(21)
  n <- 60
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$x1[1] <- 8                       # leverage row: K_X > 1
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  B.y <- 10
  eps <- 50; delta <- 1e-6
  syn <- gen.syn.dp.ols.public(d, y ~ x1 + x2, epsilon = eps, delta = delta,
                                eps.split = c(beta = 0.2, sigma = 0.8),
                                B.y = B.y, seed = 7)

  X <- model.matrix(y ~ x1 + x2, d)
  sigma2.hat <- summary(lm(y ~ x1 + x2, d))$sigma^2
  Delta   <- synthexpln:::.gs.sigma2.public(X, B.y)
  mu.v    <- synthexpln:::.gdp.mu.from.eps.delta(eps, delta) * sqrt(0.8)

  # Gaussian, not Laplace: a Laplace channel is pure eps-DP and does not enter
  # the mu quadrature the beta channel is calibrated in.
  set.seed(7 + 2L)
  expect_equal(syn$sigma2.dp,
               max(sigma2.hat + rnorm(1, 0, Delta / mu.v), 1e-6))

  # regression guard, on the scale and not on the realised draw.  A draw
  # comparison can be silently non-discriminating: when the noise is large
  # relative to sigma2.hat, both the corrected and the naive draw reach the
  # 1e-6 truncation floor and compare equal.  The claim is about the
  # calibration, so test the calibration.
  naive.scale <- (4 * B.y^2 / (n - ncol(X))) / mu.v
  expect_equal(syn$sigma.v, Delta / mu.v)
  expect_false(isTRUE(all.equal(syn$sigma.v, naive.scale)))
  # K_X > 1 on this leveraged design, so dropping it under-noises the channel
  expect_gt(Delta / (4 * B.y^2 / (n - ncol(X))), 1)
})
