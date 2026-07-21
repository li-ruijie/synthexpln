test_that("gen.syn.dp.ridge.public recovers ridge estimator", {
  set.seed(1)
  n <- 500
  X <- matrix(rnorm(n * 3), n, 3)
  beta <- c(1, -0.5, 0.3)
  y <- as.numeric(X %*% beta) + rnorm(n, sd = 0.5)
  d <- data.frame(y = y, x1 = X[, 1], x2 = X[, 2], x3 = X[, 3])

  # B.y from the DGP above, not from y. The x are iid N(0, 1) and e ~ N(0, 0.5),
  # all independent, so y is centred normal with
  # sd = sqrt(1^2 + 0.5^2 + 0.3^2 + 0.5^2) = sqrt(1.59) = 1.26, and six sd give
  # |y| <= 7.57. Round up to 10.
  syn <- gen.syn.dp.ridge.public(d, y ~ x1 + x2 + x3,
                                   lambda = 1,
                                   epsilon = 10, delta = 1e-6, B.y = 10,
                                   seed = 1)
  X.s <- model.matrix(~ x1 + x2 + x3, syn$syn)
  y.s <- syn$syn$y
  beta.ridge <- solve(crossprod(X.s) + diag(1, ncol(X.s)), crossprod(X.s, y.s))
  expect_equal(as.numeric(beta.ridge), unname(syn$beta.dp), tolerance = 1e-8)
  expect_equal(syn$variant, "V3")
})

test_that("ridge beta is calibrated to the GLOBAL sensitivity, not the ridge DFBETA", {
  set.seed(31)
  n <- 400
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  B.y <- 10; lambda <- 1
  syn <- gen.syn.dp.ridge.public(d, y ~ x1 + x2, lambda = lambda,
                                  epsilon = 10, delta = 1e-6, B.y = B.y, seed = 7)

  m.X <- model.matrix(y ~ x1 + x2, d)
  expect_equal(syn$Delta.beta,
               synthexpln:::.gs.beta.ridge.public(m.X, lambda, B.y))
  expect_equal(syn$sigma.beta, syn$Delta.beta / syn$mu.beta)

  # the standardisation and the 3-sigma truncation went with the local
  # sensitivity they were propping up
  expect_length(syn$sigma.beta, 1L)
  expect_null(syn$s.j)
  expect_null(syn$LS)

  # regression guard.  The retired calibration was the ridge DFBETA, read off the
  # realised residuals.  If anyone reinstates it, this fails.
  v.s     <- ifelse(apply(m.X, 2, sd) > 1e-8, apply(m.X, 2, sd), 1)
  m.X.std <- sweep(m.X, 2, v.s, "/")
  LS      <- synthexpln:::.sens.beta.ridge(m.X.std, d$y, lambda)$LS
  retired <- LS * synthexpln:::.dp.gm.multiplier(8, 1e-6)
  expect_false(isTRUE(all.equal(syn$sigma.beta, retired)))
})

test_that("the ridge global sensitivity reduces to the OLS one as lambda -> 0", {
  set.seed(32)
  n <- 200
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- rnorm(n)
  m.X <- model.matrix(y ~ x1 + x2, d)
  B.y <- 10

  # Delta_beta = 2 B_y max_i ||A^-1 x_i||, A = X'X + lambda I.  At lambda = 0,
  # A = X'X and this is the OLS global sensitivity.
  expect_equal(synthexpln:::.gs.beta.ridge.public(m.X, 1e-10, B.y),
               synthexpln:::.gs.beta.ols.public(m.X, B.y),
               tolerance = 1e-6)
  # and the penalty shrinks the sensitivity: more regularisation, less leverage
  # for any one row
  expect_lt(synthexpln:::.gs.beta.ridge.public(m.X, 100, B.y),
            synthexpln:::.gs.beta.ridge.public(m.X, 1, B.y))
})

test_that("ridge sigma2.dp draws at the M-matrix-corrected GAUSSIAN scale", {
  set.seed(22)
  n <- 60
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$x1[1] <- 8                       # leverage row: K > 1
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  B.y <- 10
  lambda <- 1
  eps <- 50; delta <- 1e-6
  syn <- gen.syn.dp.ridge.public(d, y ~ x1 + x2, lambda = lambda,
                                  epsilon = eps, delta = delta,
                                  eps.split = c(beta = 0.2, sigma = 0.8),
                                  B.y = B.y, seed = 7)

  X <- model.matrix(y ~ x1 + x2, d)
  A.inv <- solve(crossprod(X) + lambda * diag(ncol(X)))
  b.r <- as.numeric(A.inv %*% crossprod(X, d$y))
  sigma2.hat <- sum((d$y - X %*% b.r)^2) / (n - ncol(X))
  Delta <- synthexpln:::.gs.sigma2.ridge(X, lambda, B.y)
  mu.v  <- synthexpln:::.gdp.mu.from.eps.delta(eps, delta) * sqrt(0.8)

  # Gaussian, not Laplace: the release is stated in mu, and a Laplace channel
  # does not enter a mu quadrature at all.
  set.seed(7 + 2L)
  expect_equal(syn$sigma2.dp,
               max(sigma2.hat + rnorm(1, 0, Delta / mu.v), 1e-6))
  expect_equal(syn$mu.check, syn$mu.total)

  # regression guard, on the scale and not on the realised draw.  Comparing
  # draws is not discriminating here: sigma2.hat is about 1 and the corrected
  # noise is several times that, so both the corrected and the naive draw land
  # on the 1e-6 truncation floor and compare equal.  The claim is about the
  # calibration, so test the calibration.
  naive.scale <- (4 * B.y^2 / (n - ncol(X))) / mu.v
  expect_equal(syn$sigma.v, Delta / mu.v)
  expect_false(isTRUE(all.equal(syn$sigma.v, naive.scale)))
  # K_M > 1 on this leveraged design, so dropping it under-noises the channel
  expect_gt(Delta / (4 * B.y^2 / (n - ncol(X))), 1)
})
