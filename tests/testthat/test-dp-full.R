# V4: sensitive covariates.  The sound path is subsample-and-aggregate; the
# DFBETA local-sensitivity path is retired and reachable only as
# sens.method = "local", which warns.

# A public coefficient box and a public sigma^2 box for the fixtures below.
# These are facts about the DGP written in each test, fixed before any data are
# drawn.  They are does not read off d.
box.lo <- -5
box.hi <- 5

test_that("gen.syn.dp.full composes DP on X with DP on beta", {
  set.seed(1)
  d <- data.frame(
    y = rnorm(500),
    age = rnorm(500, 50, 10),
    sex = factor(sample(c("M", "F"), 500, replace = TRUE))
  )
  x.spec <- list(
    age = list(type = "continuous", bounds = c(0, 100)),
    sex = list(type = "categorical", levels = c("M", "F"))
  )
  # y is rnorm(500), independent of age and sex, so the population coefficients
  # are 0 and sigma^2 is 1.  Both boxes come from that DGP, not from d.
  syn <- gen.syn.dp.full(d, y ~ age + sex,
                          x.spec = x.spec,
                          epsilon = 10, delta = 1e-6,
                          clip.lo = box.lo, clip.hi = box.hi,
                          disp.lo = 0, disp.hi = 4, seed = 1)

  # structural
  expect_equal(syn$variant, "V4")
  expect_equal(syn$scope, "joint-release")
  expect_equal(nrow(syn$syn), 500)
  expect_equal(syn$sens.method, "subagg")

  # DP noise actually fired
  expect_gt(max(abs(syn$beta.dp - syn$beta.hat)), 0)

  # DP X marginals respect their declared bounds
  expect_true(all(syn$syn$age >= 0 & syn$syn$age <= 100))
  expect_true(all(syn$syn$sex %in% c("M", "F")))

  # OLS on (X_syn, y_syn) recovers beta.dp up to the small numerical slack from
  # the 1e-10 ridge regularisation on crossprod(X_syn)
  beta.analyst <- coef(lm(y ~ age + sex, syn$syn))
  v.common <- intersect(names(beta.analyst), names(syn$beta.dp))
  expect_lt(max(abs(beta.analyst[v.common] - syn$beta.dp[v.common])), 1e-5)
})

test_that("the three channels compose in ONE mu quadrature and hit the target delta", {
  set.seed(2)
  d <- data.frame(y = rnorm(400), age = rnorm(400, 50, 10),
                  income = rnorm(400, 100, 20))
  x.spec <- list(age    = list(type = "continuous", bounds = c(0, 100)),
                 income = list(type = "continuous", bounds = c(40, 160)))
  eps <- 10; delta <- 1e-6
  syn <- gen.syn.dp.full(d, y ~ age + income, x.spec = x.spec,
                          epsilon = eps, delta = delta,
                          clip.lo = box.lo, clip.hi = box.hi,
                          disp.lo = 0, disp.hi = 4, seed = 2)

  # mu^2 = mu_marg^2 + mu_corr^2 + mu_coef^2, exactly
  expect_equal(syn$mu.check, syn$mu.total)
  expect_equal(sqrt(sum(syn$mu.channels^2)), syn$mu.total)

  # and mu.total is the exact dual of (eps, delta): no slack, no accumulation
  expect_equal(synthexpln:::.gdp.delta(eps, syn$mu.total), delta,
               tolerance = 1e-9)

  # delta no longer accumulates across channels.  On the legacy path a copula
  # release reported delta.total = 2 * delta (two Gaussian mechanisms composed
  # by basic composition).  Under mu-GDP the whole release is calibrated to one
  # (eps, delta) through mu.
  expect_equal(syn$delta.total, delta)
  expect_true(syn$use.copula)
})

test_that("mu.split shifts budget between the X channels and the coefficient", {
  set.seed(21)
  d <- data.frame(y = rnorm(400), age = rnorm(400, 50, 10),
                  income = rnorm(400, 100, 20))
  x.spec <- list(age    = list(type = "continuous", bounds = c(0, 100)),
                 income = list(type = "continuous", bounds = c(40, 160)))
  args <- list(d = d, formula = y ~ age + income, x.spec = x.spec,
               epsilon = 10, delta = 1e-6, clip.lo = box.lo, clip.hi = box.hi,
               disp.lo = 0, disp.hi = 4, seed = 2)

  infer <- do.call(gen.syn.dp.full,
                   c(args, list(mu.split = c(marg = 0.05, corr = 0.05, coef = 0.90))))
  descr <- do.call(gen.syn.dp.full,
                   c(args, list(mu.split = c(marg = 0.40, corr = 0.40, coef = 0.20))))

  # the descriptive split buys X fidelity at the cost of coefficient precision
  expect_gt(descr$mu.channels[["marg"]], infer$mu.channels[["marg"]])
  expect_lt(descr$mu.channels[["coef"]], infer$mu.channels[["coef"]])
  # but the total privacy level is identical: the split moves budget, not cost
  expect_equal(descr$mu.total, infer$mu.total)
  expect_equal(descr$mu.check, descr$mu.total)

  expect_error(do.call(gen.syn.dp.full,
                       c(args, list(mu.split = c(marg = 0.5, corr = 0.5, coef = 0.5)))),
               "must sum to 1")
})

test_that("the public boxes are REQUIRED, not defaulted", {
  set.seed(3)
  d <- data.frame(y = rnorm(200), age = rnorm(200, 50, 10))
  x.spec <- list(age = list(type = "continuous", bounds = c(0, 100)))

  # This is the point.  The old API needed no bounds at all, which was the
  # symptom of the disease: with no public box there is no sensitivity bound and
  # therefore no theorem, so the call must not succeed.
  expect_error(gen.syn.dp.full(d, y ~ age, x.spec = x.spec, epsilon = 10),
               "clip.lo and clip.hi")
  expect_error(gen.syn.dp.full(d, y ~ age, x.spec = x.spec, epsilon = 10,
                               clip.lo = box.lo, clip.hi = box.hi),
               "disp.lo and disp.hi")
})

test_that("the retired DFBETA path is reachable but WARNS", {
  set.seed(4)
  d <- data.frame(y = rnorm(200), age = rnorm(200, 50, 10))
  x.spec <- list(age = list(type = "continuous", bounds = c(0, 100)))

  expect_warning(
    syn <- gen.syn.dp.full(d, y ~ age, x.spec = x.spec, epsilon = 10,
                           sens.method = "local", B.y = 6, seed = 2),
    "NO formal \\(eps, delta\\)-DP guarantee")
  expect_equal(syn$sens.method, "local")
  # it still needs a public B.y for its variance channel.  suppressWarnings
  # because the retirement warning fires before the refusal does, and it is the
  # refusal under test here.
  expect_error(
    suppressWarnings(
      gen.syn.dp.full(d, y ~ age, x.spec = x.spec, epsilon = 10,
                      sens.method = "local")),
    "B.y is required")
})

test_that("legacy path keeps the debiased residual-norm release and its eps split", {
  set.seed(41)
  n <- 300
  d <- data.frame(age = pmin(pmax(rnorm(n, 50, 10), 1), 99))
  d$y <- 1 + 0.1 * d$age + rnorm(n)
  x.spec <- list(age = list(type = "continuous", bounds = c(0, 100)))
  B.y <- 20
  suppressWarnings(
    syn <- gen.syn.dp.full(d, y ~ age, x.spec = x.spec,
                            epsilon = 10, delta = 1e-6,
                            sens.method = "local",
                            eps.split = c(X = 0.1, beta = 0.1, sigma = 0.8),
                            B.y = B.y, seed = 5))

  expect_equal(syn$epsilon.split[["X"]],     1, tolerance = 1e-8)
  expect_equal(syn$epsilon.split[["beta"]],  1, tolerance = 1e-8)
  expect_equal(syn$epsilon.split[["sigma"]], 8, tolerance = 1e-8)

  # the released quantity is the residual norm, exposed for audit
  expect_false(is.null(syn$resid.norm.dp))
  # sigma2.dp is its debiased square: b = 2 B.y / eps.sigma, E[L^2] = 2 b^2
  b <- 2 * B.y / syn$epsilon.split[["sigma"]]
  expect_equal(syn$sigma2.dp,
               max((syn$resid.norm.dp^2 - 2 * b^2) / (n - 2), 1e-6))
  g <- sqrt(sum(residuals(lm(y ~ age, d))^2))
  expect_lt(abs(syn$resid.norm.dp - g), 20 * b)

  # and it carries no mu accounting, because a Laplace channel cannot enter a
  # mu quadrature at all.  That is why the path is retired.
  expect_true(is.na(syn$mu.total))
  expect_true(is.na(syn$mu.check))
})

test_that("copula path preserves the inferential invariant and recovers X dependence", {
  set.seed(31)
  n  <- 3000
  x1 <- rnorm(n, 50, 10)
  x2 <- 0.8 * (x1 - 50) + rnorm(n, 0, 6) + 100        # correlated with x1
  d  <- data.frame(y = 2 + 0.5 * x1 - 0.3 * x2 + rnorm(n),
                   age = x1, income = x2)
  x.spec <- list(age    = list(type = "continuous", bounds = c(0, 100)),
                 income = list(type = "continuous", bounds = c(40, 160)))
  # a descriptive-fidelity split: this release wants X structure, so it buys it
  ms <- c(marg = 0.25, corr = 0.25, coef = 0.50)

  syn <- gen.syn.dp.full(d, y ~ age + income, x.spec = x.spec,
                          epsilon = 10, delta = 1e-6, mu.split = ms,
                          clip.lo = box.lo, clip.hi = box.hi,
                          disp.lo = 0, disp.hi = 4,
                          use.copula = TRUE, seed = 3)

  # invariant: refit recovers beta.dp exactly (projection, not X_syn's law)
  b.hat <- coef(lm(y ~ age + income, syn$syn))
  v.common <- intersect(names(b.hat), names(syn$beta.dp))
  expect_lt(max(abs(b.hat[v.common] - syn$beta.dp[v.common])), 1e-5)

  # the copula recovers the sign and approximate magnitude of the X correlation, which the
  # independent path cannot
  syn0 <- gen.syn.dp.full(d, y ~ age + income, x.spec = x.spec,
                           epsilon = 10, delta = 1e-6, mu.split = ms,
                           clip.lo = box.lo, clip.hi = box.hi,
                           disp.lo = 0, disp.hi = 4,
                           use.copula = FALSE, seed = 3)
  rho.true <- cor(d$age, d$income)
  err.cop  <- abs(cor(syn$syn$age,  syn$syn$income)  - rho.true)
  err.ind  <- abs(cor(syn0$syn$age, syn0$syn$income) - rho.true)
  expect_lt(err.cop, err.ind)

  expect_true(syn$use.copula)
  expect_false(is.null(syn$dp.corr))

  # with the copula OFF the correlation share is folded back into the marginals
  # in quadrature, so the total mu is preserved either way
  expect_equal(syn0$mu.total, syn$mu.total)
  expect_equal(syn0$mu.check, syn0$mu.total)
  expect_equal(syn0$mu.channels[["corr"]], 0)
})

test_that("score.method does not change V4 coefficient recovery", {
  set.seed(33)
  n  <- 2000
  x1 <- rnorm(n, 50, 10); x2 <- 0.8 * (x1 - 50) + rnorm(n, 0, 6) + 100
  d  <- data.frame(y = 2 + 0.5 * x1 - 0.3 * x2 + rnorm(n), age = x1, income = x2)
  x.spec <- list(age    = list(type = "continuous", bounds = c(0, 100)),
                 income = list(type = "continuous", bounds = c(40, 160)))
  gen <- function(sm)
    gen.syn.dp.full(d, y ~ age + income, x.spec = x.spec, epsilon = 10,
                    delta = 1e-6, clip.lo = box.lo, clip.hi = box.hi,
                    disp.lo = 0, disp.hi = 4, use.copula = TRUE,
                    score.method = sm, seed = 3)
  a  <- gen("box")
  b  <- gen("normal")
  cc <- gen("adaptive")

  # the copula does not touch the coefficient release: identical DP betas
  expect_equal(a$beta.dp, b$beta.dp)
  expect_equal(a$beta.dp, cc$beta.dp)
  # but score.method does change the X law
  expect_false(isTRUE(all.equal(cor(a$syn$age, a$syn$income),
                                cor(b$syn$age, b$syn$income))))
})

test_that("use.copula = FALSE and single-continuous inputs keep old behaviour", {
  set.seed(32)
  d <- data.frame(y = rnorm(400), age = rnorm(400, 50, 10),
                  sex = factor(sample(c("M", "F"), 400, replace = TRUE)))
  x.spec <- list(age = list(type = "continuous", bounds = c(0, 100)),
                 sex = list(type = "categorical", levels = c("M", "F")))
  # only one continuous column: copula is a no-op regardless of the flag
  syn <- gen.syn.dp.full(d, y ~ age + sex, x.spec = x.spec, epsilon = 10,
                          clip.lo = box.lo, clip.hi = box.hi,
                          disp.lo = 0, disp.hi = 4,
                          use.copula = TRUE, seed = 4)
  expect_equal(syn$use.copula, FALSE)          # auto-disabled, p_c < 2
  expect_null(syn$dp.corr)
  expect_equal(syn$delta.total, syn$delta, tolerance = 1e-12)
})
