# test-dp-soundness.R -- the properties the rest of the suite cannot fail on.
#
# WHY this FILE exists.  Before 0.2.0 this package had 229 passing tests while
# every one of its five DP generators calibrated its noise to the DFBETA local
# sensitivity, and its soundest route retained the real response in the release.
# The tests checked the ARITHMETIC of each mechanism against its own assumed
# sensitivity, with real rigour, and not once asked whether that sensitivity
# was a bound at all.  A green suite coexisted with mechanisms that had no valid
# calibration.
#
# These tests therefore assert the missing property from four directions.  Each is
# written to FAIL on the code as it stood before 0.2.0.
#
#   1. The noise calibration must be NEIGHBOUR-INVARIANT.  This is the
#      executable form of the hypothesis every DP theorem here needs: one bound,
#      binding both neighbouring datasets.
#   2. The release must NOT echo the empirical residual distribution.  (In
#      test-dp-subagg.R, where the fixture lives.)
#   3. A "public" constant must not be a sample STATISTIC.  Test 1 cannot catch
#      this: a bound hardcoded from a sample fit is a constant in the source,
#      hence trivially neighbour-invariant, and still a function of the
#      sensitive data.  The tell is statistical.
#   4. EPSILON must BITE.  If a channel reaches the release without passing
#      through the noise, reducing epsilon does not protect it, so a quantity
#      flat in epsilon is the fingerprint of an unprotected channel.  This is the
#      test that catches the defects nobody has thought of yet.

# ---------------------------------------------------------------------------
# 1. Neighbour-invariance of the noise calibration
# ---------------------------------------------------------------------------

test_that("SOUNDNESS 1: every global sensitivity is neighbour-invariant", {
  set.seed(101)
  n <- 300
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n, 0, 5))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  B.y <- 20

  # a neighbouring dataset: one row's RESPONSE replaced, adversarially
  d.nb <- d
  d.nb$y[1] <- B.y

  m.X    <- model.matrix(y ~ x1 + x2, d)
  m.X.nb <- model.matrix(y ~ x1 + x2, d.nb)   # X is public and unchanged

  # The bound must be identical on both neighbours.  It is a function of the
  # public design and the public bound, and of nothing else.
  expect_equal(synthexpln:::.gs.beta.ols.public(m.X, B.y),
               synthexpln:::.gs.beta.ols.public(m.X.nb, B.y))
  expect_equal(synthexpln:::.gs.beta.ridge.public(m.X, 1, B.y),
               synthexpln:::.gs.beta.ridge.public(m.X.nb, 1, B.y))
  expect_equal(synthexpln:::.gs.sigma2.public(m.X, B.y),
               synthexpln:::.gs.sigma2.public(m.X.nb, B.y))
  expect_equal(synthexpln:::.gs.sigma2.ridge(m.X, 1, B.y),
               synthexpln:::.gs.sigma2.ridge(m.X.nb, 1, B.y))

  # The noise scale the generators actually use must be identical too.  This
  # is the property, at the level that matters.
  s.a <- gen.syn.dp.ols.public(d,    y ~ x1 + x2, epsilon = 5, B.y = B.y, seed = 1)
  s.b <- gen.syn.dp.ols.public(d.nb, y ~ x1 + x2, epsilon = 5, B.y = B.y, seed = 1)
  expect_equal(s.a$sigma.beta, s.b$sigma.beta)
  expect_equal(s.a$sigma.v,    s.b$sigma.v)

  r.a <- gen.syn.dp.ridge.public(d,    y ~ x1 + x2, lambda = 1, epsilon = 5,
                                 B.y = B.y, seed = 1)
  r.b <- gen.syn.dp.ridge.public(d.nb, y ~ x1 + x2, lambda = 1, epsilon = 5,
                                 B.y = B.y, seed = 1)
  expect_equal(r.a$sigma.beta, r.b$sigma.beta)
  expect_equal(r.a$sigma.v,    r.b$sigma.v)
})

test_that("SOUNDNESS 1: the RETIRED local sensitivity MOVES between neighbours", {
  # The positive control.  If this ever stops moving, the test above has become
  # vacuous and is no longer checking anything.  Local sensitivity is a function
  # of the realised data: that is not a subtlety, it is the definition, and it is
  # exactly why no DP theorem can be instantiated with it.
  set.seed(102)
  n <- 300
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n, 0, 5))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  d.nb <- d
  d.nb$y[1] <- 20

  ls.a <- sensitivity.local(lm(y ~ x1 + x2, d))
  ls.b <- sensitivity.local(lm(y ~ x1 + x2, d.nb))
  expect_false(isTRUE(all.equal(ls.a, ls.b)))

  # the ridge DFBETA likewise
  m.X <- model.matrix(y ~ x1 + x2, d)
  rl.a <- synthexpln:::.sens.beta.ridge(m.X, d$y,    1)$LS
  rl.b <- synthexpln:::.sens.beta.ridge(m.X, d.nb$y, 1)$LS
  expect_false(isTRUE(all.equal(rl.a, rl.b)))
})

test_that("SOUNDNESS 1: the subagg per-coordinate sensitivity is data-independent", {
  # subsample-and-aggregate takes its sensitivity from the public box and the
  # block count, s_j = (hi_j - lo_j) / m.  No statistic of the data appears in
  # it, so it holds under full-row replacement: x, y, or both.
  set.seed(103)
  n <- 600
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)

  # a full-row neighbour: both covariates and the response moved, to the corners
  d.nb <- d
  d.nb[1, ] <- list(x1 = 5, x2 = -5, y = 50)

  disp <- dp.dispersion.public(1, source = "unit residual variance, by the DGP")
  gen <- function(dd) gen.syn.dp.projected.subagg(
    dd, y ~ x1 + x2, epsilon = 5, m = 20, family = gaussian(),
    clip.lo = -5, clip.hi = 5, dispersion = disp, seed = 1)

  a <- gen(d); b <- gen(d.nb)
  expect_equal(a$sigma,   b$sigma)
  expect_equal(a$sigma.Q, b$sigma.Q)
  expect_equal(a$s.j,     b$s.j)

  # and the bound is what the theorem claims: the aggregate moves by at most
  # s_j when one row changes.  Measured against the WORST case above.
  moved <- max(abs(a$beta.agg - b$beta.agg))
  expect_lte(moved, max(a$s.j) + 1e-9)
})

# ---------------------------------------------------------------------------
# 3. A "public" constant must not be a sample statistic
# ---------------------------------------------------------------------------

test_that("SOUNDNESS 3: no generator reads its noise scale off the sample", {
  # The bound must not track the data.  Inflate the response by a factor of 10
  # and hold the declared public bound fixed: a sound calibration does not move,
  # because B.y is an assumption about the universe, not a measurement of d.
  #
  # The removed default was B.y <- 1.5 * max(abs(y)), which would scale exactly
  # with the data below and so make the noise a function of what it protects.
  set.seed(104)
  n <- 300
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  d.big <- transform(d, y = y * 10)
  B.y <- 100                                  # public, fixed, covers both

  a <- gen.syn.dp.ols.public(d,     y ~ x1 + x2, epsilon = 5, B.y = B.y, seed = 1)
  b <- gen.syn.dp.ols.public(d.big, y ~ x1 + x2, epsilon = 5, B.y = B.y, seed = 1)
  expect_equal(a$sigma.beta, b$sigma.beta)
  expect_equal(a$sigma.v,    b$sigma.v)

  # the retired default would have moved by ~10x.  Assert that it would, so the
  # test above is known to be discriminating and not merely quiet.
  expect_equal(1.5 * max(abs(d.big$y)) / (1.5 * max(abs(d$y))), 10,
               tolerance = 0.02)

  # subagg likewise: its scale comes from the box, not from y
  disp <- dp.dispersion.public(1, source = "declared")
  sa <- function(dd) gen.syn.dp.projected.subagg(
    dd, y ~ x1 + x2, epsilon = 5, family = gaussian(),
    clip.lo = -5, clip.hi = 5, dispersion = disp, seed = 1)$sigma
  expect_equal(sa(d), sa(d.big))
})

test_that("SOUNDNESS 3: the dispersion constructors demand their provenance", {
  # `source` is the load-bearing part of the design.  The unsound value still
  # exists as a construction, but to use it you have to type a sentence that
  # indicts it, which no one will do and which any reader of the code will catch.
  expect_error(dp.dispersion.public(4), "source is REQUIRED")
  expect_error(dp.dispersion.private(0, 9), "source is REQUIRED")
  expect_error(dp.dispersion.public(4, source = ""), "source is REQUIRED")
  expect_error(dp.dispersion.public(-1, source = "x"), "positive")
  expect_error(dp.dispersion.private(9, 0, source = "x"), "hi > lo")

  # the defect now has to look like this, and it is self-indicting
  d <- data.frame(x = rnorm(50)); d$y <- d$x + rnorm(50)
  bad <- dp.dispersion.public(var(residuals(lm(y ~ x, d))),
                              source = "the sample")
  expect_equal(bad$mode, "public")
  expect_match(bad$source, "sample")
})

# ---------------------------------------------------------------------------
# 4. Epsilon must bite
# ---------------------------------------------------------------------------

test_that("SOUNDNESS 4: EPSILON BITES on every released channel", {
  # If a channel reaches the release without passing through the noise,
  # reducing epsilon does not protect it, and the fingerprint of that is a
  # quantity flat in epsilon.  Sweep epsilon over two orders of magnitude and
  # require every channel's noise to fall monotonically.
  set.seed(105)
  n <- 500
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  v.eps <- c(0.1, 0.5, 1, 5, 10, 50)

  # V2: both channels
  m2 <- vapply(v.eps, function(e) {
    s <- gen.syn.dp.ols.public(d, y ~ x1 + x2, epsilon = e, B.y = 20, seed = 1)
    c(beta = s$sigma.beta, v = s$sigma.v)
  }, numeric(2))
  expect_true(all(diff(m2["beta", ]) < 0))
  expect_true(all(diff(m2["v", ]) < 0))

  # V3: both channels
  m3 <- vapply(v.eps, function(e) {
    s <- gen.syn.dp.ridge.public(d, y ~ x1 + x2, lambda = 1, epsilon = e,
                                 B.y = 20, seed = 1)
    c(beta = s$sigma.beta, v = s$sigma.v)
  }, numeric(2))
  expect_true(all(diff(m3["beta", ]) < 0))
  expect_true(all(diff(m3["v", ]) < 0))

  # Route 2: the coefficient channel
  disp <- dp.dispersion.public(1, source = "declared")
  v.r2 <- vapply(v.eps, function(e)
    gen.syn.dp.projected.subagg(d, y ~ x1 + x2, epsilon = e, family = gaussian(),
                                clip.lo = -5, clip.hi = 5, dispersion = disp,
                                seed = 1)$sigma[[1]], numeric(1))
  expect_true(all(diff(v.r2) < 0))

  # V4: all three channels, including the marginal and correlation releases
  x.spec <- list(x1 = list(type = "continuous", bounds = c(-5, 5)),
                 x2 = list(type = "continuous", bounds = c(-5, 5)))
  m4 <- vapply(v.eps, function(e) {
    s <- gen.syn.dp.full(d, y ~ x1 + x2, x.spec = x.spec, epsilon = e,
                         clip.lo = -5, clip.hi = 5, disp.lo = 0, disp.hi = 4,
                         seed = 1)
    s$mu.channels
  }, numeric(3))
  # mu RISES with epsilon on every channel (mu is signal-to-noise, so noise falls)
  expect_true(all(diff(m4["marg", ]) > 0))
  expect_true(all(diff(m4["corr", ]) > 0))
  expect_true(all(diff(m4["coef", ]) > 0))
})

test_that("SOUNDNESS 4: mu-GDP composition is exact, not merely reported", {
  # Every generator REPORTS a mu.total.  Assert it is the real one: the channels
  # must re-compose in quadrature to exactly that number, and that number must be
  # the exact dual of the (eps, delta) the caller asked for.  A reported mu that
  # nothing checks is how a stated guarantee drifts from a delivered one.
  set.seed(106)
  n <- 400
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  x.spec <- list(x1 = list(type = "continuous", bounds = c(-5, 5)),
                 x2 = list(type = "continuous", bounds = c(-5, 5)))

  for (eps in c(0.5, 1, 4, 10)) {
    delta <- 1e-6
    mu.star <- synthexpln:::.gdp.mu.from.eps.delta(eps, delta)

    s2 <- gen.syn.dp.ols.public(d, y ~ x1 + x2, epsilon = eps, delta = delta,
                                B.y = 20, seed = 1)
    expect_equal(s2$mu.total, mu.star)
    expect_equal(sqrt(s2$mu.beta^2 + s2$mu.v^2), mu.star)
    expect_equal(synthexpln:::.gdp.delta(eps, s2$mu.total), delta,
                 tolerance = 1e-9)

    s4 <- gen.syn.dp.full(d, y ~ x1 + x2, x.spec = x.spec, epsilon = eps,
                          delta = delta, clip.lo = -5, clip.hi = 5,
                          disp.lo = 0, disp.hi = 4, seed = 1)
    expect_equal(s4$mu.total, mu.star)
    expect_equal(s4$mu.check, mu.star)
    expect_equal(synthexpln:::.gdp.delta(eps, s4$mu.total), delta,
                 tolerance = 1e-9)
    # delta does not accumulate across the three channels
    expect_equal(s4$delta.total, delta)
  }
})

test_that("SOUNDNESS 4: the classical multiplier UNDER-delivers delta above eps = 1", {
  # The reason .gdp.mu.from.eps.delta exists rather than the classical bound.
  # At eps = 10, delta = 1e-6 the Dwork-Roth Thm A.1 multiplier certifies only
  # delta = 1.9e-6, nearly twice what was asked for.  It is valid only for
  # eps < 1, and a mechanism calibrated by it above that silently under-noises.
  eps <- 10; delta <- 1e-6

  mu.classical <- 1 / (sqrt(2 * log(1.25 / delta)) / eps)
  delta.got    <- synthexpln:::.gdp.delta(eps, mu.classical)
  expect_gt(delta.got, delta)
  expect_equal(delta.got, 1.9e-6, tolerance = 0.1)

  # the exact dual delivers exactly what was asked.  uniroot's tol is 1e-12 on
  # MU, so the relative error it induces on delta is ~1e-12 rather than 1e-15,
  # which makes 1e-9 the defensible tolerance to claim here.
  mu.exact <- synthexpln:::.gdp.mu.from.eps.delta(eps, delta)
  expect_equal(synthexpln:::.gdp.delta(eps, mu.exact), delta, tolerance = 1e-9)
  expect_lt(mu.exact, mu.classical)     # exact => more noise => a real guarantee

  # below eps = 1 the classical bound is valid, and the package's multiplier
  # uses it there
  expect_equal(synthexpln:::.dp.gm.multiplier(0.5, delta),
               sqrt(2 * log(1.25 / delta)) / 0.5)
})

# ---------------------------------------------------------------------------
# The DP surface, enumerated: no generator may reach dfbeta on its sound path
# ---------------------------------------------------------------------------

test_that("SOUNDNESS: the sound paths do not call stats::dfbeta", {
  # The structural guard.  Mask dfbeta with a hard error and run every sound
  # path.  If any of them reaches the retired local sensitivity, it dies here.
  #
  # This would have caught the original defect in one line, and it is the test
  # the pre-0.2.0 suite could not have written, because every generator would
  # have failed it.
  set.seed(107)
  n <- 600
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  db <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  db$y <- rbinom(n, 1, plogis(-0.3 + 0.8 * db$x1))
  x.spec <- list(x1 = list(type = "continuous", bounds = c(-5, 5)),
                 x2 = list(type = "continuous", bounds = c(-5, 5)))
  disp <- dp.dispersion.public(1, source = "declared")

  with_mocked_bindings(
    dfbeta = function(...) stop("RETIRED: a sound path must not reach dfbeta"),
    .package = "stats",
    {
      expect_s3_class(
        gen.syn.dp.ols.public(d, y ~ x1 + x2, epsilon = 5, B.y = 20, seed = 1),
        "synthexpln")
      expect_s3_class(
        gen.syn.dp.ridge.public(d, y ~ x1 + x2, lambda = 1, epsilon = 5,
                                B.y = 20, seed = 1),
        "synthexpln")
      expect_s3_class(
        gen.syn.dp.projected.subagg(d, y ~ x1 + x2, epsilon = 5,
                                    family = gaussian(), clip.lo = -5,
                                    clip.hi = 5, dispersion = disp, seed = 1),
        "synthexpln")
      expect_s3_class(
        gen.syn.dp.full(d, y ~ x1 + x2, x.spec = x.spec, epsilon = 5,
                        clip.lo = -5, clip.hi = 5, disp.lo = 0, disp.hi = 4,
                        seed = 1),
        "synthexpln")
      expect_s3_class(
        gen.syn.dp.logistic(db, y ~ x1 + x2, epsilon = 5, clip.lo = -3,
                            clip.hi = 3, seed = 1),
        "synthexpln")
    })
})
