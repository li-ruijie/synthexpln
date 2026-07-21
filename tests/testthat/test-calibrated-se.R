test_that("calibrated.se matches the closed form and inflates the naive SE", {
    se.syn     <- c(0.10, 0.20, 0.05)
    sigma.beta <- 0.08
    s.j        <- c(0.3, 0.3, 0.3)
    out <- calibrated.se(se.syn, sigma.beta, s.j)
    expect_equal(out, sqrt(se.syn^2 + (sigma.beta / s.j)^2))
    expect_true(all(out >= se.syn))                 # DP noise only inflates
    expect_true(all(is.finite(out) & out > 0))
})

test_that("calibrated.se returns the naive SE when the DP noise term is zero", {
    se.syn <- c(0.1, 0.2)
    expect_equal(calibrated.se(se.syn, 0, c(0.3, 0.3)), se.syn)
})

test_that("s.j defaults to 1, because the sound routes carry no column scale", {
    # The correction is sqrt(SE^2 + sigma_beta^2).  There is no column scale in
    # it.  A divisor is right only when the mechanism perturbed the standardised
    # coefficients, and the global-sensitivity routes do not: they add isotropic
    # noise to the raw beta.
    se.syn <- c(0.1, 0.2)
    expect_equal(calibrated.se(se.syn, 0.05),
                 sqrt(se.syn^2 + 0.05^2))
    expect_equal(calibrated.se(se.syn, 0.05, s.j = 1),
                 calibrated.se(se.syn, 0.05))
})

test_that("the calibrated SE uses the RAW-scale noise, and cannot be shrunk by a stale column scale", {
    # the regression guard.  income has a column sd in the hundreds.  If the
    # correction divides sigma.beta by that scale, the DP variance component is
    # shrunk by the same factor and the corrected interval is barely
    # wider than the uncorrected one.  That understates the uncertainty of a
    # private release, which is the only direction this error must not take.
    set.seed(1)
    n <- 400
    d <- data.frame(age = runif(n, 20, 80), income = runif(n, 0, 4000))
    d$y <- 10 + 0.05 * d$age + 0.001 * d$income + rnorm(n)

    syn <- gen.syn.dp.ols.public(d, y ~ age + income, epsilon = 5,
                                 delta = 1e-6, B.y = 40, seed = 1)

    # the mechanism adds isotropic noise on the raw scale, so every coefficient
    # carries the same raw sd
    expect_equal(unname(syn$sigma.beta.raw),
                 rep(syn$sigma.beta, 3))
    expect_named(syn$sigma.beta.raw,
                 c("(Intercept)", "age", "income"))
    # and there is no column scale on the object to divide by
    expect_null(syn$s.j)

    m <- fidelity.inferential(syn, orig = d, formula = y ~ age + income,
                              calibrate = TRUE)
    # the DP correction actually FIRES: the interval widens, and by the right
    # amount, on every coefficient including the wide-scale one
    expect_true(all(m$se.corr > m$se.syn))
    expect_equal(unname(m$se.corr),
                 unname(sqrt(m$se.syn^2 + syn$sigma.beta^2)))

    # what the stale column-scale division would have produced: a correction
    # shrunk by income's column sd, hence invisible
    s.income <- sd(d$income)
    stale    <- sqrt(m$se.syn[["income"]]^2 + (syn$sigma.beta / s.income)^2)
    expect_gt(m$se.corr[["income"]], 10 * stale)
})

test_that("calibration REFUSES rather than silently reporting an uncorrected SE", {
    # The old code asked for $sigma.beta and $s.j and, if either was missing,
    # quietly returned se.corr = se.syn "because there is no DP noise to correct
    # for".  Once the sound routes stopped carrying a column scale that branch
    # began firing on releases that do carry noise.  An SE that is too small is
    # an inferential claim the release cannot support, so it must fail loudly.
    set.seed(2)
    d <- data.frame(x = rnorm(200))
    d$y <- 1 + 0.5 * d$x + rnorm(200)
    syn <- gen.syn.dp.ols.public(d, y ~ x, epsilon = 5, B.y = 10, seed = 1)

    # a DP object whose noise scale does not name every coefficient must ERROR
    broken <- syn
    broken$sigma.beta.raw <- syn$sigma.beta.raw["(Intercept)"]
    expect_error(
      fidelity.inferential(broken, orig = d, formula = y ~ x, calibrate = TRUE),
      "does not name every coefficient")

    # a non-DP projection genuinely has no noise, and that is the only case that
    # may fall back to the uncorrected SE
    plain <- projection(d, y ~ x)
    expect_message(
      m <- fidelity.inferential(plain, orig = d, formula = y ~ x,
                                calibrate = TRUE),
      "carries no DP noise")
    expect_equal(m$se.corr, m$se.syn)
})

test_that("every DP generator reports a named raw-scale noise sd", {
    # calibrated.se() must not have to guess which scale it is looking at, so
    # the contract is: every DP object names its own noise, per coefficient, on
    # the raw scale.
    set.seed(3)
    n <- 800
    d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
    d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
    db <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
    db$y <- rbinom(n, 1, plogis(-0.3 + 0.8 * db$x1))
    x.spec <- list(x1 = list(type = "continuous", bounds = c(-5, 5)),
                   x2 = list(type = "continuous", bounds = c(-5, 5)))
    disp <- dp.dispersion.public(1, source = "declared")

    objs <- list(
      V1 = suppressWarnings(gen.syn.dp.projected(d, y ~ x1 + x2, epsilon = 5)),
      V2 = gen.syn.dp.ols.public(d, y ~ x1 + x2, epsilon = 5, B.y = 20, seed = 1),
      V3 = gen.syn.dp.ridge.public(d, y ~ x1 + x2, lambda = 1, epsilon = 5,
                                   B.y = 20, seed = 1),
      V4 = gen.syn.dp.full(d, y ~ x1 + x2, x.spec = x.spec, epsilon = 5,
                           clip.lo = -5, clip.hi = 5, disp.lo = 0, disp.hi = 4,
                           seed = 1),
      V5 = gen.syn.dp.logistic(db, y ~ x1 + x2, epsilon = 5, clip.lo = -3,
                               clip.hi = 3, seed = 1),
      D2 = gen.syn.dp.projected.subagg(d, y ~ x1 + x2, epsilon = 5,
                                       family = gaussian(), clip.lo = -5,
                                       clip.hi = 5, dispersion = disp, seed = 1))

    for (nm in names(objs)) {
      o <- objs[[nm]]
      raw <- o$sigma.beta.raw
      expect_false(is.null(raw), info = nm)
      expect_setequal(names(raw), names(o$beta.dp))
      expect_true(all(is.finite(raw) & raw > 0), info = nm)
      # and fidelity.inferential will actually use it
      expect_equal(synthexpln:::.dp.noise.sd(o, names(o$beta.dp)),
                   raw[names(o$beta.dp)], info = nm)
    }
})
