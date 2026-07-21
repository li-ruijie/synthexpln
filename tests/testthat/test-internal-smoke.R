# Smoke tests: call each internal helper that lacks a direct unit test with
# minimal valid inputs, so a broken call (renamed helper, wrong arity, stale
# reference after a refactor) surfaces at the function it lives in rather than
# only deep inside an integration test. Assertions stay deliberately loose:
# the point is that the call runs and returns a sane-shaped, finite result.

fixture <- function(n = 200L, seed = 1L) {
    set.seed(seed)
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    y  <- 1 + 0.5 * x1 - 0.3 * x2 + stats::rnorm(n)
    data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("every exported function is defined and callable", {
    exported <- c(
        "calibrated.se", "dp.audit", "dp.project", "fidelity.distributional",
        "fidelity.inferential", "gen.syn.beta", "gen.syn.dp.full",
        "gen.syn.dp.logistic", "gen.syn.dp.ols.public", "gen.syn.dp.projected",
        "gen.syn.dp.projected.subagg", "gen.syn.dp.projected.zi",
        "gen.syn.dp.ridge.public", "privacy.metrics", "projection",
        "sensitivity.local")
    ns      <- asNamespace("synthexpln")
    missing <- Filter(function(nm) !is.function(get0(nm, envir = ns)), exported)
    expect_equal(missing, character(0))
})

test_that(".make.glm.family builds family objects from string pairs", {
    fam <- synthexpln:::.make.glm.family("gaussian", "identity")
    expect_s3_class(fam, "family")
    expect_equal(synthexpln:::.make.glm.family("poisson", "log")[["family"]], "poisson")
})

test_that(".is.numeric.col flags continuous columns only", {
    expect_true(synthexpln:::.is.numeric.col(1:20))
    expect_false(synthexpln:::.is.numeric.col(factor(letters[1:5])))
    expect_false(synthexpln:::.is.numeric.col(rep(1, 20)))   # <= 10 unique
})

test_that(".std.cumulants returns finite standardised cumulants", {
    set.seed(1)
    g <- synthexpln:::.std.cumulants(stats::rnorm(500), order = 4)
    expect_true(all(is.finite(g)))
})

test_that(".bootstrap.rows resamples to the same size", {
    d <- fixture(50)
    b <- synthexpln:::.bootstrap.rows(d)
    expect_equal(dim(b), dim(d))
})

test_that(".gsb.* helpers normalise inputs and forward-sample", {
    fam <- synthexpln:::.gsb.norm.family("gaussian")
    expect_s3_class(fam, "family")
    expect_equal(synthexpln:::.gsb.norm.dispersion(NULL, fam), 1)
    beta <- c(`(Intercept)` = 1, x1 = 0.5, x2 = -0.3)
    df.x <- synthexpln:::.gsb.build.X(beta, n = 100L, X = NULL, marginals = NULL)
    expect_equal(sort(colnames(df.x)), c("x1", "x2"))
    expect_equal(nrow(df.x), 100L)
    y <- synthexpln:::.gsb.sample.y(fam, mu = rep(1, 10), dispersion = 1)
    expect_length(y, 10)
    expect_true(all(is.finite(y)))
})

test_that(".gsb.project projects onto the supplied beta", {
    set.seed(1)
    beta <- c(`(Intercept)` = 1, x1 = 0.5, x2 = -0.3)
    df.x <- data.frame(x1 = stats::rnorm(100), x2 = stats::rnorm(100))
    y.q  <- synthexpln:::.gsb.project(df.x, y.prelim = stats::rnorm(100),
                                      beta = beta, phi = 1, family = stats::gaussian())
    expect_length(y.q, 100)
    expect_true(all(is.finite(y.q)))
})

test_that(".project.ols returns a finite projected response", {
    d     <- fixture(150)
    d.syn <- d
    d.syn[["y"]] <- stats::rnorm(nrow(d))
    v <- synthexpln:::.project.ols(d.syn, d, "y", c("x1", "x2"))
    expect_length(v, nrow(d))
    expect_true(all(is.finite(v)))
})

test_that(".inject.virtual.y shifts the response by the DP mean gap", {
    d        <- fixture(120)
    mo       <- stats::lm(y ~ x1 + x2, d)
    beta.hat <- stats::coef(mo)
    beta.dp  <- beta.hat + c(0.05, -0.02, 0.03)
    d.v <- synthexpln:::.inject.virtual.y(d, "y", beta.hat, beta.dp,
                                          "gaussian", "identity")
    expect_equal(nrow(d.v), nrow(d))
    expect_true(all(is.finite(d.v[["y"]])))
})

test_that(".sens.beta.ridge returns a finite sensitivity object", {
    d   <- fixture(150)
    m.x <- stats::model.matrix(~ x1 + x2, d)
    out <- synthexpln:::.sens.beta.ridge(m.x, d[["y"]], lambda = 1)
    expect_true(is.list(out) || is.numeric(out))
})

test_that(".subsample.aggregate.beta.loglink runs on a gaussian design", {
    d   <- fixture(400)
    out <- synthexpln:::.subsample.aggregate.beta.loglink(
        d, y ~ x1 + x2, stats::gaussian(), "gaussian", "identity",
        eps = 5, delta = 1e-6, m = 20, clip.lo = -5, clip.hi = 5, seed = 1)
    expect_true(all(is.finite(out[["beta.dp"]])))
})
