test_that("gen.syn.beta() runs with only `beta` supplied", {
    set.seed(1L)
    o <- gen.syn.beta(c(1, -0.5, 2))
    expect_s3_class(o, "synthexpln")
    expect_equal(o$variant, "from-beta")
    expect_equal(nrow(o$syn), 1000L)
    expect_named(o$beta, c("(Intercept)", "x1", "x2"))
    # refit should be in the right neighbourhood
    b.refit <- coef(lm(o$formula, o$syn))
    expect_equal(unname(b.refit), c(1, -0.5, 2), tolerance = 0.1)
})

test_that("gen.syn.beta() preserves named beta", {
    o <- gen.syn.beta(c(`(Intercept)` = 0.5, age = 0.1, sex = -1),
                      seed = 1L)
    expect_named(o$beta, c("(Intercept)", "age", "sex"))
    expect_true(all(c("age", "sex") %in% colnames(o$syn)))
})

test_that("gen.syn.beta() accepts family as a string with link", {
    o <- gen.syn.beta(c(-1, 0.3, 0.7), family = "poisson/log",
                      n = 500L, seed = 3L)
    expect_equal(o$family$family, "poisson")
    expect_equal(o$family$link,   "log")
    expect_true(all(o$syn$.y == as.integer(o$syn$.y)))
    expect_true(all(o$syn$.y >= 0))
})

test_that("gen.syn.beta() accepts a bare family string", {
    o <- gen.syn.beta(c(0, 1.2), family = "binomial", n = 200L, seed = 4L)
    expect_equal(o$family$family, "binomial")
    expect_equal(o$family$link,   "logit")
    expect_true(all(o$syn$.y %in% c(0, 1)))
})

test_that("gen.syn.beta() builds X from a named list of marginals", {
    o <- gen.syn.beta(
        beta      = c(`(Intercept)` = 30, age = 0.4, female = -2.1),
        marginals = list(age    = function(n) stats::rnorm(n, 50, 12),
                         female = function(n) stats::rbinom(n, 1L, 0.5)),
        seed = 5L
    )
    expect_true(mean(o$syn$age) > 40 && mean(o$syn$age) < 60)
    expect_true(all(o$syn$female %in% c(0, 1)))
})

test_that("gen.syn.beta() bootstraps marginals from a data.frame", {
    src <- data.frame(x1 = rnorm(50), x2 = rbinom(50, 1L, 0.3))
    o <- gen.syn.beta(c(0, 1, -1), n = 500L, marginals = src, seed = 6L)
    expect_equal(nrow(o$syn), 500L)
    # bootstrap should reproduce values from the source set
    expect_true(all(o$syn$x2 %in% src$x2))
})

test_that("gen.syn.beta() uses user-supplied X verbatim", {
    set.seed(6489L)
    X.in <- data.frame(x1 = rnorm(200), x2 = rnorm(200))
    o <- gen.syn.beta(c(0, 1, -1), X = X.in, seed = 6L)
    expect_equal(nrow(o$syn), 200L)
    expect_equal(o$syn$x1, X.in$x1)
})

test_that("gen.syn.beta() X overrides marginals with a message", {
    expect_message(
        gen.syn.beta(c(0, 1),
                     X = data.frame(x1 = rnorm(50)),
                     marginals = list(x1 = function(n) rep(0, n)),
                     seed = 1L),
        "ignored"
    )
})

test_that("project mode recovers beta exactly for gaussian/identity", {
    beta.tgt <- c(`(Intercept)` = 2, x1 = 0.8, x2 = -0.3)
    o <- gen.syn.beta(beta.tgt, method = "project", seed = 7L)
    b.refit <- coef(lm(o$formula, o$syn))
    expect_equal(unname(b.refit), unname(beta.tgt), tolerance = 1e-10)
})

test_that("project mode recovers beta to Newton tolerance for gamma/log", {
    beta.tgt <- c(`(Intercept)` = 1.5, x1 = 0.4, x2 = -0.6)
    o <- gen.syn.beta(beta.tgt, family = "gamma/log", dispersion = 0.5,
                      method = "project", n = 2000L, seed = 9L)
    b.refit <- coef(glm(o$formula, data = o$syn,
                        family = Gamma(link = "log")))
    expect_equal(unname(b.refit), unname(beta.tgt), tolerance = 1e-5)
})

test_that("project mode recovers beta within rounding for poisson/log", {
    beta.tgt <- c(`(Intercept)` = 0.5, x1 = 0.3, x2 = -0.2)
    o <- gen.syn.beta(beta.tgt, family = "poisson/log",
                      method = "project", n = 2000L, seed = 8L)
    b.refit <- coef(glm(o$formula, data = o$syn, family = poisson()))
    # Probabilistic rounding floor: paper documents a 2.3% (sd 1.1)
    # Monte Carlo mean.  Use 5% to keep the test stable across seeds.
    expect_equal(unname(b.refit), unname(beta.tgt), tolerance = 0.05)
})

test_that("gen.syn.beta() seed makes runs reproducible", {
    o1 <- gen.syn.beta(c(0, 1, -1), seed = 42L)
    o2 <- gen.syn.beta(c(0, 1, -1), seed = 42L)
    expect_equal(o1$syn, o2$syn)
})

test_that("gen.syn.beta() rejects non-numeric beta", {
    expect_error(gen.syn.beta("not numeric"), "numeric vector")
})

test_that("gen.syn.beta() rejects unknown family strings", {
    expect_error(gen.syn.beta(c(0, 1), family = "tweedie"),
                 "Unknown family")
})

test_that("gen.syn.beta() rejects X missing required columns", {
    expect_error(
        gen.syn.beta(c(`(Intercept)` = 0, x1 = 1, x2 = 1),
                     X = data.frame(x1 = rnorm(10))),
        "missing columns"
    )
})

test_that("gen.syn.beta() rejects marginals list missing entries", {
    expect_error(
        gen.syn.beta(c(`(Intercept)` = 0, x1 = 1, x2 = 1),
                     marginals = list(x1 = function(n) rnorm(n))),
        "missing entries"
    )
})

test_that("gen.syn.beta() rejects non-positive dispersion", {
    expect_error(gen.syn.beta(c(0, 1), dispersion = -1),
                 "positive scalar")
    expect_error(gen.syn.beta(c(0, 1), dispersion = 0),
                 "positive scalar")
})

test_that("gen.syn.beta() returns a synthexpln object printable by S3 method", {
    o <- gen.syn.beta(c(0, 1, -1), seed = 1L)
    out <- utils::capture.output(print(o))
    expect_true(any(grepl("from-beta", out)))
    expect_true(any(grepl("rows:", out)))
})
