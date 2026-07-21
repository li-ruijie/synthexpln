# Design-Gram projection: the copula base's exact-count permutation,
# the pinned-Gram OLS triple, the two-part corollary path, and the loud
# failure guards.

test_that("copula base with strata.exact pins every cell count", {
    set.seed(6489)
    d <- data.frame(
        y = rnorm(300, 10, 2),
        x = rnorm(300, 45, 9),
        g = factor(sample(c("a", "b", "c"), 300, TRUE)),
        h = factor(sample(c("u", "r"), 300, TRUE)))
    d.s <- synthexpln:::.copula.base(d, strata.exact = TRUE)
    expect_identical(table(d.s[["g"]], d.s[["h"]]),
                     table(d[["g"]], d[["h"]]))
    expect_true(any(as.character(d.s[["g"]]) != as.character(d[["g"]])))
})

test_that("design.exact OLS release carries the exact triple", {
    set.seed(6489)
    d <- data.frame(
        y = rnorm(300, 10, 2),
        x = rnorm(300, 45, 9),
        g = factor(sample(c("a", "b", "c"), 300, TRUE)))
    syn <- projection(d, y ~ x + g, seed = 11, design.exact = TRUE)
    expect_true(syn[["design.exact"]])
    m.o <- stats::model.matrix(~ x + g, d)
    m.s <- stats::model.matrix(~ x + g, syn[["syn"]])
    s.scl <- max(abs(crossprod(m.o)))
    expect_lt(max(abs(crossprod(m.s) - crossprod(m.o))), 1e-8 * s.scl)
    mo.o <- stats::lm(y ~ x + g, d)
    mo.s <- stats::lm(y ~ x + g, syn[["syn"]])
    expect_lt(max(abs(stats::coef(mo.s) - stats::coef(mo.o))), 1e-8)
    expect_lt(max(abs(stats::coef(summary(mo.s))[, 2] -
                      stats::coef(summary(mo.o))[, 2])), 1e-8)
    # submodel spot-check (the 2^p claim)
    mo.o2 <- stats::lm(y ~ x, d)
    mo.s2 <- stats::lm(y ~ x, syn[["syn"]])
    expect_lt(max(abs(stats::coef(mo.s2) - stats::coef(mo.o2))), 1e-8)
})

test_that("design.exact refuses the bootstrap base and unpinned draws", {
    d <- data.frame(y = rnorm(50), x = rnorm(50))
    expect_error(projection(d, y ~ x, base = "bootstrap",
                            design.exact = TRUE),
                 "design.exact")
    set.seed(2)
    d2 <- data.frame(
        y = rnorm(200), x = rnorm(200),
        g = factor(sample(c("a", "b"), 200, TRUE)))
    d.boot <- synthexpln:::.copula.base(d2, strata.exact = FALSE)
    expect_error(synthexpln:::.dg.project(d.boot, d2, "y", c("x", "g")),
                 "D'D not pinned")
})

test_that("two-part design.exact pins both block Grams and the triple", {
    set.seed(6489)
    v.x <- rnorm(400, 45, 9)
    d <- data.frame(
        y = ifelse(runif(400) < 0.35, 0,
                   rgamma(400, 2, 2 / exp(0.8 + 0.03 * v.x))),
        x = v.x,
        g = factor(sample(c("a", "b"), 400, TRUE)))
    syn <- projection.twopart(d, y ~ x + g, family = Gamma("log"),
                              seed = 21, design.exact = TRUE)
    expect_true(syn[["design.exact"]])
    d.s <- syn[["syn"]]
    m.o <- stats::model.matrix(~ x + g, d)
    m.s <- stats::model.matrix(~ x + g, d.s)
    s.scl <- max(abs(crossprod(m.o)))
    expect_lt(max(abs(crossprod(m.s) - crossprod(m.o))), 1e-8 * s.scl)
    m.o.p <- m.o[d[["y"]] > 0, , drop = FALSE]
    m.s.p <- m.s[d.s[["y"]] > 0, , drop = FALSE]
    expect_lt(max(abs(crossprod(m.s.p) - crossprod(m.o.p))), 1e-8 * s.scl)
    # positive-block refit: coefficients and standard errors
    mo.p0 <- stats::glm(y ~ x + g, d[d[["y"]] > 0, ],
                        family = stats::quasi(link = "log",
                                              variance = "mu^2"),
                        start = c(1, 0, 0))
    mo.ps <- stats::glm(y ~ x + g, d.s[d.s[["y"]] > 0, ],
                        family = stats::quasi(link = "log",
                                              variance = "mu^2"),
                        start = c(1, 0, 0),
                        control = stats::glm.control(epsilon = 1e-12,
                                                     maxit = 200))
    expect_lt(max(abs(stats::coef(mo.ps) - stats::coef(mo.p0))), 1e-4)
    expect_lt(max(abs(stats::coef(summary(mo.ps))[, 2] -
                      stats::coef(summary(mo.p0))[, 2])
                  / stats::coef(summary(mo.p0))[, 2]), 1e-3)
})

test_that("the bootstrap default is unchanged in shape and fields", {
    d <- data.frame(y = rnorm(80), x = rnorm(80))
    syn <- projection(d, y ~ x, seed = 5)
    expect_s3_class(syn, "synthexpln")
    expect_false(syn[["design.exact"]])
    expect_identical(syn[["variant"]], "projection")
    expect_identical(dim(syn[["syn"]]), dim(d))
})
