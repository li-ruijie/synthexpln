# test-dp-suffstat.R -- V6, the whitened sufficient-statistic release.
# Mirrors verify-gdp-mechanisms.R Task 11 of the replication repository: the mu
# accounting, the exact refit recovery, determinism, the dispersion modes,
# the PD repair, the lambda rules, the refusal of missing public constants,
# and the inference pair (the 2026-07-16 coverage finding).

mk.d <- function(n = 400L) {
    set.seed(42)
    d <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
    d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + stats::rnorm(n)
    d
}
W3 <- local({
    m <- diag(3)
    dimnames(m) <- list(c("(Intercept)", "x1", "x2"),
                        c("(Intercept)", "x1", "x2"))
    m
})
ARGS <- list(B.x = sqrt(1 + 2 * 6^2), y.centre = 0, y.scale = 1,
             y.clip = c(-8, 8))

gen <- function(d, ...) do.call(gen.syn.dp.suffstats,
    c(list(d = d, formula = y ~ x1 + x2, W = W3), ARGS, list(...)))

test_that("mu accounting composes in quadrature and hits the dual", {
    r <- gen(mk.d(), epsilon = 10, seed = 7)
    expect_equal(sqrt(sum(r$mu.channels^2)), r$mu.total, tolerance = 1e-9)
    expect_equal(synthexpln:::.gdp.to.delta(r$mu.total, 10), 1e-6,
                 tolerance = 1e-9)
})

test_that("the public-X refit recovers beta.dp exactly", {
    r <- gen(mk.d(), epsilon = 10, x.public = TRUE, seed = 7)
    expect_equal(unname(coef(stats::lm(y ~ x1 + x2, r$syn))),
                 unname(r$beta.dp), tolerance = 1e-8)
})

test_that("the release is deterministic under the seed", {
    r1 <- gen(mk.d(), epsilon = 5, seed = 11)
    r2 <- gen(mk.d(), epsilon = 5, seed = 11)
    expect_identical(r1$beta.dp, r2$beta.dp)
    expect_identical(r1$S.dp, r2$S.dp)
})

test_that("the public-dispersion mode drops the y'y channel", {
    disp <- dp.dispersion.public(1.0, source = "test constant")
    r <- gen(mk.d(), epsilon = 10, dispersion = disp, seed = 7)
    expect_identical(unname(r$mu.channels[["yy"]]), 0)
    expect_identical(r$sigma2.dp, 1.0)
    expect_true(is.na(r$q.dp))
    expect_equal(sqrt(sum(r$mu.channels^2)), r$mu.total, tolerance = 1e-9)
})

test_that("the PD repair keeps the solve finite at a small budget", {
    r <- gen(mk.d(), epsilon = 0.1, seed = 3)
    expect_true(all(is.finite(r$beta.dp)))
    expect_gte(r$lambda.eff, r$lambda)
})

test_that("the lambda rules behave", {
    r.e <- gen(mk.d(), epsilon = 10, seed = 7)      # released-eigenmin default
    expect_gt(r.e$mu.channels[["eig"]], 0)
    r.n <- gen(mk.d(), epsilon = 10, lambda = 40, seed = 7)
    expect_identical(unname(r.n$mu.channels[["eig"]]), 0)
    expect_identical(r.n$lambda, 40)
    expect_error(gen(mk.d(), epsilon = 10, lambda = -1), "non-negative")
})

test_that("missing public constants are refused", {
    expect_error(gen.syn.dp.suffstats(mk.d(), y ~ x1 + x2, epsilon = 1),
                 "PUBLIC")
})

test_that("the inference pair satisfies the de-shrink identity", {
    r <- gen(mk.d(), epsilon = 10, seed = 7)
    expect_equal(r$beta.infer - r$beta.dp, r$bias.fold, tolerance = 1e-9)
    # At eps = 10 the eigenvalue floor cannot bind, so the plain solve must
    # reproduce beta.infer exactly (on the whitened scale, backed out raw).
    v.w0 <- solve(r$S.dp, r$b.dp)
    v.raw <- as.numeric(W3 %*% v.w0) * ARGS$y.scale
    v.raw[[1L]] <- v.raw[[1L]] + ARGS$y.centre
    expect_equal(unname(r$beta.infer), v.raw, tolerance = 1e-8)
    expect_lt(max(abs(r$Cov.beta.infer - t(r$Cov.beta.infer))), 1e-8)
    expect_true(all(diag(r$Cov.beta.infer) > 0))
})

test_that("a row outside the public ball is norm-clipped, not trusted", {
    d <- mk.d()
    d[1L, c("x1", "x2")] <- c(50, -50)     # far outside |x_j| <= 6
    r <- gen(d, epsilon = 10, seed = 7)
    expect_true(all(is.finite(r$beta.dp)))
    # The clip bounds what enters S: every whitened row norm is at most B.x,
    # so trace(S) is at most n B.x^2, up to the released diagonal noise (three
    # N(0, sigma.vech) entries, allowed six joint standard deviations here).
    expect_lte(sum(diag(r$S.dp)),
               nrow(d) * ARGS$B.x^2 + 6 * sqrt(3) * r$sigma.channels[["vech"]])
})
