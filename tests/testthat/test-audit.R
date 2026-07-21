# The bound applied by dp.audit is (e^eps + delta) / (1 + e^eps + delta), due to
# Wasserman and Zhou (2010).  It falls as epsilon decreases and becomes
# vacuous as epsilon grows, so each test below selects epsilon accordingly.
#
# Both tests in this file were gated behind SYNTHEXPLN_FULL_TESTS=1 until
# 2026-08-06 and consequently not executed.  Runtime and stability both permit
# them to run by default.  The file completes in under one second, the V1
# verdict held on all 10 seeds tested, and the zero-noise AUC measured exactly 1
# on all 20 seeds tested.

.wz.bound <- function(eps, delta) (exp(eps) + delta) / (1 + exp(eps) + delta)

test_that("the audit bound is vacuous at large epsilon and binding at small", {
    # The zero-noise test below audits at eps = 1 rather than eps = 10 for this
    # reason.  At eps = 10 the bound exceeds every attainable AUC, so the audit
    # cannot reject a mechanism at that level, irrespective of the information
    # the mechanism discloses.
    expect_gt(.wz.bound(10, 1e-6), 0.9999)
    expect_equal(.wz.bound(1, 1e-6), 0.731059, tolerance = 1e-5)
    expect_lt(.wz.bound(0.1, 1e-6), 0.53)
})

test_that("dp.audit passes for V1 at eps=10", {
    set.seed(1)
    n <- 300
    d <- data.frame(y = rnorm(n), x = rnorm(n))
    gen <- function(data) suppressWarnings(
        gen.syn.dp.projected(data, y ~ x, epsilon = 10, delta = 1e-6))
    a <- dp.audit(gen, d, epsilon = 10, delta = 1e-6, n.trials = 50)
    expect_true(a[["pass"]])
})

test_that("dp.audit flags a zero-noise mechanism once the bound is binding", {
    set.seed(1)
    d <- data.frame(y = rnorm(50), x = rnorm(50))
    # The mechanism adds no DP noise.  The release is a deterministic function
    # of its input, so D and D' each yield a single point, the classifier
    # separates them perfectly, and the AUC is 1.
    bad.gen <- function(data) {
        structure(list(syn = data, beta.dp = stats::coef(stats::lm(y ~ x, data))),
                  class = "synthexpln")
    }
    a <- dp.audit(bad.gen, d, epsilon = 1, delta = 1e-6, n.trials = 50)

    expect_false(a[["pass"]])
    # Confirm that rejection follows from the AUC exceeding the bound itself,
    # and not solely from the Monte Carlo tolerance added to it.
    expect_gt(a[["auc"]], a[["theory.bound"]])
    expect_equal(a[["auc"]], 1)
})

test_that("the same zero-noise mechanism is NOT flagged at eps=10", {
    # The audit evaluates one specific (epsilon, delta) claim.  At eps = 10 the
    # bound is approximately 0.99995, and the rejection threshold is that value
    # plus the default tolerance of 0.02, namely 1.02.  An AUC is at most 1, so
    # the audit accepts every mechanism at this level, and a passing result
    # carries no information about the privacy of the mechanism under test.  The
    # earlier version of this file asserted rejection at eps = 10, a condition
    # the arithmetic excludes.
    set.seed(1)
    d <- data.frame(y = rnorm(50), x = rnorm(50))
    bad.gen <- function(data) {
        structure(list(syn = data, beta.dp = stats::coef(stats::lm(y ~ x, data))),
                  class = "synthexpln")
    }
    a <- dp.audit(bad.gen, d, epsilon = 10, delta = 1e-6, n.trials = 50)

    expect_true(a[["pass"]])
    expect_equal(a[["auc"]], 1)
})
