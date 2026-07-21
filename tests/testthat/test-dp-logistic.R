# V5: binary outcomes.  The sound path is subsample-and-aggregate over a public
# coefficient box; the DFBETA path is retired and warns.
#
# The box below is on the LOG-ODDS scale: effects of |beta| <= 3 are what these
# designs anticipate.  It is fixed before the data are drawn, does not read off d.
box.lo <- -3
box.hi <- 3

test_that("gen.syn.dp.logistic generates Bernoulli outcomes from DP beta", {
  set.seed(1)
  n <- 500
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
  p <- plogis(-0.5 + 0.8 * x1 - 0.3 * x2)
  y <- rbinom(n, 1, p)
  d <- data.frame(y = y, x1 = x1, x2 = factor(x2))

  syn <- gen.syn.dp.logistic(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6,
                              clip.lo = box.lo, clip.hi = box.hi, seed = 1)
  expect_true(all(syn$syn$y %in% c(0, 1)))
  expect_equal(syn$variant, "V5")
  expect_equal(syn$sens.method, "subagg")
  # sign agreement
  fit.syn <- glm(y ~ x1 + x2, family = binomial("logit"), syn$syn)
  expect_equal(unname(sign(coef(fit.syn))[2]), sign(0.8))
})

test_that("public X sends the whole mu budget to the coefficient", {
  # n = 1200 over m = 20 blocks leaves 60 rows per block for a 3-parameter
  # logistic fit.  Thinner blocks separate, and a separated block fit returns a
  # divergent coefficient that the public box then clips.  That is the mechanism
  # working as designed (the box exists for this case), but it makes the
  # fixture noisy, so give the blocks enough rows to fit.
  set.seed(12)
  n <- 1200
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$x1 - 0.4 * d$x2))
  eps <- 10; delta <- 1e-6
  syn <- gen.syn.dp.logistic(d, y ~ x1 + x2, epsilon = eps, delta = delta,
                              clip.lo = box.lo, clip.hi = box.hi, seed = 5)

  # public covariates cost nothing, so there is no X channel to fund
  expect_equal(syn$mu.channels[["coef"]], syn$mu.total)
  expect_equal(syn$mu.channels[["marg"]], 0)
  expect_equal(syn$mu.channels[["corr"]], 0)
  expect_equal(syn$mu.check, syn$mu.total)
  expect_equal(synthexpln:::.gdp.delta(eps, syn$mu.total), delta,
               tolerance = 1e-9)

  # the binomial has no free dispersion, so the aggregate stays p-dimensional
  # and there is no sigma^2 coordinate riding along
  expect_length(syn$beta.dp, 3L)
  expect_true(is.na(syn$subagg$sigma2.dp))
})

test_that("the public coefficient box is REQUIRED, not defaulted", {
  set.seed(13)
  d <- data.frame(y = rbinom(100, 1, 0.5), x1 = rnorm(100))
  expect_error(gen.syn.dp.logistic(d, y ~ x1, epsilon = 10),
               "clip.lo and clip.hi")
})

test_that("the retired DFBETA path is reachable but WARNS", {
  set.seed(14)
  d <- data.frame(y = rbinom(200, 1, 0.5), x1 = rnorm(200))
  expect_warning(
    syn <- gen.syn.dp.logistic(d, y ~ x1, epsilon = 10, sens.method = "local",
                               seed = 1),
    "NO formal \\(eps, delta\\)-DP guarantee")
  expect_equal(syn$sens.method, "local")
  expect_true(is.na(syn$mu.total))
})

test_that("sensitive-X path privatises X, keeps all coefficients, and agrees in sign", {
  set.seed(2)
  n  <- 3000
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rbinom(n, 1, 0.5)
  y  <- rbinom(n, 1, plogis(0.4 + 1.0 * x1 - 0.7 * x2 + 0.3 * x3))
  d  <- data.frame(y = y, x1 = x1, x2 = x2, x3 = factor(x3))
  x.spec <- list(
    x1 = list(type = "continuous",  breaks = seq(-4, 4, length.out = 11)),
    x2 = list(type = "continuous",  breaks = seq(-4, 4, length.out = 11)),
    x3 = list(type = "categorical", levels = c("0", "1")))
  syn <- gen.syn.dp.logistic(d, y ~ x1 + x2 + x3, epsilon = 10, delta = 1e-6,
                              x.public = FALSE, x.spec = x.spec,
                              clip.lo = box.lo, clip.hi = box.hi, seed = 5)
  expect_true(all(syn$syn$y %in% c(0, 1)))
  expect_equal(syn$variant, "V5")
  # X was privatised (not the original column)
  expect_false(isTRUE(all.equal(syn$syn$x1, d$x1)))
  # every coefficient survives the round trip (no silent factor drop)
  fit <- glm(y ~ x1 + x2 + x3, family = binomial("logit"), syn$syn)
  expect_setequal(names(coef(fit)), names(syn$beta.dp))
  # dominant slope keeps its sign
  expect_equal(unname(sign(coef(fit))["x1"]), 1)

  # default is independent marginals, and every channel is Gaussian, so the
  # whole release composes in one mu quadrature at the target delta
  expect_false(syn$use.copula)
  expect_equal(syn$mu.check, syn$mu.total)
  expect_equal(syn$delta.total, 1e-6)
  # with the copula off the correlation share is folded back into the marginals
  expect_equal(syn$mu.channels[["corr"]], 0)
  expect_gt(syn$mu.channels[["marg"]], 0)
})

test_that("the copula no longer charges a second delta", {
  set.seed(3)
  n  <- 2000
  x1 <- rnorm(n); x2 <- 0.6 * x1 + sqrt(1 - 0.36) * rnorm(n)
  y  <- rbinom(n, 1, plogis(0.3 + x1 - 0.5 * x2))
  d  <- data.frame(y = y, x1 = x1, x2 = x2)
  x.spec <- list(x1 = list(type = "continuous", breaks = seq(-4, 4, length.out = 11)),
                 x2 = list(type = "continuous", breaks = seq(-4, 4, length.out = 11)))
  syn <- gen.syn.dp.logistic(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6,
                              x.public = FALSE, x.spec = x.spec,
                              clip.lo = box.lo, clip.hi = box.hi,
                              use.copula = TRUE, seed = 1)
  expect_true(syn$use.copula)
  expect_true(all(syn$syn$y %in% c(0, 1)))

  # It used to report delta.total = 2 * delta: two Gaussian mechanisms composed
  # by BASIC composition, so delta accumulated per channel.  Under mu-GDP the
  # channels compose in quadrature and the whole release is calibrated to one
  # (eps, delta).  Delta does not accumulate.
  expect_equal(syn$delta.total, 1e-6)
  expect_equal(syn$mu.check, syn$mu.total)
})

test_that("sensitive-X without x.spec is an error", {
  d <- data.frame(y = rbinom(50, 1, 0.5), x1 = rnorm(50))
  expect_error(
    gen.syn.dp.logistic(d, y ~ x1, epsilon = 5, x.public = FALSE,
                        clip.lo = box.lo, clip.hi = box.hi, seed = 1),
    "x.spec")
})

test_that("score.method is validated and defaults to adaptive", {
  set.seed(4)
  n  <- 1500
  x1 <- rnorm(n); x2 <- 0.5 * x1 + sqrt(0.75) * rnorm(n)
  y  <- rbinom(n, 1, plogis(0.2 + x1 - 0.4 * x2))
  d  <- data.frame(y = y, x1 = x1, x2 = x2)
  x.spec <- list(x1 = list(type = "continuous", breaks = seq(-4, 4, length.out = 11)),
                 x2 = list(type = "continuous", breaks = seq(-4, 4, length.out = 11)))
  gen <- function(..., seed = 7)
    gen.syn.dp.logistic(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6,
                        x.public = FALSE, x.spec = x.spec,
                        clip.lo = box.lo, clip.hi = box.hi,
                        use.copula = TRUE, seed = seed, ...)
  expect_error(gen(score.method = "bogus"), "should be one of")

  # default equals explicit "adaptive"; "box" differs
  a  <- gen()
  b  <- gen(score.method = "adaptive")
  cc <- gen(score.method = "box")
  # x2's copula draw depends on the released correlation, so it reflects
  # score.method; x1 (the first copula column) is correlation-invariant by the
  # chol structure, so it is identical across methods and not a useful probe.
  expect_equal(a$syn$x2, b$syn$x2)
  expect_false(isTRUE(all.equal(a$syn$x2, cc$syn$x2)))

  # and score.method does not touch the coefficient channel
  expect_equal(a$beta.dp, cc$beta.dp)
})

test_that("score.method 'normal' recovers X-correlation while 'box' is noise-dominated", {
  set.seed(11)
  n  <- 2000
  x1 <- rnorm(n); x2 <- 0.6 * x1 + sqrt(1 - 0.36) * rnorm(n)
  y  <- rbinom(n, 1, plogis(0.3 + x1 - 0.5 * x2))
  d  <- data.frame(y = y, x1 = x1, x2 = x2)
  x.spec <- list(x1 = list(type = "continuous", breaks = seq(-4, 4, length.out = 11)),
                 x2 = list(type = "continuous", breaks = seq(-4, 4, length.out = 11)))
  # a descriptive-fidelity split: this release wants the X structure, so it
  # funds the correlation channel instead of starving it
  ms <- c(marg = 0.25, corr = 0.25, coef = 0.50)
  cor.of <- function(method) vapply(1:20, function(s)
    cor(gen.syn.dp.logistic(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6,
                            x.public = FALSE, x.spec = x.spec, mu.split = ms,
                            clip.lo = box.lo, clip.hi = box.hi,
                            use.copula = TRUE,
                            score.method = method, seed = s)$syn[c("x1", "x2")])[1, 2],
    numeric(1))
  v.norm <- cor.of("normal"); v.box <- cor.of("box")
  # The claim is RELATIVE, and it is about the signal-to-noise ratio: the box PIT
  # compresses concentrated data, so the second-moment signal falls far below the
  # DP noise, while the normal PIT spreads the scores and lifts it.  Assert that,
  # rather than an absolute SD threshold, which is a function of the budget split
  # and moves whenever mu.split does.
  expect_gt(stats::sd(v.box), 5 * stats::sd(v.norm))
  # normal recovers the true correlation, and closely
  expect_equal(mean(v.norm), cor(d$x1, d$x2), tolerance = 0.1)
  expect_lt(stats::sd(v.norm), 0.1)
})

test_that("score.method 'adaptive' recovers X-correlation on skewed margins where 'normal' over-clamps", {
  set.seed(21)
  n  <- 2000
  z1 <- rnorm(n); z2 <- 0.6 * z1 + sqrt(1 - 0.36) * rnorm(n)
  x1 <- exp(0.8 * z1); x2 <- exp(0.8 * z2)      # right-skewed log-normal margins
  y  <- rbinom(n, 1, plogis(0.3 + as.numeric(scale(x1)) - 0.5 * as.numeric(scale(x2))))
  d  <- data.frame(y = y, x1 = x1, x2 = x2)
  x.spec <- list(x1 = list(type = "continuous", breaks = seq(0, 10, length.out = 11)),
                 x2 = list(type = "continuous", breaks = seq(0, 10, length.out = 11)))
  ms <- c(marg = 0.25, corr = 0.25, coef = 0.50)
  clamp.of <- function(method) mean(vapply(1:20, function(s)
    abs(cor(gen.syn.dp.logistic(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6,
            x.public = FALSE, x.spec = x.spec, mu.split = ms,
            clip.lo = box.lo, clip.hi = box.hi, use.copula = TRUE,
            score.method = method, seed = s)$syn[c("x1", "x2")])[1, 2]) > 0.9,
    logical(1)))
  # the location-scale "normal" PIT mis-scales the skew and clamps to |r| ~ 1
  # often; "adaptive" maps through each margin's own DP histogram and rarely does
  expect_gt(clamp.of("normal"), 0.5)
  expect_lt(clamp.of("adaptive"), 0.2)
})
