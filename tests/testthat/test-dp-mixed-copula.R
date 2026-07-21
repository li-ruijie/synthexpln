# the mixed-margin copula.  Categorical covariates join the same Gaussian
# copula as the continuous block, each encoded by K-1 bounded indicator scores
# computed from its released DP proportions, and are reconstructed through a
# calibrated probit/logit back-map.  Descriptive fidelity only: the refit
# recovers beta.dp exactly regardless.  Mirrors lib/35-dp-full.r and its
# guard verify-dp-soundness.R Task 7.

# A public coefficient box and a public sigma^2 box for the fixtures below,
# facts about each test's DGP fixed before any data are drawn (coefficients
# 1, 0.05, -0.03, 2, 4 and sigma^2 = 1).
box.lo <- -10
box.hi <- 10

# Correlation ratio eta, the categorical-continuous association metric the
# paper's 58e MC reports (the absolute point-biserial for a binary covariate).
eta.ratio <- function(v.cat, v.num) {
  s.gm <- mean(v.num)
  l.g  <- Filter(function(g) length(g) > 0L, split(v.num, v.cat))
  ss.b <- sum(vapply(l.g, function(g) length(g) * (mean(g) - s.gm)^2,
                     numeric(1)))
  sqrt(ss.b / sum((v.num - s.gm)^2))
}

# Fixture: g's three levels sit at distinct x1 locations, so eta(g, x1) is
# strong (0.90 at this seed) and an independent-margin release destroys it.
mk.mixed.d <- function(n = 600L, seed = 42L) {
  set.seed(seed)
  g  <- factor(sample(c("a", "b", "c"), n, replace = TRUE,
                      prob = c(0.5, 0.3, 0.2)), levels = c("a", "b", "c"))
  x1 <- pmin(pmax(rnorm(n, c(a = 25, b = 50, c = 75)[as.character(g)], 10),
                  0), 100)
  x2 <- pmin(pmax(0.5 * x1 + rnorm(n, 25, 10), 0), 100)
  y  <- 1 + 0.05 * x1 - 0.03 * x2 + 2 * (g == "b") + 4 * (g == "c") + rnorm(n)
  data.frame(y = y, x1 = x1, x2 = x2, g = g)
}
mixed.spec <- list(x1 = list(type = "continuous", bounds = c(0, 100)),
                   x2 = list(type = "continuous", bounds = c(0, 100)),
                   g  = list(type = "categorical", levels = c("a", "b", "c")))
# the descriptive-fidelity split the paper's 58e MC uses
ms.descr <- c(marg = 0.25, corr = 0.25, coef = 0.50)

gen.mixed <- function(d, use.copula, back.map = "probit", seed = 7)
  gen.syn.dp.full(d, y ~ x1 + x2 + g, x.spec = mixed.spec, epsilon = 20,
                  delta = 1e-6, mu.split = ms.descr,
                  clip.lo = box.lo, clip.hi = box.hi,
                  disp.lo = 0, disp.hi = 9,
                  use.copula = use.copula, back.map = back.map, seed = seed)

test_that("categoricals join the copula: extended release, calibrated marginal", {
  d   <- mk.mixed.d()
  syn <- gen.mixed(d, use.copula = TRUE)

  # the released correlation covers p_total = 2 continuous + (K-1) = 2
  # indicator columns, not just the continuous block
  expect_true(syn$use.copula)
  expect_equal(dim(syn$dp.corr), c(4L, 4L))

  # the categorical column is returned as a factor on the declared levels
  expect_s3_class(syn$syn$g, "factor")
  expect_equal(levels(syn$syn$g), c("a", "b", "c"))
  expect_equal(nrow(syn$syn), 600L)

  # the calibrated back-map reproduces the released marginal: at eps = 20 the
  # DP noise on the proportions is small, so the synthetic proportions must
  # sit near the true ones (0.527, 0.273, 0.200 at this fixture seed)
  v.props <- as.numeric(prop.table(table(syn$syn$g)))
  expect_lt(max(abs(v.props - as.numeric(prop.table(table(d$g))))), 0.08)

  # and the back-map argument is validated
  expect_error(gen.mixed(d, use.copula = TRUE, back.map = "nonsense"))
})

test_that("mixed copula is descriptive-only: the refit recovers beta.dp exactly", {
  d <- mk.mixed.d()
  Map(function(bm) {
    syn <- gen.mixed(d, use.copula = TRUE, back.map = bm)
    b   <- coef(lm(y ~ x1 + x2 + g, syn$syn))
    v.common <- intersect(names(b), names(syn$beta.dp))
    # same 1e-5 slack as the continuous copula test: the only slack is the
    # 1e-10 ridge on crossprod(X_syn)
    expect_lt(max(abs(b[v.common] - syn$beta.dp[v.common])), 1e-5)
  }, c("probit", "logit"))
})

test_that("probit recovers the categorical-continuous association, independent destroys it", {
  d <- mk.mixed.d()
  eta.true <- eta.ratio(d$g, d$x1)                 # 0.90 at this seed
  expect_gt(eta.true, 0.8)

  syn.cop <- gen.mixed(d, use.copula = TRUE)
  syn.ind <- gen.mixed(d, use.copula = FALSE)
  eta.cop <- eta.ratio(syn.cop$syn$g, syn.cop$syn$x1)
  eta.ind <- eta.ratio(syn.ind$syn$g, syn.ind$syn$x1)

  # measured at this seed with the A6 de-attenuation: 0.67 against 0.06 (and
  # 0.61/0.10, 0.61/0.06 at the neighbouring seeds), so these margins are wide
  expect_gt(eta.cop, 0.5)
  expect_lt(eta.ind, 0.15)
  expect_gt(eta.cop, 2 * eta.ind)
})

test_that("logit back-map ATTENUATES the association (A2)", {
  d <- mk.mixed.d()
  syn.p <- gen.mixed(d, use.copula = TRUE, back.map = "probit")
  syn.l <- gen.mixed(d, use.copula = TRUE, back.map = "logit")
  eta.p <- eta.ratio(syn.p$syn$g, syn.p$syn$x1)
  eta.l <- eta.ratio(syn.l$syn$g, syn.l$syn$x1)
  # measured 0.67 against 0.27 at this seed (0.61/0.32, 0.61/0.30 nearby): the
  # softmax's independent choice noise costs about half the recovered
  # association, so "probit" is the default and "logit" is caveated
  expect_lt(eta.l, eta.p - 0.1)
})

test_that("de-attenuation (A6) rescales only the cat-cont cross entries, from released probs", {
  # a 3-column released correlation: 2 continuous + 1 indicator (K = 2 at
  # released p = c(0.5, 0.5)), with a known cross entry
  m.R <- diag(3)
  m.R[1, 2] <- m.R[2, 1] <- 0.40                 # cont-cont: must NOT move
  m.R[1, 3] <- m.R[3, 1] <- 0.30                 # cat-cont: divided by c_k
  m.R[2, 3] <- m.R[3, 2] <- 0.10                 # cat-cont: divided by c_k
  l.marg <- list(x1 = list(counts = 1, breaks = 0:1),
                 x2 = list(counts = 1, breaks = 0:1),
                 g  = list(levels = c("a", "b"), probs = c(0.5, 0.5)))
  c.k <- dnorm(qnorm(0.5)) / sqrt(0.5 * 0.5)     # 0.7979, the p = 0.5 factor
  m.out <- synthexpln:::.mixed.deattenuate(m.R, l.marg, c("x1", "x2"), "g")
  # the corrected matrix is already PSD here, so the PSD projection is the
  # identity and the entries are exact
  expect_equal(m.out[1, 2], 0.40, tolerance = 1e-8)
  expect_equal(m.out[1, 3], 0.30 / c.k, tolerance = 1e-8)
  expect_equal(m.out[2, 3], 0.10 / c.k, tolerance = 1e-8)
  expect_equal(diag(m.out), rep(1, 3), tolerance = 1e-8)
  # extreme fitted entries clamp at 0.99 before the PSD projection
  m.R2 <- m.R
  m.R2[1, 3] <- m.R2[3, 1] <- 0.95               # 0.95 / 0.7979 > 1
  m.out2 <- synthexpln:::.mixed.deattenuate(m.R2, l.marg, c("x1", "x2"), "g")
  expect_lte(max(abs(m.out2)), 1 + 1e-8)
})

test_that("V5 logistic inherits the mixed-margin copula", {
  set.seed(5)
  n  <- 600L
  g  <- factor(sample(c("a", "b", "c"), n, replace = TRUE,
                      prob = c(0.5, 0.3, 0.2)), levels = c("a", "b", "c"))
  x1 <- pmin(pmax(rnorm(n, c(a = 25, b = 50, c = 75)[as.character(g)], 10),
                  0), 100)
  x2 <- pmin(pmax(0.5 * x1 + rnorm(n, 25, 10), 0), 100)
  y  <- rbinom(n, 1, plogis(-2 + 0.03 * x1 + 0.6 * (g == "c")))
  d  <- data.frame(y = y, x1 = x1, x2 = x2, g = g)
  spec <- list(x1 = list(type = "continuous", breaks = seq(0, 100, length.out = 11)),
               x2 = list(type = "continuous", breaks = seq(0, 100, length.out = 11)),
               g  = list(type = "categorical", levels = c("a", "b", "c")))
  # suppressWarnings: the m = 20 blocks fit 5 logistic coefficients on 30 rows
  # each, so glm.fit occasionally reports separated blocks ("fitted
  # probabilities numerically 0 or 1").  That is a property of the small
  # per-block fits, orthogonal to the X-copula behaviour under test, and the
  # clipped aggregate absorbs it by construction.
  gen5 <- function(cop) suppressWarnings(gen.syn.dp.logistic(
      d, y ~ x1 + x2 + g, epsilon = 20, delta = 1e-6,
      x.public = FALSE, x.spec = spec, clip.lo = -5, clip.hi = 5,
      mu.split = c(marg = 0.25, corr = 0.25, coef = 0.50),
      use.copula = cop, seed = 7))
  syn.c <- gen5(TRUE)
  syn.i <- gen5(FALSE)
  # the release is a joint one: binary response, factor levels intact
  expect_true(all(syn.c$syn$y %in% 0:1))
  expect_equal(levels(syn.c$syn$g), c("a", "b", "c"))
  expect_equal(dim(syn.c$dp.corr), c(4L, 4L))    # p_total = 2 + (K - 1)
  # the mixed copula recovers the g-x1 association, independent margins lose it
  eta.c <- eta.ratio(syn.c$syn$g, syn.c$syn$x1)
  eta.i <- eta.ratio(syn.i$syn$g, syn.i$syn$x1)
  expect_gt(eta.c, 2 * eta.i)
})

test_that("extended vech(S) sensitivity is SOUND, LOAD-BEARING, and NOT attained (A3)", {
  # The testthat analogue of verify-dp-soundness.R Task 7.  Hands the SHIPPED
  # score construction adversarial neighbouring pairs that drive the
  # continuous rows onto the clamp AND flip a categorical level, and asserts:
  #   (a) sound: every pair moves vech(S) by <= Delta(p_total).
  #   (b) TEETH: some level-flip pair moves it PAST the continuous-only
  #       Delta(p_c), so the indicator columns genuinely enter the release
  #       and p_total, not p_c, is the width of the sensitivity.
  #   (c) A3: the max move stays strictly below Delta(p_total): a categorical
  #       column takes one of K fixed winsorised values strictly inside the
  #       clamp, so the extended bound is a sound upper bound, over-noising
  #       the categorical channel in the safe direction.
  s.tau <- 0.005
  s.M   <- qnorm(1 - s.tau)
  n     <- 400L
  v.lev <- c("a", "b", "c")                        # "c" rare, for the teeth
  brk   <- seq(0, 5, length.out = 6)
  s.tol <- 1e-9

  set.seed(57)
  v.u <- runif(n - 1L, 0.2, 4.8)
  v.v <- runif(n - 1L, 0.2, 4.8)
  v.g <- c(rep("a", round((n - 1L) * 0.60)), rep("b", round((n - 1L) * 0.38)))
  v.g <- c(v.g, rep("c", (n - 1L) - length(v.g))) # ~2% prevalence of "c"
  cnt.u <- synthexpln:::.dp.hist.continuous(v.u, brk, mu = 5)
  cnt.v <- synthexpln:::.dp.hist.continuous(v.v, brk, mu = 5)
  prb.g <- synthexpln:::.dp.cat.probs(v.g, v.lev, mu = 5)

  scores <- function(s.u, s.v, s.g) cbind(
    synthexpln:::.dp.copula.scores(c(v.u, s.u), brk[1], brk[6], s.tau,
                                   "adaptive", counts = cnt.u, breaks = brk),
    synthexpln:::.dp.copula.scores(c(v.v, s.v), brk[1], brk[6], s.tau,
                                   "adaptive", counts = cnt.v, breaks = brk),
    synthexpln:::.cat.indicator.scores(factor(c(v.g, s.g), levels = v.lev),
                                       v.lev, prb.g, s.M))
  vech.l2 <- function(m) sqrt(sum(m[upper.tri(m, diag = TRUE)]^2))
  move <- function(a, b) {
    m.A <- scores(a$u, a$v, a$g)
    m.B <- scores(b$u, b$v, b$g)
    expect_lte(max(abs(m.A)), s.M + s.tol)         # every column stays clamped
    expect_lte(max(abs(m.B)), s.M + s.tol)
    vech.l2(crossprod(m.A) / n - crossprod(m.B) / n)
  }

  p.tot  <- ncol(scores(1, 1, "a"))                # 2 continuous + 2 indicators
  expect_equal(p.tot, 4L)
  D.tot  <- synthexpln:::.dp.copula.corr.delta(n, p.tot, s.M)
  D.cont <- synthexpln:::.dp.copula.corr.delta(n, 2L, s.M)

  l.corner <- list(
    list(a = list(u = +1e6, v = +1e6, g = "a"),
         b = list(u = +1e6, v = -1e6, g = "c")),   # clamp corners + flip to rare c
    list(a = list(u = +1e6, v = +1e6, g = "c"),
         b = list(u = -1e6, v = +1e6, g = "b")),   # clamp corners + flip c to b
    list(a = list(u = 2.5, v = 2.5, g = "a"),
         b = list(u = 2.5, v = 2.5, g = "c")))     # flip only, continuous centred
  v.corner <- vapply(l.corner, function(p) move(p$a, p$b), numeric(1))
  set.seed(571)
  v.rand <- replicate(100, {
    g2 <- sample(v.lev, 2L, replace = TRUE)
    move(list(u = runif(1, -8, 13), v = runif(1, -8, 13), g = g2[1]),
         list(u = runif(1, -8, 13), v = runif(1, -8, 13), g = g2[2]))
  })
  s.max <- max(c(v.corner, v.rand))

  expect_lte(s.max, D.tot + s.tol)                 # (a) THE GUARD
  expect_gt(max(v.corner), D.cont)                 # (b) TEETH: p_total load-bearing
  expect_lt(s.max, D.tot - s.tol)                  # (c) A3: strict upper bound
})

test_that("mixed calibration matches the source tree's fixed-point sequence", {
  # .mixed.cat.calibrate is the source loop of lib/35-dp-full.r expressed as
  # tail recursion: same init, same damped update, same convergence test, same
  # cap.  Fed the same pool and target it must return the same shifts as an
  # inline transliteration of the source loop.
  set.seed(11)
  m.Z <- matrix(rnorm(20000 * 2), 20000, 2)
  v.p <- c(0.5, 0.3, 0.2)
  src.loop <- function(m.Z.cal, v.probs, back.map) {
    v.shift <- log(pmax(v.probs[-1L], 1e-6) / max(v.probs[1L], 1e-6))
    for (it in seq_len(100L)) {
      v.pr <- synthexpln:::.mixed.cat.props(m.Z.cal, v.shift, back.map)
      if (max(abs(v.pr - v.probs)) < 1e-4) break
      v.shift <- v.shift + 0.7 * log(pmax(v.probs[-1L], 1e-6) /
                                     pmax(v.pr[-1L], 1e-6))
    }
    v.shift
  }
  Map(function(bm) expect_equal(
        synthexpln:::.mixed.cat.calibrate(m.Z, v.p, bm),
        src.loop(m.Z, v.p, bm)),
      c("probit", "logit"))
})
