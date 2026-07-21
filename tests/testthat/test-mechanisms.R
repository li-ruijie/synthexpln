test_that(".dp.gm.multiplier: classical branch for eps < 1, analytic for eps >= 1", {
  delta <- 1e-6
  # classical branch (Dwork-Roth 2014 Thm A.1), exact formula
  expect_equal(synthexpln:::.dp.gm.multiplier(0.5, delta),
               sqrt(2 * log(1.25 / delta)) / 0.5, tolerance = 1e-12)
  # analytic branch: the returned multiplier certifies exactly delta at eps
  # (Balle-Wang 2018; .gdp.delta is the mechanism's (eps, delta) profile)
  for (eps in c(1, 4, 10)) {
    m <- synthexpln:::.dp.gm.multiplier(eps, delta)
    expect_equal(synthexpln:::.gdp.delta(eps, 1 / m), delta, tolerance = 1e-9)
  }
  # the classical multiplier under-delivers delta for eps >= 1 (the bug the
  # analytic branch fixes): at eps = 10 it certifies only ~1.9e-6, not 1e-6
  m.classical <- sqrt(2 * log(1.25 / delta)) / 10
  expect_gt(synthexpln:::.gdp.delta(10, 1 / m.classical), delta)
})

test_that(".dp.gaussian.noise has correct standard deviation", {
  set.seed(1)
  noise <- replicate(5000, synthexpln:::.dp.gaussian.noise(
    s.p = 1, eps = 1, delta = 1e-6, s.Delta = 1))
  # eps = 1 sits on the analytic branch of the calibrated multiplier
  expected.sd <- synthexpln:::.dp.gm.multiplier(1, 1e-6)
  expect_equal(sd(noise), expected.sd, tolerance = 0.02 * expected.sd)
})

test_that(".dp.gaussian applies noise to every element", {
  set.seed(1)
  v <- rep(0, 10)
  out <- synthexpln:::.dp.gaussian(v, eps = 1, delta = 1e-6, s.Delta = 1)
  expect_equal(length(out), 10)
  expect_gt(sum(abs(out)), 0)  # noise added
})

test_that(".dp.hist.continuous adds Laplace noise per bin", {
  set.seed(1)
  x <- rnorm(200)
  breaks <- seq(-3, 3, length.out = 11)
  counts <- synthexpln:::.dp.hist.continuous(x, breaks = breaks, eps = 1)
  expect_equal(length(counts), 10)
  expect_true(all(counts >= 0))  # Laplace with rounding
})

test_that(".dp.cat perturbs categorical counts via Laplace mechanism", {
  set.seed(1)
  x <- sample(letters[1:3], 100, replace = TRUE)
  out <- synthexpln:::.dp.cat(x, levels = letters[1:3], eps = 1)
  expect_equal(sort(unique(out)), letters[1:3])
  expect_equal(length(out), 100)
})

test_that(".qdphist inverts the DP histogram CDF (monotone, in range, matches sampler)", {
  set.seed(1)
  breaks <- seq(0, 10, length.out = 6)          # 5 bins on [0, 10]
  counts <- c(10, 20, 40, 20, 10)               # unimodal, sum 100
  # in range and monotone non-decreasing
  u <- seq(0.001, 0.999, length.out = 50)
  q <- synthexpln:::.qdphist(u, counts, breaks)
  expect_true(all(q >= 0 & q <= 10))
  expect_true(all(diff(q) >= -1e-9))
  # median near the middle bin [4, 6]
  expect_gt(synthexpln:::.qdphist(0.5, counts, breaks), 3)
  expect_lt(synthexpln:::.qdphist(0.5, counts, breaks), 7)
  # inverse-CDF of uniforms reproduces the forward sampler's distribution
  set.seed(2); a <- synthexpln:::.qdphist(runif(20000), counts, breaks)
  set.seed(2); b <- synthexpln:::.sample.from.dp.hist(counts, breaks, 20000)
  expect_lt(abs(mean(a) - mean(b)), 0.1)
  expect_lt(abs(sd(a) - sd(b)), 0.1)
})

test_that(".qdphist handles guard paths: zero counts, interior empty bin, boundary u", {
  breaks <- seq(0, 10, length.out = 6)          # 5 bins on [0, 10]
  # (a) all-zero counts fall back to a uniform histogram, output in range
  q0 <- synthexpln:::.qdphist(seq(0.01, 0.99, length.out = 20), rep(0, 5), breaks)
  expect_true(all(q0 >= 0 & q0 <= 10))
  expect_true(all(diff(q0) >= -1e-9))
  # (b) an interior zero-probability bin is skipped, not rested inside
  br3 <- seq(0, 9, length.out = 4)              # 3 bins: [0,3] [3,6] [6,9]
  qz  <- synthexpln:::.qdphist(seq(0.001, 0.999, length.out = 200), c(10, 0, 10), br3)
  expect_false(any(qz > 3 + 1e-6 & qz < 6 - 1e-6))
  expect_true(all(diff(qz) >= -1e-9))
  # (c) literal u = 0 and u = 1 are clamped (not errors) and stay in range
  qb <- synthexpln:::.qdphist(c(0, 1), c(10, 20, 40, 20, 10), breaks)
  expect_equal(length(qb), 2L)
  expect_true(all(qb >= 0 & qb <= 10))
})

test_that(".dp.copula.scores are bounded, monotone, and row-local", {
  x   <- c(-5, 0, 3, 5, 20)                     # includes out-of-box values
  lo  <- 0; hi <- 10; tau <- 1 / (2 * 100)
  w   <- synthexpln:::.dp.copula.scores(x, lo, hi, tau)
  M   <- qnorm(1 - tau)
  expect_true(all(abs(w) <= M + 1e-9))          # bounded by M
  expect_true(all(diff(w) >= -1e-9))            # monotone in x
  expect_equal(w[1], qnorm(tau), tolerance = 1e-9)      # x = -5 clips to lo
  expect_equal(w[5], qnorm(1 - tau), tolerance = 1e-9)  # x = 20 clips to hi
  x2 <- x; x2[3] <- 9
  w2 <- synthexpln:::.dp.copula.scores(x2, lo, hi, tau)
  expect_equal(w[-3], w2[-3])                    # row-local: only w[3] changes
})

test_that(".pdphist is the forward CDF: in [0,1], monotone, inverts .qdphist", {
  breaks <- seq(0, 10, length.out = 6)          # 5 bins on [0, 10]
  counts <- c(10, 20, 40, 20, 10)               # unimodal, sum 100
  x <- seq(-2, 12, length.out = 60)             # includes out-of-range values
  Fx <- synthexpln:::.pdphist(x, counts, breaks)
  expect_true(all(Fx >= 0 & Fx <= 1))           # a CDF
  expect_true(all(diff(Fx) >= -1e-9))           # non-decreasing
  expect_equal(synthexpln:::.pdphist(-100, counts, breaks), 0)    # flat below range
  expect_equal(synthexpln:::.pdphist(100, counts, breaks), 1)     # flat above range
  # round trip F(Q(u)) = u on the DP histogram
  u <- seq(0.01, 0.99, length.out = 50)
  q <- synthexpln:::.qdphist(u, counts, breaks)
  expect_lt(max(abs(synthexpln:::.pdphist(q, counts, breaks) - u)), 1e-6)
})

test_that(".dp.copula.scores 'adaptive' is bounded and needs counts + breaks", {
  set.seed(3)
  breaks <- seq(-4, 4, length.out = 11)
  x      <- rnorm(2000)
  counts <- synthexpln:::.dp.hist.continuous(x, breaks, eps = 5)
  tau <- 0.005; M <- qnorm(1 - tau)
  w   <- synthexpln:::.dp.copula.scores(x, -4, 4, tau, method = "adaptive",
                                     counts = counts, breaks = breaks)
  expect_true(all(abs(w) <= M + 1e-9))          # bounded by the score clamp
  # on N(0,1) data the adaptive PIT tracks the standardised x closely
  expect_gt(cor(w, pmin(pmax(x, -M), M)), 0.99)
  # counts/breaks are required for the adaptive transform
  expect_error(synthexpln:::.dp.copula.scores(x, -4, 4, tau, method = "adaptive"),
               "requires counts and breaks")
})

test_that(".dp.copula.corr returns a valid correlation matrix and tracks eps", {
  set.seed(11)
  n <- 4000; tau <- 1 / (2 * n); M <- qnorm(1 - tau)
  # correlated latent: two columns with true rho ~ 0.6
  z1 <- rnorm(n); z2 <- 0.6 * z1 + sqrt(1 - 0.36) * rnorm(n)
  W  <- cbind(z1, z2)
  R  <- synthexpln:::.dp.copula.corr(W, eps = 50, delta = 1e-6, s.M = M, seed = 1)
  # valid correlation matrix: symmetric, unit diagonal, PSD
  expect_equal(dim(R), c(2L, 2L))
  expect_equal(diag(R), c(1, 1), tolerance = 1e-8)
  expect_equal(R, t(R), tolerance = 1e-10)
  expect_gte(min(eigen(R, only.values = TRUE)[["values"]]), -1e-8)
  # at large eps the off-diagonal recovers the latent correlation
  expect_equal(R[1, 2], 0.6, tolerance = 0.05)
  # smaller eps -> more noise -> off-diagonal farther from truth on average
  err.hi <- mean(replicate(40, abs(synthexpln:::.dp.copula.corr(
    W, 50, 1e-6, M)[1, 2] - 0.6)))
  err.lo <- mean(replicate(40, abs(synthexpln:::.dp.copula.corr(
    W, 1, 1e-6, M)[1, 2] - 0.6)))
  expect_lt(err.hi, err.lo)
})

test_that(".sample.from.dp.copula reproduces marginals and the target correlation", {
  set.seed(21)
  n <- 8000
  br1 <- seq(0, 10, length.out = 11); br2 <- seq(0, 20, length.out = 11)
  m1 <- synthexpln:::.dp.hist.continuous(runif(n, 0, 10),  br1, 5)
  m2 <- synthexpln:::.dp.hist.continuous(runif(n, 0, 20),  br2, 5)
  l.marg <- list(a = list(counts = m1, breaks = br1),
                 b = list(counts = m2, breaks = br2))
  R <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
  d <- synthexpln:::.sample.from.dp.copula(R, l.marg, n, seed = 7)
  expect_equal(names(d), c("a", "b"))
  expect_equal(nrow(d), n)
  # marginals stay within their public boxes
  expect_true(all(d$a >= 0 & d$a <= 10))
  expect_true(all(d$b >= 0 & d$b <= 20))
  # copula induces the requested rank dependence
  expect_gt(cor(d$a, d$b, method = "spearman"), 0.5)
  # independence target (R = I) gives near-zero dependence
  d0 <- synthexpln:::.sample.from.dp.copula(diag(2), l.marg, n, seed = 7)
  expect_lt(abs(cor(d0$a, d0$b, method = "spearman")), 0.1)
})

test_that(".dp.copula.corr noise sd matches the calibrated sensitivity constant", {
  set.seed(41)
  n <- 3000; p <- 2L; tau <- 1 / (2 * n); M <- qnorm(1 - tau)
  eps <- 4; delta <- 1e-6
  # calibrated per-entry Gaussian sd = Delta_2 * .dp.gm.multiplier(eps, delta).
  # Read Delta from the mechanism's own definition rather than re-deriving it: a
  # test that recomputes the constant it guards tests only its own arithmetic.
  sig.expect <- synthexpln:::.dp.copula.corr.delta(n, p, M) *
    synthexpln:::.dp.gm.multiplier(eps, delta)
  # Fix W so the only variability across reps is the DP noise, not the
  # empirical correlation's 1/sqrt(n) sampling variability; the released
  # off-diagonal sd then equals the calibrated per-entry sigma.
  W   <- matrix(qnorm(runif(n * p, tau, 1 - tau)), n, p)
  r12 <- replicate(400, synthexpln:::.dp.copula.corr(W, eps, delta, M)[1, 2])
  # a wrong constant (missing factor 2, or p vs p^2) shifts this sd out of band
  # (tolerance is relative: 0.15 = 15%, a factor-2 error is 50%+ and fails)
  expect_equal(sd(r12), sig.expect, tolerance = 0.15)
})

test_that(".dp.hist.continuous uses the replace-one 2/eps Laplace scale (review C1)", {
  set.seed(101)
  # Large, known per-bin counts so the >= 0 clip effectively does not fire,
  # leaving (released - true) equal to the raw Laplace noise.
  breaks      <- 0:10                                   # 10 unit-width bins
  x           <- rep(seq(0.5, 9.5, by = 1), each = 100) # exactly 100 per bin
  true.counts <- rep(100, 10)
  eps         <- 1
  # The replace-one (Hamming-1) count sensitivity is 2, so the calibrated
  # Laplace scale is 2/eps and the per-count noise sd is (2/eps) * sqrt(2).
  # .dp.cat shares this scale.  The old add/remove bug used 1/eps (sd sqrt(2)),
  # 50% low, which this 5% band rejects.
  expect.sd <- (2 / eps) * sqrt(2)
  noise <- replicate(
    2000, synthexpln:::.dp.hist.continuous(x, breaks, eps) - true.counts)
  expect_equal(sd(as.numeric(noise)), expect.sd, tolerance = 0.05 * expect.sd)
})

test_that("the count mechanisms take EXACTLY ONE of eps and mu", {
  # A channel calibrated in eps (Laplace) and a channel calibrated in mu
  # (Gaussian) cannot be composed with one another: the Laplace release is pure
  # eps-DP and does not enter a mu quadrature at all.  Supplying both, or
  # neither, is therefore a category error and not a defaultable convenience.
  x <- runif(200, 0, 10)
  breaks <- 0:10
  expect_error(synthexpln:::.dp.hist.continuous(x, breaks),
               "exactly one of eps")
  expect_error(synthexpln:::.dp.hist.continuous(x, breaks, eps = 1, mu = 1),
               "exactly one of eps")
  expect_error(synthexpln:::.dp.cat(letters[1:3], letters[1:3]),
               "exactly one of eps")
  expect_error(synthexpln:::.dp.cat(letters[1:3], letters[1:3], eps = 1, mu = 1),
               "exactly one of eps")
})

test_that(".dp.hist.continuous's mu path uses the sqrt(2) L2 sensitivity", {
  # A replace-one row swap decreases one bin count and increases another, so
  # the count vector moves by (-1, +1) and its L2 sensitivity is sqrt(2), not 2
  # (that is the L1 sensitivity, which is what the Laplace path is calibrated
  # to).  The Gaussian noise sd is therefore sqrt(2) / mu.
  set.seed(202)
  breaks      <- 0:10
  x           <- rep(seq(0.5, 9.5, by = 1), each = 100)
  true.counts <- rep(100, 10)
  mu          <- 0.5

  expect.sd <- sqrt(2) / mu
  noise <- replicate(
    2000, synthexpln:::.dp.hist.continuous(x, breaks, mu = mu) - true.counts)
  expect_equal(sd(as.numeric(noise)), expect.sd, tolerance = 0.05 * expect.sd)

  # and it is Gaussian, not Laplace: excess kurtosis ~ 0, where a Laplace would
  # sit at 3
  z <- as.numeric(noise)
  expect_equal(mean((z - mean(z))^4) / stats::sd(z)^4, 3, tolerance = 0.25)
})

test_that(".dp.cat.probs is the RELEASE, and .dp.cat only post-processes it", {
  # The mechanism must be observable.  .dp.cat used to return only a sample, so
  # the released probability vector could not be looked at: at any usable mu the
  # multinomial sampling noise swamps the DP noise, and a test written against
  # the sample measures the sampler, not the channel.
  lv <- letters[1:4]
  x  <- rep(lv, each = 250)

  # the release is a probability vector over the declared levels
  p <- synthexpln:::.dp.cat.probs(x, lv, mu = 1)
  expect_length(p, 4L)
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0))

  # the sample is drawn from it and nothing else, so it is post-processing
  set.seed(301); a <- synthexpln:::.dp.cat(x, lv, mu = 1)
  set.seed(301)
  b <- sample(lv, length(x), replace = TRUE,
              prob = synthexpln:::.dp.cat.probs(x, lv, mu = 1))
  expect_equal(a, b)
})

test_that(".dp.cat's mu path uses the sqrt(2) L2 sensitivity", {
  # Same argument as the histogram: a replace-one swap moves one unit between two
  # levels, so the count vector moves by (-1, +1) and its L2 sensitivity is
  # sqrt(2).  Read the noise straight off the released counts.
  set.seed(203)
  lv          <- letters[1:4]
  x           <- rep(lv, each = 2500)      # large counts: the >= 0 clip does not fire
  true.counts <- rep(2500, 4)
  mu          <- 0.5
  N           <- length(x)

  # invert the renormalisation: probs * N recovers the noised counts, since the
  # noise is small relative to N and the total is preserved to first order
  noise <- replicate(2000, {
    p <- synthexpln:::.dp.cat.probs(x, lv, mu = mu)
    p * sum(true.counts) - true.counts
  })
  expect_equal(sd(as.numeric(noise)), sqrt(2) / mu,
               tolerance = 0.1 * sqrt(2) / mu)

  # and it is Gaussian, not Laplace: excess kurtosis ~ 0, where Laplace sits at 3
  z <- as.numeric(noise)
  expect_equal(mean((z - mean(z))^4) / stats::sd(z)^4, 3, tolerance = 0.3)
})

test_that("EPSILON BITES on the categorical channel", {
  # A channel flat in mu is a channel the noise did not reach.  Measured on the
  # release, not on the sample, so the sampler's own noise cannot mask it.
  lv <- letters[1:4]
  x  <- rep(lv, each = 250)
  spread <- vapply(c(0.05, 0.2, 1, 5, 25), function(m) {
    set.seed(204)
    mean(apply(replicate(300, synthexpln:::.dp.cat.probs(x, lv, mu = m)),
               1, stats::sd))
  }, numeric(1))
  expect_true(all(diff(spread) < 0))
})

test_that(".dp.copula.corr's mu path calibrates sigma = Delta / mu", {
  set.seed(205)
  n <- 500; p <- 3
  m.W <- matrix(rnorm(n * p), n, p)
  s.M <- qnorm(1 - 0.005)
  mu  <- 2

  # Delta_2 = M^2 p / n (the L2 sensitivity of the released vech(S)), so
  # sigma = Delta_2 / mu.  A large mu drives the noise to zero and the release must
  # converge on the (rescaled) second-moment matrix.
  R.exact <- synthexpln:::.dp.copula.corr(m.W, NULL, 1e-6, s.M, seed = 1,
                                          mu = 1e8)
  m.S <- crossprod(m.W) / n
  v.d <- 1 / sqrt(diag(m.S))
  R.true <- v.d * m.S * rep(v.d, each = p)
  expect_equal(R.exact, (R.true + t(R.true)) / 2, tolerance = 1e-5)

  # and a small mu is visibly noisier than a large one
  R.noisy <- synthexpln:::.dp.copula.corr(m.W, NULL, 1e-6, s.M, seed = 1,
                                          mu = 0.05)
  expect_gt(max(abs(R.noisy - R.true)), max(abs(R.exact - R.true)))

  # the result is always a valid correlation matrix: PSD, unit diagonal
  expect_equal(diag(R.noisy), rep(1, p), tolerance = 1e-8)
  expect_gte(min(eigen(R.noisy, symmetric = TRUE, only.values = TRUE)$values),
             -1e-8)
})

test_that("the mu-GDP layer composes in quadrature and inverts exactly", {
  expect_equal(synthexpln:::.gdp.compose(c(3, 4)), 5)
  expect_equal(synthexpln:::.gdp.compose(c(1, 1, 1)), sqrt(3))

  for (eps in c(0.25, 1, 5, 20)) {
    delta <- 1e-6
    mu <- synthexpln:::.gdp.mu.from.eps.delta(eps, delta)
    # the dual is exact: mu is the multiplier reaching (eps, delta)
    expect_equal(synthexpln:::.gdp.delta(eps, mu), delta, tolerance = 1e-9)
    # and .gdp.sigma is Delta / mu
    expect_equal(synthexpln:::.gdp.sigma(2.5, eps, delta), 2.5 / mu)
  }

  # .gdp.to.delta is .gdp.delta with the arguments swapped
  expect_equal(synthexpln:::.gdp.to.delta(1.5, 3),
               synthexpln:::.gdp.delta(3, 1.5))
})

test_that(".dp.hist.continuous releases the HISTOGRAM, so its sensitivity is 2", {
  # The test above checks the noise scale against an assumed sensitivity of 2.
  # It cannot check that 2 is the right sensitivity, and its own fixture hides the
  # question: x puts 100 rows in every bin, so no bin is ever empty.
  #
  # That gap shipped a defect.  Until 2026-07-13 the counts came from table() on
  # the bin indices, and table() drops the bins no row occupies, so the occupied
  # bins were left-packed and the zeros padded onto the end.  An empty interior
  # bin slid every count above it down one slot.  The released vector was
  # therefore not the histogram, and 2 is a HISTOGRAM's sensitivity.
  #
  # Assert both halves: the vector is the histogram, and one row moving cannot
  # shift it by more than the mechanism is calibrated for.
  breaks <- 0:5                                     # 5 unit-width bins
  # Bin 3 ([2,3)) is EMPTY, and bins 4 and 5 are not.  The old code returned
  # c(2, 1, 1, 1, 0) here, relocating both upper bins one slot down.
  x <- c(0.5, 0.5, 1.5, 3.5, 4.5)
  expect_equal(synthexpln:::.dp.hist.continuous(x, breaks, eps = 1e12),
               c(2, 1, 0, 1, 1), tolerance = 1e-6)

  # Sensitivity, measured on the neighbouring pair an adversary would pick: one
  # row moves INTO the empty interior bin.  eps very large, so the release is the
  # count vector itself and the distance below is the count function's own.
  l1 <- function(a, b) sum(abs(a - b))
  k  <- 100L
  a  <- c(rep(0.5, k), rep(2.5, k))                 # bin 2 empty
  b  <- c(rep(0.5, k), rep(2.5, k - 1L), 1.5)       # one row moved into it
  s.sens <- l1(synthexpln:::.dp.hist.continuous(a, breaks, eps = 1e12),
               synthexpln:::.dp.hist.continuous(b, breaks, eps = 1e12))
  # The Laplace scale is 2/eps, i.e. an L1 sensitivity of 2.  The old code moved
  # by 2 * (k - 1) = 198 here, under-noising this channel by 99x, and by O(n) in
  # general.
  expect_lte(s.sens, 2 + 1e-6)
})

test_that(".dp.copula.corr releases vech(S), so its sensitivity is M^2 p / n", {
  # Same defect class as the test above, same channel family, found a day later.
  #
  # The mechanism draws p(p+1)/2 iid normals, lays them on the upper triangle
  # INCLUDING the diagonal, and mirrors.  The mirroring is post-processing, so the
  # released vector is vech(S), not the p x p matrix S.  A Gaussian mechanism is
  # mu-GDP at mu = Delta / sigma with Delta the L2 sensitivity OF the released
  # vector, so Delta must be the sensitivity of vech(S) and of nothing else.
  #
  # Until 2026-07-14 Delta was 2 M^2 p / n, the triangle-inequality bound on the
  # FROBENIUS norm of the full matrix.  Frobenius counts each off-diagonal entry
  # twice where vech counts it once, so it is sound (it can only over-noise) and it
  # is not the sensitivity: this channel drew 2x the noise its own guarantee needed.
  #
  # Assert both halves, because soundness alone is what let it hide:
  #   (a) no neighbouring pair moves the release by more than Delta   [guarantee]
  #   (b) the worst pair ATTAINS Delta                                [sharpness]
  # A bound nothing attains was derived from some other object.
  n <- 400L; p <- 2L
  tau <- 0.005
  M   <- qnorm(1 - tau)
  lo  <- 0; hi <- 5

  set.seed(53L)
  m.bg <- matrix(runif((n - 1L) * p, 0.2, 4.8), n - 1L, p)
  scores.with <- function(v.last) vapply(seq_len(p), function(j)
    synthexpln:::.dp.copula.scores(c(m.bg[, j], v.last[j]), lo, hi, tau, method = "box"),
    numeric(n))
  vech.l2 <- function(m.D) sqrt(sum(m.D[upper.tri(m.D, diag = TRUE)]^2))
  move.of <- function(v.a, v.b) {
    m.A <- scores.with(v.a); m.B <- scores.with(v.b)
    expect_lte(max(abs(m.A)), M + 1e-6)          # the clamp the bound rests on
    expect_lte(max(abs(m.B)), M + 1e-6)
    vech.l2(crossprod(m.A) / n - crossprod(m.B) / n)
  }

  s.Delta <- synthexpln:::.dp.copula.corr.delta(n, p, M)

  # (a) soundness, over adversarial corners and a random sweep of the box
  v.corner <- c(move.of(c(+1e6, +1e6), c(+1e6, -1e6)),
                move.of(c(2.5, 2.5),   c(-1e6, +1e6)),
                move.of(c(+1e6, +1e6), c(-1e6, +1e6)))
  set.seed(431L)
  v.rand <- replicate(200L, move.of(runif(p, -8, 13), runif(p, -8, 13)))
  expect_lte(max(c(v.corner, v.rand)), s.Delta + 1e-9)

  # (b) sharpness.  Both canary rows clamp to the corner (a record outside the
  # public box lands exactly on +-M) and the opposite sign in the second column
  # makes them orthogonal, which zeroes the cross term and the diagonal of the
  # move.  That is the maximiser, and it attains Delta exactly.  On the old
  # constant this is 0.5 * Delta and the test fails, reporting the 2x over-noising.
  expect_equal(move.of(c(+1e6, +1e6), c(+1e6, -1e6)), s.Delta, tolerance = 1e-9)
})
