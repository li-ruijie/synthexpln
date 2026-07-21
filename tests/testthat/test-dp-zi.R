# gen.syn.dp.projected.zi is the sound subsample-and-aggregate zero-inflated
# generator: a mu_Z-GDP logistic zero-model and a mu_P-GDP positive component,
# composed in quadrature.  Both coefficient channels clip each per-block fit into
# a public box, so the boxes are required arguments and no bound on |y| is needed.

make.zi <- function(n, seed) {
    set.seed(seed)
    x1 <- rnorm(n); x2 <- rnorm(n)
    z  <- rbinom(n, 1, plogis(0.3 + 0.9 * x1 - 0.5 * x2))
    mu <- exp(0.6 + 0.4 * x1 - 0.2 * x2)
    y  <- ifelse(z == 1, rgamma(n, shape = 2, scale = mu / 2), 0)
    data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("gen.syn.dp.projected.zi releases a sound mu-GDP ZI dataset", {
    d   <- make.zi(1500, 6489)
    syn <- gen.syn.dp.projected.zi(
        d, y ~ x1 + x2, epsilon = 10, m = 20,
        clip.lo = -8, clip.hi = 8, clip.z.lo = -8, clip.z.hi = 8,
        dispersion = dp.dispersion.public(1, source = "design residual variance"),
        family = quasi(link = "log", variance = "mu^2"), seed = 1)

    expect_s3_class(syn, "synthexpln")
    expect_identical(syn[["variant"]], "ZI-subagg")
    expect_identical(syn[["scope"]], "analysis-level")

    # the two channels compose in quadrature, and the total is the exact dual of
    # the requested (eps, delta)
    expect_equal(syn[["mu.total"]], sqrt(syn[["mu.Z"]]^2 + syn[["mu.P"]]^2))
    expect_equal(unname(syn[["mu.total"]]),
                 synthexpln:::.gdp.mu.from.eps.delta(10, 1e-6),
                 tolerance = 1e-6)

    # the release carries both a zero mass and a positive part
    expect_true(any(syn[["syn"]][["y"]] == 0))
    expect_true(any(syn[["syn"]][["y"]] > 0))
    expect_length(syn[["Z.syn"]], nrow(d))

    # zero-model decision agreement: the synthetic zero-model refit recovers the
    # private gamma_DP slope signs
    d.zr       <- syn[["syn"]]
    d.zr[[".Z"]] <- as.integer(syn[["syn"]][["y"]] > 0)
    g.syn      <- coef(glm(.Z ~ x1 + x2, d.zr, family = binomial(link = "logit")))
    v.slope    <- c("x1", "x2")
    expect_identical(unname(sign(g.syn[v.slope])),
                     unname(sign(syn[["gamma.dp"]][v.slope])))

    # public boxes are required (the box width sets the noise scale)
    expect_error(
        gen.syn.dp.projected.zi(d, y ~ x1 + x2, epsilon = 10,
            dispersion = dp.dispersion.public(1, source = "x"),
            family = quasi(link = "log", variance = "mu^2")),
        "clip")
    # a public dispersion is required for a free-dispersion family
    expect_error(
        gen.syn.dp.projected.zi(d, y ~ x1 + x2, epsilon = 10,
            clip.lo = -8, clip.hi = 8, clip.z.lo = -8, clip.z.hi = 8,
            family = quasi(link = "log", variance = "mu^2")),
        "dispersion")
})

test_that("the positive block stays mu_P-GDP under a positivity toggle", {
    # Adversarial corner-clipping positives: slopes far beyond the +-8 box so
    # each thin block clips to a corner by its composition.  A one-row
    # replacement that toggles positivity changes the sensitive count |P|.  If
    # the partition is taken over |P| it reshuffles and the released aggregate
    # moves by more than the calibrated (hi - lo)/m sensitivity; the fixed-N
    # partition keeps it within sqrt(p).  This is verify-dp-soundness.R Task 6.
    m <- 20L; sqrt.p <- sqrt(3); box.lo <- -8; box.hi <- 8

    d.pos <- local({
        n.pos <- 110L
        set.seed(3001L)
        typ <- rep(c(1, -1), length.out = n.pos)
        x1  <- runif(n.pos, -0.6, 0.6); x2 <- runif(n.pos, -0.6, 0.6)
        eta <- pmin(pmax(1 + typ * 12 * x1 + typ * 10 * x2, -5), 6)
        data.frame(y = exp(eta), x1 = x1, x2 = x2)
    })

    # scaled aggregate Q = a / s from the shipped internal, under partition `blk`
    # (NULL reproduces the pre-fix partition over the positive subset).  The
    # per-block quasi-gamma fits diverge on this deliberately extreme data (that
    # drives them onto the clip box), so their convergence warnings are expected
    # and suppressed; the clipping is the mechanism under test.
    Q.of <- function(dd, blk) {
        l <- suppressWarnings(synthexpln:::.subsample.aggregate.beta.loglink(
            dd, y ~ x1 + x2, quasi(link = "log", variance = "mu^2"),
            "quasi-gamma", "log", eps = NA_real_, delta = 1e-6, m = m,
            clip.lo = box.lo, clip.hi = box.hi, seed = 6489L, mu = 1,
            block = blk))
        a <- l[["beta.agg"]]
        a / l[["s.j"]][names(a)]
    }
    l2 <- function(u, v) sqrt(sum((u - v)^2))

    # partition the fixed full N (positives + 40 structural zeros) and restrict
    # to the positive rows, which occupy the first nrow(d.pos) slots here.
    n.all <- nrow(d.pos) + 40L
    set.seed(909L); blk.full <- sample(rep_len(seq_len(m), n.all))
    blk.pos <- blk.full[seq_len(nrow(d.pos))]

    # toggle the three highest-leverage positive rows in turn, take the worst
    v.k <- order(abs(d.pos[["x1"]]) + abs(d.pos[["x2"]]), decreasing = TRUE)[1:3]

    move.sound <- max(vapply(v.k, function(k) {
        l2(Q.of(d.pos, blk.pos),
           Q.of(d.pos[-k, , drop = FALSE], blk.pos[-k]))
    }, numeric(1)))
    move.prefix <- max(vapply(v.k, function(k) {
        l2(Q.of(d.pos, NULL), Q.of(d.pos[-k, , drop = FALSE], NULL))
    }, numeric(1)))

    # the fix keeps every toggle within sqrt(p); the pre-fix |P| partition does
    # not, so the guard has teeth
    expect_lte(move.sound, sqrt.p + 1e-6)
    expect_gt(move.prefix, sqrt.p)
})

test_that("gen.syn.dp.projected.zi scales coefficient noise with epsilon", {
    d    <- make.zi(1200, 11)
    args <- list(formula = y ~ x1 + x2, m = 20,
                 clip.lo = -8, clip.hi = 8, clip.z.lo = -8, clip.z.hi = 8,
                 dispersion = dp.dispersion.public(1, source = "design variance"),
                 family = quasi(link = "log", variance = "mu^2"))
    s.lo <- do.call(gen.syn.dp.projected.zi,
                    c(list(d = d, epsilon = 1,  seed = 2), args))
    s.hi <- do.call(gen.syn.dp.projected.zi,
                    c(list(d = d, epsilon = 20, seed = 2), args))
    # both the positive and zero-model noise scales shrink as epsilon grows
    expect_gt(mean(s.lo[["sigma.beta"]]),  mean(s.hi[["sigma.beta"]]))
    expect_gt(mean(s.lo[["sigma.gamma"]]), mean(s.hi[["sigma.gamma"]]))
})
