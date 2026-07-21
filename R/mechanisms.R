# mechanisms.R -- Internal DP mechanism helpers
# Ported from lib/28-dp-projection.r (lines 184-198, plus the
# dp.gm.multiplier calibration added there by d26d6ca) and
# lib/35-dp-full.r (lines 33-78).
#
# API NOTE: Source returns structured lists; these helpers return bare vectors
# for cleaner downstream composition (see function-level docs).

# ---------------------------------------------------------------------------
# Private RNG helper
# ---------------------------------------------------------------------------

#' Sample from the Laplace distribution
#'
#' @param n Number of draws.
#' @param location Location parameter (default 0).
#' @param scale Scale parameter (default 1).
#' @return Numeric vector of length \code{n}.
#' @noRd
#' @importFrom stats runif rnorm
.rlaplace <- function(n, location = 0, scale = 1) {
    v.u <- runif(n) - 0.5
    location - scale * sign(v.u) * log(1 - 2 * abs(v.u))
}

# ---------------------------------------------------------------------------
# Gaussian mechanism helpers
# ---------------------------------------------------------------------------

#' (eps, delta) profile of a Gaussian mechanism
#'
#' The delta certified at privacy level \code{eps} by a Gaussian mechanism
#' with noise multiplier \code{mu = s.Delta / sigma} (Balle-Wang 2018,
#' equivalently the mu-GDP profile of Dong-Roth-Su 2022).  Monotone
#' increasing in \code{mu}.  The e^eps term is evaluated in log-space so
#' the profile stays numerically stable at large \code{eps}.
#'
#' @param eps Privacy budget (epsilon > 0).
#' @param mu Noise multiplier s.Delta / sigma (> 0).
#' @return The certified delta.
#' @noRd
.gdp.delta <- function(eps, mu) {
    stats::pnorm(mu / 2 - eps / mu) -
        exp(eps + stats::pnorm(-mu / 2 - eps / mu, log.p = TRUE))
}

#' Gaussian-mechanism noise multiplier (per unit sensitivity)
#'
#' Returns sigma / Delta for an (eps, delta)-DP Gaussian mechanism at
#' L2-sensitivity Delta.  For eps < 1 the classical Dwork-Roth (2014,
#' Thm A.1) bound sqrt(2 log(1.25/delta)) / eps applies.  That bound is
#' valid only for eps < 1: for eps >= 1 it can under-noise (at eps = 10,
#' delta = 1e-6 it certifies only delta = 1.9e-6), so for eps >= 1 the
#' multiplier is calibrated by the exact analytic Gaussian mechanism of
#' Balle and Wang (2018), valid for all eps > 0.  Every Gaussian release
#' in the package draws its noise scale from this helper.
#'
#' Ported verbatim from \code{dp.gm.multiplier} in
#' \file{28-dp-projection.r} (code-tree commit d26d6ca), so package and
#' code-tree releases share one calibration.
#'
#' @param eps Privacy budget (epsilon > 0).
#' @param delta Privacy parameter (delta in (0, 1)).
#' @return Scalar noise multiplier sigma / Delta.
#' @noRd
.dp.gm.multiplier <- function(eps, delta) {
    if (eps < 1)                                  # classical, valid eps < 1
        return(sqrt(2 * log(1.25 / delta)) / eps)
    mu <- stats::uniroot(function(m) .gdp.delta(eps, m) - delta,
                         interval = c(1e-8, 1e3), tol = 1e-12)[["root"]]
    1 / mu
}

# ---------------------------------------------------------------------------
# mu-GDP composition and calibration
# ---------------------------------------------------------------------------
# Ported from lib/28-dp-projection.r lines 347-361.
#
# A release built from several Gaussian channels has to state one privacy
# level, and the only way to get there is to compose the channels.  Gaussian
# DP (Dong, Roth and Su 2022) composes in quadrature, which the classical
# (eps, delta) pair does not: two (eps, delta) mechanisms give (2 eps, 2
# delta) by basic composition, a bound that degrades linearly and that no
# amount of care recovers from, so every channel is calibrated in mu, the
# total is sqrt(sum mu_i^2), and the (eps, delta) label is read off the exact
# dual at the very end.
#
# The package was missing this.  V2 and V3 spent an EPS budget across a
# Gaussian beta channel and a Laplace sigma^2 channel, and a Laplace channel
# does not enter a mu quadrature at all, so the two could not be composed into
# the mu-GDP statement the replication repository and the paper both make.

#' Compose mu-GDP channels in quadrature (internal)
#'
#' k mechanisms that are mu_1, ..., mu_k-GDP compose to a
#' sqrt(sum mu_i^2)-GDP mechanism.
#'
#' @param v.mu Numeric vector of per-channel mu.
#' @return Scalar total mu.
#' @noRd
.gdp.compose <- function(v.mu) sqrt(sum(v.mu^2))

#' (eps, delta) profile at multiplier mu, argument order (mu, eps) (internal)
#' @noRd
.gdp.to.delta <- function(mu, eps) .gdp.delta(eps, mu)

#' Largest mu admissible at a target (eps, delta) (internal)
#'
#' Inverts \code{.gdp.delta} exactly, for all eps > 0.  Unlike
#' \code{.dp.gm.multiplier} this does NOT fall back to the classical
#' Dwork-Roth bound below eps = 1: the exact dual is valid everywhere and is
#' smaller, so it is the right calibration whenever the release is stated in
#' mu.  The two therefore agree exactly for eps >= 1 and the exact form gives
#' less noise for eps < 1.
#'
#' @param eps Privacy budget (epsilon > 0).
#' @param delta Privacy parameter (delta in (0, 1)).
#' @return Scalar mu.
#' @noRd
.gdp.mu.from.eps.delta <- function(eps, delta) {
    stopifnot(delta > 0, delta < 1, eps > 0)
    stats::uniroot(function(m) .gdp.delta(eps, m) - delta,
                   interval = c(1e-8, 1e3), tol = 1e-12)[["root"]]
}

#' Gaussian sd reaching (eps, delta) at L2 sensitivity s.Delta (internal)
#'
#' sigma = s.Delta / mu*, with mu* the exact dual of (eps, delta).
#'
#' @param s.Delta L2 sensitivity of the query.
#' @param eps,delta Target privacy parameters.
#' @return Scalar noise sd.
#' @noRd
.gdp.sigma <- function(s.Delta, eps, delta)
    s.Delta / .gdp.mu.from.eps.delta(eps, delta)

#' Add Gaussian noise to a numeric vector (calibrated Gaussian mechanism)
#'
#' Implements the Gaussian mechanism at sigma = s.Delta *
#' \code{.dp.gm.multiplier(eps, delta)}: the classical Dwork-Roth (2014,
#' Thm A.1) bound for eps < 1 and the exact analytic Balle-Wang (2018)
#' calibration for eps >= 1, where the classical bound under-delivers
#' delta.
#'
#' Differs from source \code{dp.gaussian} in \file{28-dp-projection.r} (which
#' returns a list containing \code{v.x.dp}, \code{sigma}, and \code{eta});
#' this version returns the noised vector directly.
#'
#' @param v.x Numeric vector to privatise.
#' @param eps Privacy budget (epsilon > 0).
#' @param delta Privacy parameter (delta in (0, 1)).
#' @param s.Delta L2 sensitivity bound on the query.
#' @param seed Optional integer RNG seed; \code{NA} to skip (default).
#' @return Numeric vector of the same length as \code{v.x} with i.i.d.
#'   N(0, sigma^2) noise added.
#' @noRd
.dp.gaussian <- function(v.x, eps, delta, s.Delta, seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    s.sigma <- s.Delta * .dp.gm.multiplier(eps, delta)
    v.x + rnorm(length(v.x), 0, s.sigma)
}

#' Draw Gaussian noise vector (calibrated Gaussian mechanism)
#'
#' Like \code{.dp.gaussian} but returns only the noise term rather than
#' the privatised query value.  Useful for pre-generating noise pools.
#'
#' Differs from source \code{dp.gaussian.noise} in \file{28-dp-projection.r}
#' (which returns a list \code{list(v.eta, sigma)}); this version returns the
#' noise vector directly so that \code{replicate()} composition is natural.
#'
#' @param s.p Length of the noise vector to generate.
#' @param eps Privacy budget (epsilon > 0).
#' @param delta Privacy parameter (delta in (0, 1)).
#' @param s.Delta L2 sensitivity bound on the query.
#' @param seed Optional integer RNG seed; \code{NA} to skip (default).
#' @return Numeric vector of length \code{s.p} drawn from
#'   N(0, sigma^2) where sigma = s.Delta * \code{.dp.gm.multiplier(eps,
#'   delta)}.
#' @noRd
.dp.gaussian.noise <- function(s.p, eps, delta, s.Delta, seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    s.sigma <- s.Delta * .dp.gm.multiplier(eps, delta)
    rnorm(s.p, 0, s.sigma)
}

# ---------------------------------------------------------------------------
# Continuous histogram mechanism
# ---------------------------------------------------------------------------

#' DP histogram for a continuous variable (Laplace mechanism)
#'
#' Bins \code{x} using \code{breaks}, adds i.i.d. Laplace(0, 2/eps) noise to
#' each bin count, and clips to non-negative values.  Scale 2/eps is the
#' replace-one (Hamming-1) count sensitivity of 2 (review C1).
#'
#' The count vector is built with \code{factor(levels = seq_len(B))} so that
#' bins no row occupies are retained as zeros.  This is load-bearing for the
#' privacy claim, not a presentational detail.  Until 2026-07-13 the counts came
#' from \code{table()} on the raw bin indices, which drops empty bins, and the
#' remainder was zero-padded on the right.  An empty interior bin therefore slid
#' every count above it down one slot, and the released vector was not the
#' histogram.  The sensitivity of 2 is the histogram's (a replace-one swap
#' decreases one bin's count and increases another's); the left-packed vector
#' has L1 sensitivity O(n), since moving one row into an empty interior bin
#' shifts every count above it.  The mechanism was therefore under-noised by a
#' factor that grows with the sample size.  Covered by the sensitivity test in
#' \file{tests/testthat/test-mechanisms.R}.
#'
#' Differs from source \code{dp.hist.continuous} in \file{35-dp-full.r} (which
#' returns \code{list(breaks, probs, B)}); this version returns the noised
#' count vector directly.  Renormalisation and bin-sampling are the caller's
#' responsibility.
#'
#' Two calibrations, exactly one of \code{eps} or \code{mu}:
#' \describe{
#'   \item{\code{mu}}{Gaussian mechanism, mu-GDP.  A replace-one row swap
#'     decreases one count and increases another, so the count vector moves by
#'     (-1, +1) and its L2 sensitivity is sqrt(2); the noise sd is
#'     sqrt(2) / mu.  This lets the marginal channel join the mu-GDP quadrature
#'     with the correlation and coefficient releases.}
#'   \item{\code{eps}}{Legacy Laplace(2 / eps).  The L1 sensitivity is 2 under
#'     replace-one (review C1).  Pure eps-DP, so it does NOT enter a mu-GDP
#'     quadrature, and it is therefore confined to the retired
#'     \code{sens.method = "local"} path.}
#' }
#'
#' @param x Numeric vector of observations.
#' @param breaks Numeric vector of public bin boundaries (length B + 1).
#' @param eps Privacy budget (epsilon > 0) for this marginal, Laplace path.
#' @param mu Per-channel mu for this marginal, Gaussian mu-GDP path.  Supply
#'   exactly one of \code{eps} and \code{mu}.
#' @return Non-negative numeric vector of length \code{length(breaks) - 1}
#'   containing DP-noised bin counts.
#' @noRd
.dp.hist.continuous <- function(x, breaks, eps = NULL, mu = NULL) {
    if (is.null(eps) == is.null(mu))
        stop("supply exactly one of eps (Laplace) or mu (Gaussian, mu-GDP)")
    B <- length(breaks) - 1L
    v.counts <- as.numeric(table(factor(
        findInterval(x, breaks[-c(1L, B + 1L)]) + 1L, levels = seq_len(B))))
    v.noise <- if (is.null(mu)) .rlaplace(B, 0, 2 / eps)
               else stats::rnorm(B, 0, sqrt(2) / mu)
    pmax(0, v.counts + v.noise)
}

#' Quantile function (inverse CDF) of a DP histogram
#'
#' Inverts the piecewise-linear CDF implied by DP-noised bin counts, mapping
#' uniforms in [0, 1] to values in [breaks[1], breaks[B + 1]].  The forward
#' companion is \code{.sample.from.dp.hist}: \code{.qdphist(runif(n), ...)}
#' and \code{.sample.from.dp.hist(..., n)} share the same law.  Used by the
#' Gaussian-copula sampler to give each covariate its DP marginal.
#'
#' @param u Numeric vector of probabilities in [0, 1].
#' @param counts Non-negative DP-noised bin counts (length B).
#' @param breaks Numeric bin boundaries (length B + 1).
#' @return Numeric vector of the same length as \code{u}.
#' @noRd
.qdphist <- function(u, counts, breaks) {
    B      <- length(counts)
    probs  <- if (sum(counts) > 0) counts / sum(counts) else rep(1 / B, B)
    cum    <- cumsum(probs)                             # length B, cum[B] = 1
    u      <- pmin(pmax(u, 1e-12), 1 - 1e-12)
    b      <- pmin(findInterval(u, cum) + 1L, B)        # first bin with cum >= u
    lo.cum <- c(0, cum)[b]                              # cumulative prob at left edge
    frac   <- ifelse(probs[b] > 0, (u - lo.cum) / probs[b], 0.5)
    breaks[b] + frac * (breaks[b + 1L] - breaks[b])
}

#' Forward CDF of a DP histogram (piecewise-linear inverse of \code{.qdphist})
#'
#' Maps a value to \eqn{\hat F(x) \in [0, 1]} from the DP-noised bin counts,
#' flat outside \code{[breaks[1], breaks[B + 1]]}.  Backs the data-adaptive
#' copula PIT, which transforms each covariate through its own already-released
#' DP marginal (\code{score.method = "adaptive"}).
#'
#' @param x Numeric values to evaluate.
#' @param counts Non-negative DP-noised bin counts (length B).
#' @param breaks Numeric bin boundaries (length B + 1).
#' @return Numeric vector the same length as \code{x}, in [0, 1].
#' @noRd
.pdphist <- function(x, counts, breaks) {
    B     <- length(counts)
    probs <- if (sum(counts) > 0) counts / sum(counts) else rep(1 / B, B)
    cum   <- c(0, cumsum(probs))                        # left-edge cum prob, length B+1
    x.c   <- pmin(pmax(x, breaks[1]), breaks[B + 1L])
    b     <- pmax(pmin(findInterval(x.c, breaks, rightmost.closed = TRUE), B), 1L)
    frac  <- (x.c - breaks[b]) / (breaks[b + 1L] - breaks[b])
    pmin(pmax(cum[b] + probs[b] * frac, 0), 1)   # a CDF: clamp floating-point overshoot to [0, 1]
}

# ---------------------------------------------------------------------------
# Categorical mechanism
# ---------------------------------------------------------------------------

#' DP perturbation and sampling for a categorical variable (Laplace mechanism)
#'
#' Computes DP-noisy counts over \code{levels} using Laplace(0, 2/eps) noise
#' (the replace-one count sensitivity of 2, review C1), converts to a
#' probability vector, and samples \code{length(x)} values.
#'
#' Note: the mechanism is Laplace (not exponential); the exponential mechanism
#' label used in some documentation is aspirational.  For the sensitivity and
#' privacy-budget accounting used here, Laplace on counts is equivalent.
#'
#' Combines the functionality of source \code{dp.cat} and
#' \code{sample.from.dp.cat} in \file{35-dp-full.r}: this version returns a
#' sampled character vector directly instead of an intermediate list.
#'
#' Same two calibrations as \code{.dp.hist.continuous}: a replace-one row swap
#' moves one unit between two levels, so the count vector has L1 sensitivity 2
#' (Laplace, \code{eps}) and L2 sensitivity sqrt(2) (Gaussian mu-GDP,
#' \code{mu}).
#'
#' the mechanism is separated from its post-processing, and that is deliberate.
#' \code{.dp.cat.probs} is the DP release: the noised, clipped, renormalised
#' probability vector.  \code{.dp.cat} draws a sample from it, which is
#' post-processing and spends no budget.  Until 0.2.0 the two were fused and
#' only the sample was returned, so the released probability vector could not be
#' observed at all: at any usable \code{mu} the multinomial sampling noise
#' swamps the DP noise, and a test written against the sample is measuring the
#' sampler.  A privacy channel nobody can look at is a privacy channel nobody
#' audits, which is how the left-packed-histogram defect survived a green
#' suite in the continuous mechanism next door.
#'
#' @param x Character or factor vector of observed values.
#' @param levels Character vector of all possible levels.
#' @param eps Privacy budget (epsilon > 0) for this marginal, Laplace path.
#' @param mu Per-channel mu for this marginal, Gaussian mu-GDP path.  Supply
#'   exactly one of \code{eps} and \code{mu}.
#' @return \code{.dp.cat.probs}: the DP probability vector over \code{levels}.
#' @noRd
.dp.cat.probs <- function(x, levels, eps = NULL, mu = NULL) {
    if (is.null(eps) == is.null(mu))
        stop("supply exactly one of eps (Laplace) or mu (Gaussian, mu-GDP)")
    L <- length(levels)
    v.counts <- as.numeric(table(factor(x, levels = levels)))
    v.noise <- if (is.null(mu)) .rlaplace(L, 0, 2 / eps)
               else stats::rnorm(L, 0, sqrt(2) / mu)
    v.counts.dp <- pmax(0, v.counts + v.noise)
    if (sum(v.counts.dp) > 0) v.counts.dp / sum(v.counts.dp) else rep(1 / L, L)
}

#' Sample a categorical column from its DP marginal (post-processing)
#'
#' @return Character vector of length \code{length(x)} sampled from
#'   \code{levels} according to the DP probability vector.
#' @noRd
.dp.cat <- function(x, levels, eps = NULL, mu = NULL) {
    sample(levels, length(x), replace = TRUE,
           prob = .dp.cat.probs(x, levels, eps = eps, mu = mu))
}

# ---------------------------------------------------------------------------
# Gaussian-copula helpers (continuous covariate dependence)
# ---------------------------------------------------------------------------

#' Bounded normal scores for the Gaussian copula (public PIT variants)
#'
#' Maps each observation through a public probability integral transform then the
#' normal quantile, giving row-local scores bounded by \eqn{\pm \Phi^{-1}(1 -
#' \tau)}.  The boundedness makes the copula correlation release's sensitivity
#' finite, and the transform is public, so one row's score does not move another's.
#'
#' @param x Numeric observations.
#' @param lo,hi Public lower and upper bounds of the covariate's box.
#' @param tau Tail clip in (0, 0.5), typically \code{1 / (2 * n)} for \code{"box"}
#'   or \code{0.005} for \code{"normal"} and \code{"adaptive"}.
#' @param method Public PIT variant.  \code{"box"} (box-uniform, legacy),
#'   \code{"normal"} (location-scale standardisation with public scale
#'   \code{sigma_pub = (hi - lo) / 8}, spreading the scores and lifting the
#'   correlation signal-to-noise at no privacy cost), or \code{"adaptive"}
#'   (normal-score PIT through the covariate's own released DP histogram,
#'   \code{qnorm(F_hat(x))}, dropping the location-scale assumption for skewed or
#'   bounded marginals).  The sensitivity tracks the score bound, so it is
#'   unchanged.  See \code{score.method} on the generators.
#' @param counts,breaks DP-noised bin counts and boundaries of the covariate's
#'   released histogram.  Required for \code{method = "adaptive"}, ignored
#'   otherwise.
#' @return Numeric vector of bounded normal scores, same length as \code{x}.
#' @noRd
.dp.copula.scores <- function(x, lo, hi, tau, method = c("box", "normal", "adaptive"),
                              counts = NULL, breaks = NULL) {
    method <- match.arg(method)
    if (method == "adaptive") {
        if (is.null(counts) || is.null(breaks))
            stop("score.method = \"adaptive\" requires counts and breaks")
        return(stats::qnorm(pmin(pmax(.pdphist(x, counts, breaks), tau), 1 - tau)))
    }
    if (method == "normal") {
        s.M <- stats::qnorm(1 - tau)
        return(pmin(pmax((x - (lo + hi) / 2) / ((hi - lo) / 8), -s.M), s.M))
    }
    x.clip <- pmin(pmax(x, lo), hi)
    u      <- pmin(pmax((x.clip - lo) / (hi - lo), tau), 1 - tau)
    stats::qnorm(u)
}

#' L2 sensitivity of the released vector vech(S) under one-row replacement
#'
#' the single definition.  \code{.dp.copula.corr} calibrates to it and the
#' sensitivity test in \code{test-mechanisms.R} measures against it, so the
#' calibration and the test guarding it cannot drift apart.
#'
#' Sharp: attained when both swapped rows sit on the clamp and are orthogonal.
#'
#' @param n Number of rows.
#' @param p Number of continuous covariates.
#' @param s.M Score bound (\code{qnorm(1 - tau)}).
#' @return The L2 sensitivity of \eqn{vech(S)}.
#' @noRd
.dp.copula.corr.delta <- function(n, p, s.M) s.M^2 * p / n

#' DP release of a correlation matrix via the second-moment Gaussian mechanism
#'
#' Releases the correlation matrix of bounded normal scores under
#' \eqn{(\varepsilon, \delta)}-DP.  The query is the second-moment matrix
#' \eqn{S = W^T W / n}, but what is released is \eqn{vech(S)}, its upper
#' triangle including the diagonal (the noise is drawn on \code{upper.tri(diag
#' = TRUE)} and mirrored).  A Gaussian mechanism is mu-GDP at
#' \eqn{\mu = \Delta / \sigma} with \eqn{\Delta} the L2 sensitivity of the
#' released vector, so the constant is the sensitivity of \eqn{vech(S)}:
#' replacing one row moves it by at most \eqn{\Delta_2 = M^2 p / n}, and that
#' is SHARP (attained when both swapped rows sit on the clamp and are
#' orthogonal).  Noise sd is \eqn{\Delta_2} times
#' \code{.dp.gm.multiplier(eps, delta)}.  Post-processing (eigen-clip to
#' PSD, rescale to unit diagonal) turns the noised matrix into a valid
#' correlation matrix at no privacy cost.
#'
#' Until 2026-07-14 this used \eqn{\Delta_2 = 2 M^2 p / n}, the Frobenius-norm
#' bound on the move of the full matrix.  Frobenius counts each off-diagonal
#' entry twice where \eqn{vech} counts it once, so that bound is sound but can
#' only over-noise: it over-noised this channel by a factor of 2.  Raw
#' second-moment rather than
#' Fisher-z: the latter's sensitivity is unbounded as the correlation
#' approaches +-1.  \eqn{S} is not mean-centred (privatising the means would
#' cost budget), so the rescaled entries are cosine similarities and the
#' recovered correlation is biased toward zero for off-centre score columns.
#'
#' Calibration: supply \code{mu} for the mu-GDP form (sigma = Delta_2 / mu),
#' which is what composes in quadrature with the marginal and coefficient
#' channels.  Supplying \code{(eps, delta)} instead reproduces the legacy
#' classical calibration and is confined to the retired
#' \code{sens.method = "local"} path.
#'
#' @param m.W Numeric \code{n x p} matrix of bounded normal scores.
#' @param eps Privacy budget for the correlation release (legacy path;
#'   ignored when \code{mu} is supplied).
#' @param delta DP delta for this Gaussian mechanism.
#' @param s.M Score bound (\code{qnorm(1 - tau)}), matching the scores.
#' @param seed Optional integer RNG seed; \code{NA} to skip (default).
#' @param mu Per-channel mu for the mu-GDP path.
#' @return A \code{p x p} symmetric PSD correlation matrix (unit diagonal).
#' @noRd
.dp.copula.corr <- function(m.W, eps, delta, s.M, seed = NA, mu = NULL) {
    if (!is.na(seed)) set.seed(seed)
    n <- nrow(m.W)
    p <- ncol(m.W)
    m.S     <- crossprod(m.W) / n
    s.Delta <- .dp.copula.corr.delta(n, p, s.M)
    s.sigma <- if (is.null(mu)) s.Delta * .dp.gm.multiplier(eps, delta)
               else s.Delta / mu
    # symmetric Gaussian noise on the upper triangle including the diagonal
    m.E   <- matrix(0, p, p)
    v.idx <- upper.tri(m.S, diag = TRUE)
    m.E[v.idx]           <- stats::rnorm(sum(v.idx), 0, s.sigma)
    m.E[lower.tri(m.E)]  <- t(m.E)[lower.tri(m.E)]
    m.S.dp <- m.S + m.E
    # nearest PSD by eigenvalue clipping, then rescale to unit diagonal
    l.eig <- eigen((m.S.dp + t(m.S.dp)) / 2, symmetric = TRUE)
    v.lam <- pmax(l.eig[["values"]], 1e-8)
    m.psd <- l.eig[["vectors"]] %*% (v.lam * t(l.eig[["vectors"]]))
    v.d   <- 1 / sqrt(diag(m.psd))
    m.R   <- v.d * m.psd * rep(v.d, each = p)
    (m.R + t(m.R)) / 2
}

#' Sample a covariate block from a DP Gaussian copula
#'
#' Draws correlated latent normals with covariance \code{m.R}, maps them to
#' uniforms, then through each column's DP-histogram quantile function.  The
#' marginals are exactly the DP histograms; the dependence is \code{m.R}.
#'
#' @param m.R A \code{p x p} correlation matrix (from \code{.dp.copula.corr}).
#' @param l.marg Named length-p list; entry j is
#'   \code{list(counts = <DP counts>, breaks = <bin edges>)}.  Order matches
#'   the rows and columns of \code{m.R}.
#' @param n Number of rows to draw.
#' @param seed Optional integer RNG seed; \code{NA} to skip (default).
#' @return A data frame with \code{n} rows and \code{p} named numeric columns.
#' @noRd
.sample.from.dp.copula <- function(m.R, l.marg, n, seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    p   <- ncol(m.R)
    m.L <- chol(m.R + 1e-10 * diag(p))                # R = t(L) %*% L
    m.Z <- matrix(stats::rnorm(n * p), n, p) %*% m.L   # rows ~ N(0, R)
    m.U <- stats::pnorm(m.Z)
    l.col <- Map(function(j, mg) .qdphist(m.U[, j], mg[["counts"]], mg[["breaks"]]),
                 seq_len(p), l.marg)
    names(l.col) <- names(l.marg)
    as.data.frame(l.col)
}

# ---------------------------------------------------------------------------
# Mixed-margin copula: categorical covariates
# ---------------------------------------------------------------------------

#' Bounded indicator scores for a categorical covariate (mixed-margin copula)
#'
#' Encodes a categorical covariate by its K-1 non-reference level indicators
#' (the first declared level is the reference), each mapped to a bounded van
#' der Waerden normal score computed from the released DP proportions.  The
#' scores are functions of public quantities only (the released probs and the
#' clamp), so they add no privacy cost, and every score lands in
#' \code{[-s.M, s.M]}, so appending these columns to the continuous score
#' matrix leaves the \code{.dp.copula.corr.delta} sensitivity form unchanged
#' at \code{p = p_total}.  A categorical column takes one of K fixed
#' winsorised values strictly inside the clamp, so the extended bound is a
#' sound upper bound NOT attained in the categorical sub-block (assumption A3
#' on \code{\link{gen.syn.dp.full}}).
#'
#' Mirrors \code{cat.indicator.scores} in \file{35-dp-full.r}; this version
#' takes the released \code{levels} and \code{probs} directly rather than the
#' source's \code{dp.cat} list.
#'
#' @param x Character or factor vector of observed values.
#' @param levels Character vector of all possible levels (level 1 is the
#'   reference).
#' @param probs Released DP probability vector over \code{levels} (from
#'   \code{.dp.cat.probs}).
#' @param s.M Score bound (\code{qnorm(1 - tau)}), matching the continuous
#'   scores.
#' @return Numeric \code{length(x) x (K - 1)} matrix of bounded scores, one
#'   column per non-reference level; zero columns when \code{K < 2}.
#' @noRd
.cat.indicator.scores <- function(x, levels, probs, s.M) {
    K <- length(levels)
    if (K < 2L) return(matrix(numeric(0), length(x), 0L))
    v.xc <- as.character(x)
    m.sc <- vapply(seq(2L, K), function(k) {
        p.k  <- probs[k]
        s.hi <- min(max(stats::qnorm(1 - p.k / 2),   -s.M), s.M)  # rows AT level k
        s.lo <- min(max(stats::qnorm((1 - p.k) / 2), -s.M), s.M)  # rows NOT at level k
        ifelse(v.xc == as.character(levels[k]), s.hi, s.lo)
    }, numeric(length(x)))
    matrix(m.sc, nrow = length(x))
}

#' Level proportions induced by a shifted back-map (mixed-margin copula)
#'
#' Evaluates the level proportions the back-map produces on a latent pool
#' \code{m.Z} at shifts \code{v.shift}: the argmax proportions under
#' \code{"probit"}, the mean softmax probabilities under \code{"logit"}.  The
#' reference level is the leading utility-zero column.
#'
#' @param m.Z Numeric \code{n x (K - 1)} matrix of latent normals.
#' @param v.shift Numeric length K-1 shift vector.
#' @param back.map \code{"probit"} or \code{"logit"}.
#' @return Numeric length-K probability vector over the levels.
#' @noRd
.mixed.cat.props <- function(m.Z, v.shift, back.map) {
    m.U <- cbind(0, sweep(m.Z, 2L, v.shift, "+"))
    if (back.map == "probit") {
        tabulate(max.col(m.U, ties.method = "first"), ncol(m.U)) / nrow(m.U)
    } else {
        m.P <- exp(m.U - apply(m.U, 1L, max))
        m.P <- m.P / rowSums(m.P)
        colMeans(m.P)
    }
}

#' Calibrate back-map shifts to the released DP marginal (mixed-margin copula)
#'
#' Fixed-point iteration on the K-1 shifts so that the back-mapped level
#' proportions reproduce the released DP proportions \code{v.probs} under the
#' chosen back-map, evaluated on a calibration pool \code{m.Z.cal} drawn from
#' the released sub-correlation.  Damped multiplicative updates (factor 0.7),
#' convergence at max abs deviation 1e-4, capped at 100 iterations: the same
#' sequence as the source's loop in \file{35-dp-full.r}, expressed as tail
#' recursion.
#'
#' @param m.Z.cal Numeric calibration pool, \code{n.cal x (K - 1)} latents.
#' @param v.probs Released DP probability vector over the K levels.
#' @param back.map \code{"probit"} or \code{"logit"}.
#' @return Numeric length K-1 shift vector.
#' @noRd
.mixed.cat.calibrate <- function(m.Z.cal, v.probs, back.map) {
    step <- function(v.shift, it) {
        if (it > 100L) return(v.shift)
        v.p <- .mixed.cat.props(m.Z.cal, v.shift, back.map)
        if (max(abs(v.p - v.probs)) < 1e-4) return(v.shift)
        step(v.shift + 0.7 * log(pmax(v.probs[-1L], 1e-6) /
                                 pmax(v.p[-1L], 1e-6)), it + 1L)
    }
    step(log(pmax(v.probs[-1L], 1e-6) / max(v.probs[1L], 1e-6)), 1L)
}

#' Assign levels from copula latents through the calibrated back-map
#'
#' \code{"probit"}: \code{argmax(0, Z_k + shift_k)}, deterministic, so the
#' assignment transfers the full released latent correlation; the natural
#' inverse of a Gaussian copula and consistent with the normal scores used to
#' fit it.  \code{"logit"}: one draw from \code{softmax(0, Z_k + shift_k)}
#' per row, whose independent choice noise attenuates the recovered
#' association (assumption A2 on \code{\link{gen.syn.dp.full}}).  Both are
#' post-processing of the released noisy correlation and spend no budget.
#'
#' @param m.Z Numeric \code{n x (K - 1)} matrix of copula latents.
#' @param v.shift Calibrated shifts (from \code{.mixed.cat.calibrate}).
#' @param v.levels Character vector of the K levels (level 1 = reference).
#' @param back.map \code{"probit"} or \code{"logit"}.
#' @return Factor of length \code{nrow(m.Z)} with levels \code{v.levels}.
#' @noRd
.mixed.cat.assign <- function(m.Z, v.shift, v.levels, back.map) {
    m.U <- cbind(0, sweep(m.Z, 2L, v.shift, "+"))
    v.k <- if (back.map == "probit") {
        max.col(m.U, ties.method = "first")
    } else {
        m.P   <- exp(m.U - apply(m.U, 1L, max))
        m.P   <- m.P / rowSums(m.P)
        m.cum <- t(apply(m.P, 1L, cumsum))
        1L + rowSums(m.cum < stats::runif(nrow(m.U)))
    }
    factor(v.levels[v.k], levels = v.levels)
}

#' De-attenuate the categorical-continuous entries of the released correlation
#'
#' The indicator-score fit estimates the POINT-BISERIAL correlation of a
#' two-point score against the continuous score, while the argmax back-map
#' needs the LATENT (biserial) correlation.  For a normal threshold at
#' released probability \code{p_k} the two differ by exactly
#' \deqn{c_k = \phi(\Phi^{-1}(p_k)) / \sqrt{p_k (1 - p_k)},}
#' the truncated-first-moment identity \eqn{E[1\{U > \tau\} Y] = \rho
#' \phi(\tau)}, so the categorical-continuous cross entries are divided by
#' \code{c_k} and the
#' matrix re-projected to a PSD correlation.  \code{c_k} reads the released
#' proportions only, so this is public post-processing: the released
#' \code{vech(S)} and its sensitivity are untouched (assumption A6 on
#' \code{\link{gen.syn.dp.full}}).  Within- and cross-categorical blocks are
#' left as fitted; their exact inversion is the multinomial-probit MLE.
#'
#' @param m.R The released correlation matrix, continuous columns first.
#' @param l.marg Named marginals list (categorical entries carry
#'   \code{levels} and \code{probs}).
#' @param v.cont.nms,v.cat.nms Covariate names in score-matrix order.
#' @return The de-attenuated PSD correlation matrix, unit diagonal.
#' @noRd
.mixed.deattenuate <- function(m.R, l.marg, v.cont.nms, v.cat.nms) {
    n.cont <- length(v.cont.nms)
    if (n.cont == 0L || length(v.cat.nms) == 0L) return(m.R)
    v.c <- unlist(Map(function(nm) {
        v.p <- l.marg[[nm]][["probs"]]
        if (length(v.p) < 2L) return(numeric(0))
        vapply(seq_len(length(v.p) - 1L) + 1L, function(k) {
            p.k <- min(max(v.p[k], 0.01), 0.99)
            stats::dnorm(stats::qnorm(p.k)) / sqrt(p.k * (1 - p.k))
        }, numeric(1))
    }, v.cat.nms))
    m.R <- Reduce(function(m, j) {
        col <- n.cont + j
        m[seq_len(n.cont), col] <- pmin(pmax(
            m[seq_len(n.cont), col] / v.c[[j]], -0.99), 0.99)
        m[col, seq_len(n.cont)] <- m[seq_len(n.cont), col]
        m
    }, seq_along(v.c), init = m.R)
    # nearest-PSD + unit diagonal, the same post-processing .dp.copula.corr ends on
    p <- ncol(m.R)
    l.eig <- eigen((m.R + t(m.R)) / 2, symmetric = TRUE)
    v.lam <- pmax(l.eig[["values"]], 1e-8)
    m.psd <- l.eig[["vectors"]] %*% (v.lam * t(l.eig[["vectors"]]))
    v.d   <- 1 / sqrt(diag(m.psd))
    m.R2  <- v.d * m.psd * rep(v.d, each = p)
    (m.R2 + t(m.R2)) / 2
}

#' Sample the whole covariate block from a mixed-margin DP Gaussian copula
#'
#' Draws correlated latents from the released extended correlation, then maps
#' the continuous columns through their DP-histogram quantiles and each
#' categorical's K-1 latent columns through the calibrated probit/logit
#' back-map.  Column order matches the extended score matrix: continuous in
#' \code{v.cont.nms} order, then each categorical's K-1 columns in
#' \code{v.cat.nms} order.  The back-map shifts are calibrated per categorical
#' on a 20000-row pool drawn from the released sub-correlation, so the
#' calibration reads public quantities only.  The latents are drawn from the
#' DE-ATTENUATED correlation (\code{.mixed.deattenuate} above, assumption A6).
#'
#' @param m.R The released \code{p_total x p_total} correlation matrix (from
#'   \code{.dp.copula.corr} on the extended score matrix).
#' @param l.marg Named marginals list covering every covariate: continuous
#'   entries are \code{list(counts, breaks)}, categorical entries are
#'   \code{list(levels, probs)}.
#' @param v.cont.nms Continuous covariate names, in score-matrix order.
#' @param v.cat.nms Categorical covariate names, in score-matrix order.
#' @param back.map \code{"probit"} or \code{"logit"}.
#' @param n Number of rows to draw.
#' @param seed Optional integer RNG seed; \code{NA} to skip (default).
#' @return A data frame with \code{n} rows: numeric columns for
#'   \code{v.cont.nms}, factor columns for \code{v.cat.nms}.
#' @noRd
.sample.from.dp.mixed.copula <- function(m.R, l.marg, v.cont.nms, v.cat.nms,
                                         back.map, n, seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    m.R <- .mixed.deattenuate(m.R, l.marg, v.cont.nms, v.cat.nms)
    p.total <- ncol(m.R)
    m.L <- chol(m.R + 1e-10 * diag(p.total))
    m.Z <- matrix(stats::rnorm(n * p.total), n, p.total) %*% m.L
    l.cont <- Map(function(nm, j) {
        mg <- l.marg[[nm]]
        .qdphist(stats::pnorm(m.Z[, j]), mg[["counts"]], mg[["breaks"]])
    }, v.cont.nms, seq_along(v.cont.nms))
    v.K     <- vapply(v.cat.nms,
                      function(nm) length(l.marg[[nm]][["levels"]]), integer(1))
    v.start <- length(v.cont.nms) + 1L +
        cumsum(c(0L, v.K - 1L))[seq_along(v.cat.nms)]
    l.cat <- Map(function(nm, s.start, K) {
        v.cols  <- seq.int(s.start, s.start + K - 2L)
        # calibrate shifts on a pool drawn from the released sub-correlation
        m.R.cc  <- m.R[v.cols, v.cols, drop = FALSE]
        m.Zc    <- matrix(stats::rnorm(20000L * (K - 1L)), 20000L, K - 1L) %*%
            chol(m.R.cc + 1e-10 * diag(K - 1L))
        v.shift <- .mixed.cat.calibrate(m.Zc, l.marg[[nm]][["probs"]], back.map)
        .mixed.cat.assign(m.Z[, v.cols, drop = FALSE], v.shift,
                          l.marg[[nm]][["levels"]], back.map)
    }, v.cat.nms, v.start, v.K)
    as.data.frame(c(l.cont, l.cat))
}
