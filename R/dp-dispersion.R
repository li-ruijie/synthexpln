# dp-dispersion.R -- where does sigma^2 come from, and .synth.response.
# Public functions:  dp.dispersion.public(), dp.dispersion.private(),
#                    print.dp.dispersion()
# Internal helpers:  .dp.has.dispersion(), .synth.response()
# Ported from lib/28-dp-projection.r lines 507-654.

# ---------------------------------------------------------------------------
# Dispersion provenance
# ---------------------------------------------------------------------------
# A GLM's likelihood has two parameters, the coefficient and the dispersion, and
# DP must cover both.  Until this file existed, only the coefficient passed through
# a DP channel here and the dispersion reached the release untouched.  Fixing
# that means saying, explicitly, where the dispersion came from, because the two
# defensible answers carry different privacy claims.
#
# `source` is required and is printed.  That is the point of the constructor.  A
# dispersion silently lifted from the sample -- which is what shipped -- now has
# to be written as
#
#     dp.dispersion.public(var(resid(fit)), source = "the sample")
#
# which is self-indicting on the page.  A design that makes the defect
# unwriteable without saying so is worth more than any single fix.

#' Declare a public dispersion (spends no privacy budget)
#'
#' The dispersion is already public: either external knowledge (a design
#' constant, prior literature) or a value the data owner has \emph{already
#' published}.  Journal articles routinely report \eqn{R^2} or a residual SD
#' alongside a coefficient table.  Post-processing permits any function of (DP
#' output, public information), so re-using such a value costs no budget and
#' discloses nothing new.  Where the owner published it themselves, the claim is
#' mu-GDP \emph{conditional} on that publication: the synthetic data adds nothing
#' the table did not already give away, and the paper must say so.
#'
#' This mode spends no privacy budget.  It leaves
#' \eqn{\hat\beta_\mathrm{DP}} as noisy as it would be with no dispersion
#' channel at all.
#'
#' @param sigma2 The public dispersion, a single positive finite number.
#' @param source \strong{Required.}  A string naming where this value comes
#'   from.  It is printed by \code{\link{print.dp.dispersion}} and is the
#'   point of the constructor: if you cannot write the sentence, the value is
#'   not public.
#' @return An object of class \code{dp.dispersion}.
#' @seealso \code{\link{dp.dispersion.private}} when nothing has been published.
#' @export
#' @examples
#' dp.dispersion.public(4.0, source = "residual SD of 2.0 cm reported in
#'                       the cohort's published baseline table")
dp.dispersion.public <- function(sigma2, source) {
    if (missing(source) || !is.character(source) || !nzchar(source))
        stop("source is REQUIRED: name where this PUBLIC dispersion comes from ",
             "(external constant, prior literature, or the owner's own ",
             "published table). If it came from the sample, it is not public.",
             call. = FALSE)
    if (!is.numeric(sigma2) || length(sigma2) != 1L || !is.finite(sigma2) ||
        sigma2 <= 0)
        stop("sigma2 must be a single positive finite number", call. = FALSE)
    structure(list(mode = "public", sigma2 = as.numeric(sigma2),
                   source = source),
              class = "dp.dispersion")
}

#' Declare a private dispersion, to be paid for from the budget
#'
#' Nothing has been published, so the dispersion is spent from the budget: it
#' rides along as a \eqn{(p+1)}-th aggregated coordinate in
#' subsample-and-aggregate, clipped into the public box \code{[lo, hi]}.  This
#' costs.  The L2 sensitivity of the scaled aggregate becomes \eqn{\sqrt{p+1}}
#' rather than \eqn{\sqrt{p}}, so \eqn{\hat\beta_\mathrm{DP}} is about 8\% noisier
#' at \eqn{p = 6}, and the standard errors inherit the noise.
#'
#' The box \emph{width} sets the noise scale, so a box read off the sample makes
#' the noise scale a function of the private data, which is the same defect one
#' level down.  \code{lo} and \code{hi} must be public.
#'
#' @param lo,hi Public bounds on \eqn{\sigma^2}, single finite numbers with
#'   \code{hi > lo}.
#' @param source \strong{Required.}  A string naming where the box comes from.
#' @return An object of class \code{dp.dispersion}.
#' @seealso \code{\link{dp.dispersion.public}} when the value is already public.
#' @export
#' @examples
#' dp.dispersion.private(lo = 0, hi = 25,
#'                       source = "instrument precision caps the residual SD
#'                                 at 5 units, so sigma^2 <= 25")
dp.dispersion.private <- function(lo, hi, source) {
    if (missing(source) || !is.character(source) || !nzchar(source))
        stop("source is REQUIRED: name where this PUBLIC BOX on the dispersion ",
             "comes from. The box width sets the noise scale, so a box read off ",
             "the sample makes the noise scale a function of the private data.",
             call. = FALSE)
    if (!is.numeric(lo) || !is.numeric(hi) || length(lo) != 1L ||
        length(hi) != 1L || !is.finite(lo) || !is.finite(hi) || hi <= lo)
        stop("lo and hi must be single finite numbers with hi > lo",
             call. = FALSE)
    structure(list(mode = "private", lo = as.numeric(lo), hi = as.numeric(hi),
                   source = source),
              class = "dp.dispersion")
}

#' Print a dispersion declaration
#'
#' @param x A \code{dp.dispersion} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.  Called for the side effect of printing the
#'   declared mode, value or box, and the \code{source} string.
#' @export
print.dp.dispersion <- function(x, ...) {
    cat(sprintf("<dp.dispersion: %s>\n", x[["mode"]]))
    if (identical(x[["mode"]], "public"))
        cat(sprintf("  sigma^2 = %g (no budget spent)\n", x[["sigma2"]]))
    else
        cat(sprintf("  public box [%g, %g] (spends budget: sqrt(p+1) sensitivity)\n",
                    x[["lo"]], x[["hi"]]))
    cat(sprintf("  source: %s\n", x[["source"]]))
    invisible(x)
}

# Families with a free dispersion parameter.  Poisson and binomial have none:
# their variance is determined by the mean, so they need no dispersion channel
# and cannot leak one.
#' Whether this family carries a free dispersion (internal)
#' @noRd
.dp.has.dispersion <- function(family)
    family %in% c("gaussian", "gamma", "quasi-gamma")

# ---------------------------------------------------------------------------
# .synth.response -- build the release from the DP output only
# ---------------------------------------------------------------------------
# this replaces .inject.virtual.y in the sound DP routes, and the difference
# governs whether the post-processing clause applies at all.
#
# .inject.virtual.y sets  y_virt = y + mu(beta_DP) - mu(beta_hat).  It retains
# the real y.  The projection then reproduces the real data's residuals in the
# release, scale AND shape.  Measured on the replication repository:
# sigma_syn = sigma_real to four significant figures across a ninefold range of
# true sigma, with the DP noise scale held constant, so beta_DP carried no
# information about sigma.  The shape flowed through too (heavy-tailed input
# gave kurtosis 7.0 -> 8.6), so the release did NOT depend on the data "only
# through beta_DP", which is the hypothesis the subsample-and-aggregate
# proposition needs for its post-processing clause.  A canary distinguishing
# game scored the resulting release at AUC 1.000 against a mu-GDP bound of 0.566
# at eps = 1: a perfect distinguisher, and flat in epsilon, because the leak
# did not go through the noise.
#
# .synth.response instead draws the response from the model's own family, at the
# mean implied by beta_DP and at the declared dispersion.  The result is a
# function of (beta_DP, dispersion, X, seed) and nothing else, so post-processing
# genuinely applies.  The real y appears nowhere.
#
# This costs the release its distributional fidelity to the empirical residuals
# and keeps its INFERENTIAL fidelity, which is the trade this package argues for
# throughout.  Note that even under a public dispersion the residuals must still
# be drawn from the family: the scale may be public, the empirical shape is a
# separate undeclared leak.

#' Draw the synthetic response from (beta_DP, dispersion) alone (internal)
#'
#' @param d Data frame supplying the covariates (and receiving the response).
#' @param y.nm Response column name.
#' @param x.nms Predictor names.
#' @param family Lowercase family string.
#' @param link Link string ("identity", "log", "logit").
#' @param v.beta.dp The DP coefficient: the only data-dependent input, and it
#'   is itself DP.
#' @param s.sigma2 The dispersion, from a \code{dp.dispersion} object.  NA for
#'   poisson and binomial, which have none.
#' @param seed RNG seed.
#' @return \code{d} with the response column replaced.
#' @noRd
.synth.response <- function(d, y.nm, x.nms, family, link,
                            v.beta.dp, s.sigma2, seed) {
    m.X <- stats::model.matrix(stats::reformulate(x.nms), d)
    v.b <- `names<-`(rep(0, ncol(m.X)), colnames(m.X))
    v.common <- intersect(names(v.beta.dp), colnames(m.X))
    v.b[v.common] <- v.beta.dp[v.common]

    v.eta <- as.numeric(m.X %*% v.b)

    # Guard the log link against overflow.  At a small epsilon the public box is
    # wide and beta_DP is very noisy, so exp(X beta_DP) can reach Inf and the
    # draw below then yields NaN.  Clamping the SYNTHETIC linear predictor is a
    # function of beta_DP and X only, so it stays post-processing and costs no
    # budget.  It exists to keep the arithmetic finite and NOT to restore
    # utility.  A release whose coefficients are this extreme still fails the
    # analyst's refit, and that must be reported rather than concealed.
    v.eta <- pmin(pmax(v.eta, -500), 500)

    v.mu <- switch(link,
        identity = v.eta,
        log      = exp(v.eta),
        logit    = stats::plogis(v.eta),
        stop("unsupported link in .synth.response: ", link))

    s.n <- nrow(m.X)
    s.p <- ncol(m.X)
    set.seed(seed + 11L)

    v.y <- switch(family,
        # Gaussian: draw, orthogonalise against the design, rescale to the target
        # dispersion.  The refit then recovers beta_DP AND sigma^2 exactly, so
        # the exact-fidelity property survives the change intact.
        gaussian = {
            v.e  <- stats::lm.fit(m.X, stats::rnorm(s.n))[["residuals"]]
            s.ss <- sum(v.e^2)
            if (s.ss <= 0) stop("degenerate residual basis in .synth.response")
            v.mu + v.e * sqrt(s.sigma2 * (s.n - s.p) / s.ss)
        },
        # Poisson and binomial carry no free dispersion: the family fixes the
        # variance given the mean, so the draw is complete on beta_DP alone.
        poisson  = stats::rpois(s.n, pmax(v.mu, .Machine[["double.eps"]])),
        binomial = stats::rbinom(s.n, 1L, pmin(pmax(v.mu, 1e-9), 1 - 1e-9)),
        # Gamma with dispersion phi: shape = 1/phi, rate = 1/(phi mu).
        gamma = ,
        `quasi-gamma` = {
            s.phi <- max(s.sigma2, 1e-8)
            stats::rgamma(s.n, shape = 1 / s.phi,
                          rate = 1 / (s.phi * pmax(v.mu,
                                                   .Machine[["double.eps"]])))
        },
        stop("unsupported family in .synth.response: ", family))

    `[[<-`(d, y.nm, value = v.y)
}
