# Sufficient-statistics projection generator for synthexpln.
# Public functions: projection(), projection.twopart().
# Internal helpers: .bootstrap.rows, .project.ols, .make.glm.family,
#                   .project.glm.
# Ported from lib/11-projection-generator.r (read-only source).

# +------------------------------------------------------------------------+
# | projection() -- public API                                             |
# +------------------------------------------------------------------------+

#' Sufficient-statistics projection generator
#'
#' Generates synthetic data whose refitted regression coefficients recover
#' the original estimates (the sufficient statistics projection theorem
#' of Li, 2026).  For Gaussian with the identity link the recovery is at
#' IEEE 754 double-precision machine epsilon
#' (\eqn{2^{-52} \approx 2.2 \times 10^{-16}}).  Other families are looser
#' by orders of magnitude and improve with \eqn{n}, Poisson being the
#' loosest since integer rounding limits it.  See the Recovery accuracy
#' section for measured figures.
#'
#' @section Recovery accuracy:
#' Measured over 24,000 runs spanning \eqn{n} from 100 to 5000, 2 to 8
#' predictors, independent and strongly correlated designs
#' (\eqn{\rho = 0.9}), and two response levels, 100 replicates per cell.
#' The quantity reported is the deviation relative to the size of the
#' coefficient vector being recovered,
#' \eqn{\max_j |\hat\beta_j^{\mathrm{refit}} - \hat\beta_j| / \max_j |\hat\beta_j|},
#' which, unlike a per-coefficient relative error, remains bounded when a
#' fitted coefficient approaches zero.
#'
#' \tabular{lrrr}{
#'   \strong{Family} \tab \strong{median} \tab \strong{95th pct} \tab \strong{worst} \cr
#'   gaussian/identity \tab 3.9e-16 \tab 1.3e-15 \tab 3.1e-15 \cr
#'   binomial/logit    \tab 1.4e-12 \tab 8.1e-09 \tab 4.2e-08 \cr
#'   Gamma/log         \tab 4.0e-07 \tab 7.0e-06 \tab 4.6e-05 \cr
#'   quasi-gamma/log   \tab 6.7e-07 \tab 9.7e-06 \tab 1.6e-02 \cr
#'   poisson/log       \tab 7.5e-03 \tab 1.4e-01 \tab 6.6e-01 \cr
#' }
#'
#' The Gaussian worst case is 13.9 units in the last place of double
#' precision.  It creeps up with \eqn{n} as rounding accumulates, the median
#' running 2.1e-16, 2.7e-16, 4.1e-16, and 5.8e-16 at \eqn{n} of 100, 500,
#' 2000, and 5000, which is approximately \eqn{n^{0.22}} and negligible either way.
#' The other families run the other way and improve as \eqn{n} grows, the
#' pooled Gamma and quasi-Gamma median falling from 2.8e-06 at
#' \eqn{n = 100} to 9.5e-08 at \eqn{n = 5000}, so no single fixed tolerance
#' describes them.  Poisson is governed by the count magnitude, as integer
#' rounding implies: at \eqn{n = 100} the median is 1.0e-01 at a mean count
#' of 2.3 and 4.0e-03 at a mean count of 10.2.  Quote the figure for your own
#' \eqn{n} and response scale rather than the median here, and see paper
#' Section 6.
#'
#' The preliminary synthetic data is produced by row-resampling; a
#' copula-based variant may replace this in a future release.  The
#' response is then projected onto the sufficient-statistics constraint
#' surface via \code{.project.ols} (Gaussian/identity) or
#' \code{.project.glm} (all other GLMs).
#'
#' @param d Data frame with the response and predictors.
#' @param formula Model formula.
#' @param family A GLM family object (default \code{gaussian()}).
#' @param seed Optional RNG seed.
#' @param base Preliminary synthetic data mechanism. \code{"bootstrap"}
#'   (row resampling, the default and the pre-0.9.0 behaviour) or
#'   \code{"copula"} (the stratified Gaussian-copula generator with
#'   Cornish-Fisher margins on the support-transformed scale, ported
#'   from the paper's analysis code).
#' @param design.exact Logical. Pin the synthetic predictor Gram
#'   \eqn{X^{*\prime}X^* = X'X} (the design-Gram projection of Li,
#'   2026). Requires \code{base = "copula"} and draws the
#'   stratum labels as an exact-count permutation, so for
#'   Gaussian-identity analysts the refitted coefficients, standard
#'   errors, and derived statistics are exact for every main-effects
#'   submodel, and for canonical GLM analysts the Fisher-weighted Gram
#'   is pinned to a fixed-point tolerance.
#' @param supports Optional named character vector of support tokens
#'   (\code{"real"}, \code{"positive"}, \code{"zeroinfl"},
#'   \code{"count"}, \code{"bounded"}) for the continuous variables,
#'   auto-detected when \code{NULL}. Copula base only.
#' @return An object of class \code{synthexpln}: a list with
#'   \code{$syn} (synthetic data frame), \code{$beta} (coefficient
#'   vector from the original fit), \code{$variant = "projection"},
#'   \code{$scope = "none"}, and \code{$design.exact} (logical,
#'   whether the design Gram was pinned).
#' @export
#' @importFrom MASS negative.binomial
#' @examples
#' d <- data.frame(y = rnorm(100), x = rnorm(100))
#' syn <- projection(d, y ~ x)
#' coef(lm(y ~ x, syn$syn))
#'
#' # design-Gram projection: standard errors exact as well
#' d2 <- data.frame(y = rnorm(120), x = rnorm(120),
#'                  g = factor(sample(c("a", "b"), 120, TRUE)))
#' syn2 <- projection(d2, y ~ x + g, seed = 1, design.exact = TRUE)
#' summary(lm(y ~ x + g, syn2$syn))$coef[, 2]
#' summary(lm(y ~ x + g, d2))$coef[, 2]
projection <- function(d, formula, family = stats::gaussian(), seed = NULL,
                       base = c("bootstrap", "copula"),
                       design.exact = FALSE, supports = NULL) {
    b.base.given <- !missing(base)
    base <- match.arg(base)
    if (design.exact && base == "bootstrap") {
        if (b.base.given)
            stop("design.exact = TRUE requires base = \"copula\".",
                 call. = FALSE)
        base <- "copula"
    }
    if (!is.null(seed)) set.seed(seed)

    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)

    # Extract and normalise family/link strings.  R family objects return
    # mixed-case names (e.g. "Gamma", "Negative Binomial") while the
    # internal helpers use lowercase snake strings from the source code.
    # Map: "Gamma" -> "gamma", "Negative Binomial" -> "neg.binomial",
    # "inverse.gaussian" / "InverseGaussian" -> "inv.gaussian".
    fam.nm.raw <- family[["family"]]
    fam.nm <- switch(fam.nm.raw,
        Gamma                     = "gamma",
        `Negative Binomial`       = "neg.binomial",
        `inverse.gaussian`        = "inv.gaussian",
        `Inverse Gaussian`        = "inv.gaussian",
        tolower(fam.nm.raw))      # gaussian, poisson, binomial, quasi*
    link.nm <- family[["link"]]

    # fit original to extract beta_hat (for the $beta slot)
    fit.orig <- if (fam.nm == "gaussian" && link.nm == "identity") {
        stats::lm(formula, d)
    } else {
        stats::glm(formula, family = family, data = d)
    }
    beta.hat <- stats::coef(fit.orig)

    # preliminary d.syn via the requested base
    d.syn <- if (base == "copula") {
        .copula.base(d, supports = supports, strata.exact = design.exact)
    } else {
        .bootstrap.rows(d)
    }

    # design-Gram projection: pin the predictor Gram before the response
    if (design.exact)
        d.syn <- .dg.project(d.syn, d, y.nm, x.nms, fam.nm, link.nm)

    # project the response onto the sufficient-statistics constraint
    d.syn[[y.nm]] <- if (fam.nm == "gaussian" && link.nm == "identity") {
        .project.ols(d.syn, d, y.nm, x.nms)
    } else {
        .project.glm(d.syn, d, y.nm, x.nms, fam.nm, link.nm)
    }

    structure(
        list(syn = d.syn, beta = beta.hat,
             variant = "projection", scope = "none",
             design.exact = design.exact),
        class = "synthexpln"
    )
}

# +------------------------------------------------------------------------+
# | projection.twopart() -- public API for hurdle (zero-inflated) data     |
# +------------------------------------------------------------------------+

#' Two-part (hurdle) sufficient-statistics projection
#'
#' Generates synthetic data for a zero-inflated (hurdle) response while
#' releasing the covariate design and the zero/positive indicator
#' unchanged.  Only the positive response values are synthesised, by
#' projecting the positive component onto the sufficient-statistics
#' constraint surface (the two-part projection proposition of supplement
#' S13 in Li, 2026).  Refitting the positive-component model
#' to the release reproduces the original coefficients, and, because the
#' positive design is released unchanged, their standard errors as well.
#' The positive block runs through the same GLM projection as
#' \code{\link{projection}} and inherits its accuracy.  This is not the
#' Gaussian identity case, the link being required to have positive support,
#' so recovery is approximate rather than at machine epsilon.  Over 800 runs
#' at Gamma/log, with \eqn{n} from 200 to 5000 and zero fractions of 0.3 and
#' 0.6, the deviation relative to the coefficient scale had a median of
#' 1.0e-06 and a worst case of 1.3e-05, decreasing from 1.6e-06 at
#' \eqn{n = 200} to around 6e-07 by \eqn{n = 2000}.  The zero pattern was
#' reproduced row for row in every one of those runs.
#'
#' The zeros are released exactly, so the synthetic zero pattern matches
#' the original row for row.  The positive block is given its residual
#' novelty by bootstrapping the positive responses among themselves,
#' after which the projection re-imposes the score constraint.  This is
#' the two-part analogue of \code{\link{projection}}; use
#' \code{\link{projection}} when there is no zero mass to hold fixed.
#'
#' @param d Data frame with the response and predictors.  The response
#'   must contain both zeros and positive values.
#' @param formula Model formula for the positive component.
#' @param family A GLM family object for the positive component.  Its
#'   link must have positive support (\code{"log"}, \code{"inverse"}, or
#'   \code{"1/mu^2"}) so the projection subsets to the positive rows.
#'   Default \code{stats::Gamma(link = "log")}.
#' @param seed Optional RNG seed.
#' @param design.exact Logical. Full synthesis with the design pinned
#'   per block (the two-part corollary of the design-Gram projection):
#'   the copula base draws an exact-count permutation jointly with the
#'   zero indicator, the positive rows are recoloured to the original
#'   positive-block Gram and the zero rows to the zero-row Gram, and
#'   the positive component is then projected on the pinned positive
#'   design.  The release is FULLY synthetic. The zero pattern is
#'   preserved per stratification cell rather than row for row, and
#'   the positive-block refit reproduces coefficients AND standard
#'   errors.  Default \code{FALSE}, the released-design behaviour.
#' @param supports Optional named character vector of support tokens
#'   for the continuous variables, auto-detected when \code{NULL}.
#'   \code{design.exact} only.
#' @return An object of class \code{synthexpln}: a list with \code{$syn}
#'   (the synthetic data frame, sharing the original design and zero
#'   pattern unless \code{design.exact} is set), \code{$beta} (the
#'   positive-component coefficient vector from the original fit),
#'   \code{$variant = "twopart"}, \code{$scope = "none"}, and
#'   \code{$design.exact} (logical).
#' @seealso \code{\link{projection}}
#' @export
#' @examples
#' set.seed(1)
#' n  <- 150
#' x  <- rnorm(n)
#' mu <- exp(0.5 + 0.4 * x)
#' y  <- ifelse(rbinom(n, 1, 0.6) == 1, rgamma(n, shape = 2, rate = 2 / mu), 0)
#' d  <- data.frame(y = y, x = x)
#' syn <- projection.twopart(d, y ~ x, family = Gamma("log"))
#' # the zero pattern is released unchanged
#' all((syn$syn$y > 0) == (d$y > 0))
#' # the positive-component coefficients are reproduced on a refit
#' d.syn.pos <- syn$syn[syn$syn$y > 0, ]
#' coef(glm(y ~ x, family = Gamma("log"), data = d.syn.pos))
projection.twopart <- function(d, formula, family = stats::Gamma(link = "log"),
                               seed = NULL, design.exact = FALSE,
                               supports = NULL) {
    if (!is.null(seed)) set.seed(seed)

    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)

    # normalise family/link to the lowercase snake strings the internal
    # helpers use, the same mapping as projection()
    fam.nm.raw <- family[["family"]]
    fam.nm <- switch(fam.nm.raw,
        Gamma               = "gamma",
        `Negative Binomial` = "neg.binomial",
        `inverse.gaussian`  = "inv.gaussian",
        `Inverse Gaussian`  = "inv.gaussian",
        tolower(fam.nm.raw))
    link.nm <- family[["link"]]

    # the positive block is projected only when the link subsets to the
    # positive rows; a real-support link would draw the zeros into the fit
    if (!link.nm %in% c("log", "inverse", "1/mu^2"))
        stop("projection.twopart() needs a positive-support link ",
             "(\"log\", \"inverse\", or \"1/mu^2\"); got \"", link.nm, "\".",
             call. = FALSE)

    v.y   <- d[[y.nm]]
    b.pos <- v.y > 0
    if (!any(b.pos) || all(b.pos))
        stop("projection.twopart() needs a response with both zeros and ",
             "positive values (a hurdle response).", call. = FALSE)

    # positive-component fit, the quantity a refit of the release reproduces
    beta.hat <- stats::coef(stats::glm(formula, family = family,
                                       data = d[b.pos, , drop = FALSE]))

    # design-Gram mode (two-part corollary): fully synthetic release with
    # the design pinned per block on the exact-count copula permutation
    if (design.exact) {
        d.syn <- .copula.base(d, supports = supports, strata.exact = TRUE)
        b.pos.syn <- d.syn[[y.nm]] > 0
        if (sum(b.pos.syn) != sum(b.pos))
            stop("design-gram twopart: positive count not preserved",
                 call. = FALSE)
        d.syn <- .dg.twopart(d.syn, d, y.nm, x.nms)
        d.syn[[y.nm]] <- .project.glm(d.syn, d, y.nm, x.nms, fam.nm, link.nm)
        stopifnot(all((d.syn[[y.nm]] > 0) == b.pos.syn))
        return(structure(
            list(syn = d.syn, beta = beta.hat,
                 variant = "twopart", scope = "none", design.exact = TRUE),
            class = "synthexpln"))
    }

    # release the design and the zero pattern unchanged; give the positive
    # block residual novelty by bootstrapping the positive responses, then
    # force the zeros to exactly zero
    d.syn <- d
    d.syn[[y.nm]][b.pos]  <- sample(v.y[b.pos], sum(b.pos), replace = TRUE)
    d.syn[[y.nm]][!b.pos] <- 0

    # project the positive component; .project.glm subsets to the positive
    # rows (zeros present) and writes the result back, leaving zeros untouched
    d.syn[[y.nm]] <- .project.glm(d.syn, d, y.nm, x.nms, fam.nm, link.nm)

    stopifnot(all((d.syn[[y.nm]] > 0) == b.pos))   # zero pattern preserved

    structure(
        list(syn = d.syn, beta = beta.hat,
             variant = "twopart", scope = "none", design.exact = FALSE),
        class = "synthexpln"
    )
}

# +------------------------------------------------------------------------+
# | .bootstrap.rows -- preliminary synthetic data via row resampling       |
# +------------------------------------------------------------------------+

#' Row bootstrap for preliminary synthetic data
#'
#' Resamples rows of \code{d} with replacement.  Preserves the joint
#' covariate distribution and handles factor columns correctly.
#'
#' @param d Data frame.
#' @return Data frame with \code{nrow(d)} rows, row names reset.
#' @noRd
.bootstrap.rows <- function(d) {
    d.out <- d[sample(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]
    `rownames<-`(d.out, NULL)
}

# +------------------------------------------------------------------------+
# | .project.ols -- OLS projection (closed-form Gaussian corollary)        |
# +------------------------------------------------------------------------+
# Closed-form Gaussian projection.  Decomposes y* into:
#   fitted:    X* beta_hat  (from original data -- encodes X'y)
#   residual:  c (I - H*) y_cop  (from bootstrap -- n-p-1 df of novelty)
# where c = sqrt(phi_orig * df / RSS_cop) rescales so that RSS* = RSS.
# This satisfies both score (C1) and saturated LL (C2) exactly.
#
# QR decomposition is used for the hat matrix because the synthetic
# design X* may be rank-deficient (sparse strata can drop factor levels).

#' OLS projection (internal)
#' @noRd
# disposition (3).  Theorem 'sufficient statistics projection' requires
# the (weighted) design to have full column rank.  A rank-deficient design is
# outside its hypotheses, so the response is a named refusal and not a ridge:
# regularising it would emit a release the theorem does not cover, quietly, and
# which branch was taken would depend on the units of the covariates.  The
# 1e-10 added below the call sites stays absolute, a numerical nicety for the
# admissible case.  Scaling it to max(abs(.)) was measured to cost the
# projection its machine-zero exactness on the binomial branch and on the
# zero-inflated quasi-gamma designs.
.assert.full.rank <- function(m.wX, s.where) {
    s.rk <- qr(m.wX)[["rank"]]
    if(s.rk < ncol(m.wX))
        stop(sprintf(paste0("projection (%s): the weighted design is rank ",
                            "deficient, rank %d of %d columns.  Theorem ",
                            "'sufficient statistics projection' requires full ",
                            "column rank."), s.where, s.rk, ncol(m.wX)))
    invisible(TRUE)
}

.project.ols <- function(d.syn, d.orig, y.nm, x.nms) {
    fml <- stats::reformulate(x.nms, response = y.nm)

    # fit the original model to obtain beta_hat and dispersion phi
    mo.orig <- stats::lm(fml, d.orig)
    v.beta  <- stats::coef(mo.orig)
    s.phi   <- sum(stats::residuals(mo.orig)^2) /
               (nrow(d.orig) - length(v.beta))

    # synthetic design matrix and bootstrap response
    m.Xs    <- stats::model.matrix(fml, d.syn)
    v.y.cop <- d.syn[[y.nm]]

    # align beta_hat to synthetic X* columns.  The row bootstrap
    # may generate strata with fewer factor levels than the original,
    # so some columns in X* may be absent from beta_hat.  Coefficients
    # for missing levels are set to zero (conservative).
    v.nms.Xs <- colnames(m.Xs)
    v.nms.b  <- names(v.beta)
    v.common <- intersect(v.nms.Xs, v.nms.b)
    v.beta.aligned <- `names<-`(rep(0, ncol(m.Xs)), v.nms.Xs)
    v.beta.aligned[v.common] <- v.beta[v.common]

    # fitted component: X* beta_hat (carries the X'y information)
    v.fit.s <- as.numeric(m.Xs %*% v.beta.aligned)

    # residual component: (I - H*) y_cop via QR factorisation.
    # H* = X* (X*'X*)^{-1} X*' is the hat matrix of the synthetic design.
    # QR handles rank deficiency (rank < ncol) gracefully.
    l.qr    <- qr(m.Xs)
    s.rank  <- l.qr[["rank"]]
    v.Hy    <- qr.fitted(l.qr, v.y.cop)   # qr.fitted is base, not stats::
    v.e.cop <- v.y.cop - v.Hy              # (I - H*) y_cop

    # scale factor c so that RSS* = phi * df  (constraint C2).
    s.rss.cop <- sum(v.e.cop^2)
    s.df      <- nrow(m.Xs) - s.rank
    s.scale   <- `if`(s.rss.cop > 0 && s.df > 0,
                       sqrt(s.phi * s.df / s.rss.cop), 1)

    v.fit.s + v.e.cop * s.scale
}

# +------------------------------------------------------------------------+
# | .make.glm.family -- GLM family object constructor                      |
# +------------------------------------------------------------------------+
# Maps family/link string pairs to R family objects.
# The fallback for unrecognised families is quasi(log, mu^2), which
# gives quasi-likelihood with V(mu) = mu^2 (gamma-like variance).

#' GLM family object constructor (internal)
#' @noRd
.make.glm.family <- function(family, link) {
    switch(family,
        gaussian       = stats::gaussian(link = link),
        poisson        = stats::poisson(link = link),
        binomial       = stats::binomial(link = link),
        gamma          = stats::Gamma(link = link),
        `neg.binomial` = MASS::negative.binomial(theta = 1, link = link),
        `inv.gaussian` = stats::inverse.gaussian(link = link),
        stats::quasi(link = link, variance = "mu^2"))
}

# +------------------------------------------------------------------------+
# | .project.glm -- GLM projection (general case)                          |
# +------------------------------------------------------------------------+
# Core projection engine.  Given bootstrap output y_cop and the original
# data, find y* satisfying:
#   (C1) X*' w (y* - mu_hat) = 0     (score equation at beta_hat)
#   (C2) dispersion phi* = phi        (saturated LL / Pearson scaling)
#
# Three code paths, selected by family and link:
#   1. Real support (Gaussian non-identity link): weighted linear projection
#      on y-scale + RSS rescaling.
#   2. Positive support, quasi-gamma: alternating y-scale projection with
#      positivity clipping (Newton on log-scale would collapse Pearson phi).
#   3. Positive support, standard families: Newton on z = log(y) scale.
#      y = exp(z) > 0 by construction.  Poisson/NB get probabilistic
#      rounding afterwards.
#
# Ported verbatim from project.glm() in lib/11-projection-generator.r
# (lines 192-455).  Internal function names changed to dotted-private;
# stats:: prefix added throughout; source() calls removed.

#' GLM projection (internal)
#' @noRd
.project.glm <- function(d.syn, d.orig, y.nm, x.nms, family, link) {
    fml <- stats::reformulate(x.nms, response = y.nm)

    # dispatch Gaussian/identity to the closed-form OLS path
    if (family == "gaussian" && link == "identity")
        return(.project.ols(d.syn, d.orig, y.nm, x.nms))

    # for links that require y > 0 (log, inverse, 1/mu^2), exclude
    # zero observations from the original fit (ZI zeros are handled
    # separately by the copula's categorical bootstrap)
    b.has.zeros <- any(d.orig[[y.nm]] == 0)
    d.orig.fit <- `if`(link %in% c("log", "inverse", "1/mu^2") && b.has.zeros,
                       d.orig[d.orig[[y.nm]] > 0, ], d.orig)

    # fit the original model to extract beta_hat and phi.
    # quasi families need an explicit start (intercept = 1, slopes = 0)
    # because glm() cannot auto-initialise quasi variance functions.
    # standard families try default initialisation first, then fall back
    # to an eta-based start (intercept = g(y_bar), slopes = 0).
    mo.orig <- tryCatch({
        fam.obj <- .make.glm.family(family, link)
        b.quasi <- inherits(fam.obj, "family") &&
                   fam.obj[["family"]] == "quasi"
        if (b.quasi) {
            s.np <- ncol(stats::model.matrix(fml, d.orig.fit))
            stats::glm(fml, d.orig.fit, family = fam.obj,
                start = c(1, rep(0, s.np - 1)),
                control = stats::glm.control(maxit = 200))
        } else {
            tryCatch(
                stats::glm(fml, d.orig.fit, family = fam.obj,
                    control = stats::glm.control(maxit = 200)),
                error = function(e) {
                    s.np  <- ncol(stats::model.matrix(fml, d.orig.fit))
                    s.eta0 <- fam.obj[["linkfun"]](
                        mean(d.orig.fit[[y.nm]]))
                    stats::glm(fml, d.orig.fit, family = fam.obj,
                        start = c(s.eta0, rep(0, s.np - 1)),
                        control = stats::glm.control(maxit = 200))
                })
        }
    }, error = function(e) NULL)

    # if the original fit fails, return bootstrap values unmodified
    if (is.null(mo.orig)) return(d.syn[[y.nm]])

    # extract beta_hat (p coefficients) and dispersion phi from
    # the original fit.  These are the target quantities: projection
    # constructs y* so that refitting gives these back exactly.
    v.beta     <- stats::coef(mo.orig)
    s.phi.orig <- stats::summary.glm(mo.orig)[["dispersion"]]

    # build synthetic design matrix.  For ZI data with log/inverse link,
    # subset to positive observations (zeros handled by the bootstrap).
    # v.idx.fit is reused by the final write-back, which must target
    # exactly the rows the projection ran on.
    v.idx.fit <- `if`(link %in% c("log", "inverse", "1/mu^2") && b.has.zeros,
                      which(d.syn[[y.nm]] > 0), seq_len(nrow(d.syn)))
    d.syn.fit <- d.syn[v.idx.fit, , drop = FALSE]
    m.Xs    <- stats::model.matrix(fml, d.syn.fit)
    v.y.cop <- d.syn.fit[[y.nm]]
    s.n     <- nrow(m.Xs)
    s.p     <- ncol(m.Xs)

    # support flags: b.positive means y > 0 required (Gamma, Poisson, etc.);
    # b.integer means y must be a non-negative integer (Poisson, NB).
    b.positive <- family %in% c("gamma", "inv.gaussian", "poisson",
                                "neg.binomial") ||
                  (grepl("quasi", family) && link %in% c("log", "inverse"))
    b.integer  <- family %in% c("poisson", "neg.binomial")
    # b.binary means y in [0, 1] required (binomial); projected onto its own
    # [0, 1] box path below, not the real-support or positive-support paths.
    b.binary   <- family == "binomial"

    # compute mu_hat(X*) = g^{-1}(X* beta_hat) -- the fitted values
    # from the original model evaluated at synthetic covariates.
    # Fallback: if predict() fails (e.g. missing factor levels), compute
    # mu manually via aligned beta and the inverse link.
    v.mu.s <- tryCatch(
        stats::predict(mo.orig, newdata = d.syn.fit, type = "response"),
        error = function(e) {
            v.nms.Xs <- colnames(m.Xs)
            v.nms.b  <- names(v.beta)
            v.common <- intersect(v.nms.Xs, v.nms.b)
            v.beta.a <- `names<-`(rep(0, ncol(m.Xs)), v.nms.Xs)
            v.beta.a[v.common] <- v.beta[v.common]
            family(mo.orig)[["linkinv"]](as.numeric(m.Xs %*% v.beta.a))
        })

    # Numerical conditioning: for positive-support families,
    # the inverse link can produce non-positive mu when extrapolated
    # covariates push eta outside the training range (e.g. Gamma with
    # inverse link when eta crosses zero).  Replace invalid mu with
    # the median of the valid fitted values.
    b.valid <- rep(TRUE, s.n)
    if (b.positive) {
        b.valid <- is.finite(v.mu.s) & v.mu.s > 1e-6
        if (!all(b.valid)) {
            v.mu.s[!b.valid] <- stats::median(v.mu.s[b.valid])
        }
    }

    # linear predictor eta = g(mu) for score weight computation
    v.eta.s <- family(mo.orig)[["linkfun"]](pmax(v.mu.s, 1e-10))

    # score weights w_i = (dmu/deta)_i / V(mu_i), evaluated at beta_hat.
    # The GLM score equation is X' diag(w) (y - mu) = 0.
    v.dmu.deta <- family(mo.orig)[["mu.eta"]](v.eta.s)
    v.V        <- family(mo.orig)[["variance"]](pmax(v.mu.s, 1e-10))
    v.w        <- v.dmu.deta / v.V

    # quasi families (V(mu) = mu^2 with log link) give w_i = 1/mu_i and
    # score = X'(y-mu)/mu.  For these, Newton on log-scale collapses
    # Pearson dispersion, so we use y-scale alternating projection instead.
    # Standard families (Gamma, Poisson, NB) use log-scale Newton, which
    # respects the support constraint y > 0 automatically.
    b.quasi.gamma <- grepl("quasi", family)

    # ---- PATH B: BINOMIAL (project onto [0, 1] with exact score match) ----
    # For binomial, matching the score X*' W (y* - mu_hat) = 0 recovers
    # beta_hat exactly, and the standard errors follow from the Fisher
    # information X*' W X* alone, so no dispersion rescaling is needed
    # (phi = 1 is fixed).  The response must stay in [0, 1] (a valid binomial
    # proportion).  One weighted score correction zeroes the score, leaving a
    # score-orthogonal residual r.  Scaling r by any scalar keeps the score at
    # zero, so shrink it by the largest s in [0, 1] that keeps mu + s r inside
    # [0, 1] element-wise.  The result is an exact, valid fractional response;
    # a Bernoulli-rounded {0, 1} release would trade this exactness for integer
    # support (the irreducible sampling variance of the V5 DP-logistic mode).
    if (b.binary) {
        m.wX    <- sweep(m.Xs, 1, v.w, `*`)
        .assert.full.rank(m.wX, "binomial")
        m.XtwX  <- crossprod(m.wX)
        diag(m.XtwX) <- diag(m.XtwX) + 1e-10
        v.score <- as.numeric(crossprod(m.wX, v.y.cop - v.mu.s))
        v.corr  <- as.numeric(m.wX %*% solve(m.XtwX, v.score))
        v.r     <- (v.y.cop - v.corr) - v.mu.s        # score-orthogonal residual
        # largest s in [0, 1] keeping mu + s r within [0, 1] element-wise
        v.s.up  <- ifelse(v.r >  1e-12, (1 - v.mu.s) / v.r, Inf)
        v.s.lo  <- ifelse(v.r < -1e-12, (0 - v.mu.s) / v.r, Inf)
        s.shr   <- max(0, min(1, min(v.s.up), min(v.s.lo)))
        v.y.proj <- pmin(pmax(v.mu.s + s.shr * v.r, 0), 1)
    # ---- PATH 1: real support (Gaussian with non-identity link) ----------
    # Score constraint X*' W (y* - mu) = 0 is linear in y*, so one
    # weighted least-squares correction suffices.  Then rescale
    # residuals to match RSS (constraint C2: saturated LL).
    } else if (!b.positive) {
        # weighted design matrix and score correction
        m.wX    <- sweep(m.Xs, 1, v.w, `*`)
        .assert.full.rank(m.wX, "non-canonical")
        m.XtwX  <- crossprod(m.wX)
        diag(m.XtwX) <- diag(m.XtwX) + 1e-10  # Tikhonov regularisation
        v.score <- as.numeric(crossprod(m.wX, v.y.cop - v.mu.s))
        # minimum-norm correction: project score residual back onto col(wX)
        v.corr  <- as.numeric(m.wX %*% solve(m.XtwX, v.score))
        v.y.proj <- v.y.cop - v.corr

        # constraint C2: rescale residuals so RSS* = RSS_orig
        v.e.proj   <- v.y.proj - v.mu.s
        s.rss.orig <- sum(stats::residuals(mo.orig)^2)
        s.rss.proj <- sum(v.e.proj^2)
        if (s.rss.proj > 0) {
            s.scale  <- sqrt(s.rss.orig / s.rss.proj)
            v.y.proj <- v.mu.s + v.e.proj * s.scale
        }

    # ---- PATH 2: QUASI-GAMMA (y-scale alternating projection) -----------
    # For quasi(log, mu^2), the score X' diag(1/mu) (y - mu) = 0 is
    # linear in y but positivity y > 0 is not guaranteed after a single
    # projection step.  Alternate: (a) project score to zero, (b) clip
    # negatives to 1e-6, until all y > 0.  This preserves Pearson
    # dispersion better than Newton on log-scale for quasi families.
    } else if (b.quasi.gamma) {
        m.wX    <- sweep(m.Xs, 1, v.w, `*`)
        .assert.full.rank(m.wX, "quasi-gamma")
        m.XtwX  <- crossprod(m.wX)
        diag(m.XtwX) <- diag(m.XtwX) + 1e-10
        m.XtwX.inv <- solve(m.XtwX)
        v.y.cur <- v.y.cop

        s.maxiter <- 50
        s.iter    <- 0
        repeat {
            s.iter <- s.iter + 1
            # score correction step: project onto C1
            v.score <- as.numeric(crossprod(m.wX, v.y.cur - v.mu.s))
            v.corr  <- as.numeric(m.wX %*% (m.XtwX.inv %*% v.score))
            v.y.cur <- v.y.cur - v.corr

            # positivity enforcement: clip negatives, then re-check score
            s.n.neg <- sum(v.y.cur <= 0)
            if (s.n.neg == 0 || s.iter >= s.maxiter) break
            v.y.cur <- pmax(v.y.cur, 1e-6)
        }

        # constraint C2: rescale Pearson residuals to match original phi.
        v.pearson  <- (v.y.cur - v.mu.s) /
                      sqrt(family(mo.orig)[["variance"]](pmax(v.mu.s, 1e-10)))
        s.phi.cur  <- sum(v.pearson^2) / max(s.n - s.p, 1)
        if (s.phi.cur > 1e-10) {
            s.c <- sqrt(s.phi.orig / s.phi.cur)
            v.y.cur <- v.mu.s + (v.y.cur - v.mu.s) * s.c
            v.y.cur <- pmax(v.y.cur, 1e-6)
        }
        v.y.proj <- v.y.cur

    # ---- PATH 3: POSITIVE support (log-scale Newton iteration) ----------
    # The log-scale Newton step.  Work on z = log(y), solve
    #   F(z) := X*' w (exp(z) - mu_hat) = 0
    # by Newton's method.  Since y = exp(z) > 0 for all real z, the
    # support constraint is satisfied automatically at every iterate.
    #
    # Jacobian: J = X*' diag(w * exp(z)) = X*' diag(w * y).
    # Minimum-norm update: dz = diag(wy) X (X' diag(wy)^2 X)^{-1} score.
    # Convergence is fast (typically 3-5 iterations).
    } else {
        v.z <- log(pmax(v.y.cop, 1e-10))

        s.maxiter <- 30
        s.tol     <- 1e-8
        s.iter    <- 0
        repeat {
            s.iter <- s.iter + 1
            v.y <- exp(v.z)

            # evaluate the score: X*' w (y - mu_hat)
            v.score <- as.numeric(crossprod(m.Xs, v.w * (v.y - v.mu.s)))
            if (max(abs(v.score)) < s.tol || s.iter > s.maxiter) break

            # Gauss-Newton step: solve J J' lambda = score for the
            # minimum-norm update dz = J' lambda
            v.wy  <- v.w * v.y
            m.WYX <- sweep(m.Xs, 1, v.wy, `*`)
            m.JJt <- crossprod(m.WYX)
            diag(m.JJt) <- diag(m.JJt) * (1 + 1e-10)  # regularise

            v.lambda <- tryCatch(solve(m.JJt, v.score),
                                 error = function(e) rep(0, s.p))
            v.dz <- as.numeric(m.WYX %*% v.lambda)

            # damped step: cap max |dz| at 3 to prevent overshooting
            s.step <- min(1, 3 / (max(abs(v.dz)) + 1e-16))
            v.z <- v.z - s.step * v.dz
        }

        # constraint C2: rescale Pearson residuals to match phi_orig.
        v.y.newton   <- exp(v.z)
        v.sd.s       <- sqrt(family(mo.orig)[["variance"]](pmax(v.mu.s, 1e-10)))
        v.pearson    <- (v.y.newton - v.mu.s) / v.sd.s
        s.phi.newton <- sum(v.pearson^2) / max(s.n - s.p, 1)
        if (s.phi.newton > 1e-10) {
            s.c <- sqrt(s.phi.orig / s.phi.newton)
            v.y.proj <- v.mu.s + v.sd.s * v.pearson * s.c
            # enforce y > 0 after Pearson rescaling (rescaling can push
            # observations with small mu below zero)
            if (b.positive) v.y.proj <- pmax(v.y.proj, 1e-10)
        } else {
            v.y.proj <- v.y.newton
        }

        # Poisson / NB: probabilistic rounding to integers.
        # P(ceil) = frac(y), P(floor) = 1-frac(y), so E[round(y)] = y.
        # This introduces ~2.3% residual error in the score constraint
        # (unavoidable for integer-valued families).
        # NOTE: source uses set.seed(42) here for reproducibility; we
        # rely instead on the caller's seed (passed via projection(seed=))
        # so that output is deterministic without a hidden global seed
        # mutation -- a cleaner contract for a package function.
        if (b.integer) {
            v.y.cont <- v.y.proj
            v.frac   <- v.y.cont - floor(v.y.cont)
            v.y.proj <- pmax(ifelse(stats::runif(s.n) < v.frac,
                                    ceiling(v.y.cont),
                                    floor(v.y.cont)), 0)
        }

        # for observations with invalid mu (numerical conditioning),
        # keep the original bootstrap values unchanged
        if (!all(b.valid))
            v.y.proj[!b.valid] <- v.y.cop[!b.valid]
    }

    # write projected values back into the full synthetic vector, onto
    # the rows the projection ran on (v.idx.fit).  For ZI data the rows
    # outside v.idx.fit keep their zeros from the bootstrap.
    v.result <- d.syn[[y.nm]]
    v.result[v.idx.fit] <- v.y.proj
    v.result
}
