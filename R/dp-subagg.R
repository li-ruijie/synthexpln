# dp-subagg.R -- worst-case log-link DP via subsample-and-aggregate.
# Public function:  gen.syn.dp.projected.subagg()
# Internal helpers: .subsample.aggregate.beta.loglink(),
#                   .smooth.sens.beta.loglink(),
#                   .smooth.sens.beta.loglink.std()
# Ported from lib/28-dp-projection.r (subsample.aggregate.beta.loglink,
# gen.syn.dp.projected.subagg, smooth.sens.beta.loglink,
# smooth.sens.beta.loglink.std).

# +------------------------------------------------------------------------+
# | .subsample.aggregate.beta.loglink -- worst-case DP beta (internal)     |
# +------------------------------------------------------------------------+
# Partition rows into m disjoint blocks by data-independent position, fit
# the GLM per block, aggregate each coordinate by a Widened Winsorized Mean
# (clip to the public box [lo, hi], average), and add Gaussian noise
# calibrated to the worst-case per-coordinate sensitivity (hi - lo)/m via the
# scaled query, whose L2 sensitivity is sqrt(q).
#
# Why this is sound under full-row replacement, which is the relation that
# matters once the design is sensitive: the block partition is by row position
# under a fixed public seed, so it is data-independent.  Changing one row
# perturbs exactly one block.  That block's coefficient is clipped into the
# public box, so it moves by at most hi_j - lo_j, and the mean over m blocks
# moves by at most s_j = (hi_j - lo_j)/m.  No design constant and no
# data-dependent quantity enters, and it holds whether the neighbour changed x,
# y or both.  This is why subsample-and-aggregate needs neither a public design
# nor a bounded response, and it is the reason the mechanism exists
# (Nissim, Raskhodnikova and Smith 2007 designed it for statistics with no
# bounded global sensitivity).
#
# The partition is data-independent only when the row set is fixed.  A caller
# that filters to a DATA-DEPENDENT subset first (the zero-inflated positive
# block, whose row count |P| = #positives is itself sensitive) must instead
# partition the fixed full N and pass it as `block`; partitioning the filtered
# subset makes the partition a function of |P|, and a one-row replacement that
# toggles a row's positivity then reshuffles the whole partition and breaks the
# s_j bound (attains ~3x it; test-dp-zi.R).
#
# disp = TRUE appends sigma^2 as a (p+1)-th aggregated coordinate with its own
# public box, so q = p + 1 and the dispersion is paid for out of the same
# channel.  disp = FALSE is the default, so the families with no free dispersion
# stay p-dimensional.

#' Subsample-and-aggregate worst-case DP coefficient (internal)
#'
#' @param mu Per-channel mu for the mu-GDP path.  When supplied, the noise is
#'   sigma.Q = sqrt(q) / mu and the release composes in quadrature with the
#'   other channels.  When NULL the (eps, delta) pair is calibrated through the
#'   exact Gaussian-DP dual instead.
#' @param block Optional length-\code{nrow(d)} block assignment (values in
#'   1..m) that replaces the internal seeded partition, letting the caller
#'   partition a fixed superset of the rows.  The zero-inflated generator
#'   partitions the full N and restricts the assignment to the positive rows,
#'   which keeps the partition data-independent even when a one-row replacement
#'   toggles a row's positivity and so changes |P|.
#' @param disp Append sigma^2 as an aggregated coordinate (gaussian only).
#' @param disp.lo,disp.hi public box on sigma^2.  Required when disp = TRUE.
#' @noRd
.subsample.aggregate.beta.loglink <- function(d, formula, family,
                                              fam.nm, link.nm,
                                              eps, delta, m = NULL,
                                              clip.lo, clip.hi, seed,
                                              mu = NULL, block = NULL,
                                              disp = FALSE,
                                              disp.lo = NULL, disp.hi = NULL) {
    if (disp && (is.null(disp.lo) || is.null(disp.hi)))
        stop("disp.lo and disp.hi (public bounds on sigma^2) are required ",
             "when disp = TRUE")
    if (disp && fam.nm != "gaussian")
        stop("disp = TRUE is defined only for the gaussian family ",
             "(OLS dispersion)")
    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)
    fml   <- stats::reformulate(x.nms, response = y.nm)
    n     <- nrow(d)

    fit.block <- function(dd) tryCatch({
        if (fam.nm == "gaussian" && link.nm == "identity") {
            stats::lm(fml, dd)
        } else if (grepl("quasi", fam.nm)) {
            np <- ncol(stats::model.matrix(fml, dd))
            stats::glm(fml, dd, family = family, start = c(1, rep(0, np - 1)),
                       control = stats::glm.control(maxit = 200))
        } else {
            stats::glm(fml, dd, family = family)
        }
    }, error = function(e) NULL)

    mo.ref <- fit.block(d)
    if (is.null(mo.ref)) stop("reference fit failed")
    beta.hat <- stats::coef(mo.ref)
    v.names  <- names(beta.hat)
    p        <- length(v.names)
    v.qnames <- if (disp) c(v.names, ".sigma2") else v.names
    q        <- length(v.qnames)

    align.clip <- function(cl) {
        if (is.null(names(cl))) return(`names<-`(rep_len(cl, p), v.names))
        v <- cl[v.names]
        if (any(is.na(v))) stop("named clip bounds must cover every coefficient")
        v
    }
    v.lo <- align.clip(clip.lo); v.hi <- align.clip(clip.hi)
    if (disp) {
        v.lo <- `names<-`(c(v.lo, disp.lo), v.qnames)
        v.hi <- `names<-`(c(v.hi, disp.hi), v.qnames)
    }
    if (any(v.hi <= v.lo)) stop("clip.hi must exceed clip.lo for every coordinate")
    v.center <- (v.lo + v.hi) / 2

    # Block-count rule.  Self-default only on a self-partitioned call: a
    # supplied block assignment means the caller partitioned a fixed superset,
    # so nrow(d) here is a restricted, sensitive row count the rule must not
    # read.  On that path m stays a required argument, computed upstream from
    # published constants.
    if (is.null(m)) {
        if (!is.null(block))
            stop("m is required when block is supplied.  The block-count rule ",
                 "reads nrow(d), which on a restricted row set is sensitive; ",
                 "compute m upstream from published constants ",
                 "(dp.subagg.blocks) and pass it.", call. = FALSE)
        m <- dp.subagg.blocks(n, p)
    }

    # Partition the rows into m blocks.  When the caller supplies `block` (the
    # zero-inflated generator computes it over the fixed full N and restricts it
    # to the positive rows), use it verbatim, so a one-row replacement changes
    # exactly one block even when it toggles a row's positivity and changes |P|.
    # Otherwise partition by row position, which is data-independent only when
    # the row set is fixed (V4/V5 on the full n).
    if (is.null(block)) {
        set.seed(seed)
        v.block <- sample(rep_len(seq_len(m), n))
    } else {
        if (length(block) != n) stop("block must have length nrow(d) = ", n)
        v.block <- block
    }

    # per-block fit -> (beta^(b), sigma2^(b)).  non-private intermediates.
    l.beta <- Map(function(b) {
        mo.b <- fit.block(d[v.block == b, , drop = FALSE])
        v.out <- `names<-`(rep(NA_real_, q), v.qnames)
        if (!is.null(mo.b)) {
            cf <- stats::coef(mo.b)
            v.c <- intersect(v.names, names(cf))
            v.out[v.c] <- cf[v.c]
            # residual df can vanish in a thin block; NA then falls through to
            # the public box centre below, exactly as a failed fit does.
            if (disp) v.out[[".sigma2"]] <- tryCatch(stats::sigma(mo.b)^2,
                                                     error = function(e) NA_real_)
        }
        v.out
    }, seq_len(m))
    m.beta <- do.call(rbind, l.beta)                    # m x q

    # Widened Winsorized Mean: clip to the public box, impute failed blocks with
    # the (public) box centre, average coordinate-wise.
    m.clip <- mapply(function(j) {
        v <- pmin(pmax(m.beta[, j], v.lo[j]), v.hi[j])
        v[is.na(v) | !is.finite(v)] <- v.center[j]
        v
    }, seq_len(q))
    if (is.null(dim(m.clip))) m.clip <- matrix(m.clip, nrow = m)
    v.agg <- `names<-`(colMeans(m.clip), v.qnames)

    # Gaussian mechanism on the scaled aggregate.  The scaled aggregate has L2
    # sensitivity sqrt(q), so the release is (sqrt(q) / sigma.Q)-GDP.  Supplying
    # mu directly lets it compose in quadrature with the other channels; without
    # it, calibrate through the exact dual.
    #
    # This used to read sqrt(p) * .dp.gm.multiplier(eps, delta), which falls back
    # to the classical Dwork-Roth bound below eps = 1 and so over-noised there.
    # The exact dual is valid at every eps and is what the mu accounting needs.
    v.s     <- (v.hi - v.lo) / m
    sigma.Q <- if (!is.null(mu)) sqrt(q) / mu else .gdp.sigma(sqrt(q), eps, delta)
    set.seed(seed + 1L)
    v.rel   <- `names<-`(v.agg + v.s * sigma.Q * stats::rnorm(q), v.qnames)

    list(beta.hat = beta.hat,
         beta.dp  = v.rel[v.names],
         sigma2.dp = if (disp) max(v.rel[[".sigma2"]], 1e-6) else NA_real_,
         beta.agg = v.agg[v.names],
         sigma = v.s * sigma.Q, sigma.Q = sigma.Q, s.j = v.s, m = m,
         mu.total = sqrt(q) / sigma.Q,
         n.failed = sum(is.na(m.beta[, 1])), clip.lo = v.lo, clip.hi = v.hi,
         block = v.block)
}

# +------------------------------------------------------------------------+
# | .smooth.sens.beta.loglink -- analytic smooth sensitivity (D1, "try")   |
# +------------------------------------------------------------------------+
# Fisher-information analog of the OLS smooth sensitivity: X'WX in place of
# X'X and a public score-residual bound B.r in place of B.y.  Heuristic
# (assumes W and lambda_min stable across neighbours); carried for the
# smaller-when-assumptions-hold comparison against the worst-case route.

#' Analytic log-link smooth sensitivity (internal, comparison route)
#' @noRd
.smooth.sens.beta.loglink <- function(mo, eps, delta, B.x = NULL, B.r = NULL) {
    m.X   <- stats::model.matrix(mo)
    fam   <- stats::family(mo)
    v.mu  <- stats::fitted(mo)
    v.eta <- fam[["linkfun"]](v.mu)
    v.dmu.deta <- fam[["mu.eta"]](v.eta)
    v.V   <- fam[["variance"]](v.mu)
    v.w   <- v.dmu.deta^2 / v.V
    m.XtWX <- crossprod(sweep(m.X, 1, sqrt(v.w), "*"))
    LS <- 2 * max(sqrt(rowSums(stats::dfbeta(mo)^2)))
    if (is.null(B.x)) B.x <- sqrt(max(rowSums(m.X^2)))
    v.y <- stats::model.response(stats::model.frame(mo))
    v.g <- (v.dmu.deta / v.V) * (v.y - v.mu)
    if (is.null(B.r)) B.r <- max(abs(v.g))
    lambda.min <- tryCatch(
        min(eigen(m.XtWX, symmetric = TRUE, only.values = TRUE)[["values"]]),
        error = function(e) NA_real_)
    if (is.na(lambda.min) || lambda.min <= 0)
        return(list(S.smooth = NA_real_, S.local = LS, dstep = NA_real_,
                    ratio = NA_real_, lambda.min = lambda.min))
    dstep  <- 2 * B.x * B.r / lambda.min
    beta.p <- eps / (2 * log(2 / delta))
    k.star <- 1 / beta.p - LS / dstep
    S.smooth <- if (k.star > 0)
        (dstep / (beta.p * exp(1))) * exp(beta.p * LS / dstep) else LS
    list(S.smooth = S.smooth, S.local = LS, dstep = dstep,
         ratio = S.smooth / LS, lambda.min = lambda.min, B.x = B.x, B.r = B.r)
}

# +------------------------------------------------------------------------+
# | .smooth.sens.beta.loglink.std -- D1 made viable by public standardis.  |
# +------------------------------------------------------------------------+
# The raw analytic envelope is unusable on a design with a large-scale
# covariate (B.x ~ max||x_i|| large, lambda_min(X'WX) small).  Standardising
# by a public per-coefficient scale c_j (X_std = X diag(c)^{-1}) equalises the
# columns, so the Fisher information is well-conditioned and the envelope is
# finite.  The scales are public, so beta_std = diag(c) beta is a public
# reparametrisation and beta_DP = diag(c)^{-1} beta_std_DP is post-processing;
# c and a public score-residual bound B.r are the route's public inputs.

#' Public-standardised log-link smooth sensitivity (internal, D1)
#' @noRd
.smooth.sens.beta.loglink.std <- function(mo, eps, delta, x.scale,
                                          B.x = NULL, B.r = NULL) {
    v.beta <- stats::coef(mo)
    m.X    <- stats::model.matrix(mo)
    v.c    <- x.scale[names(v.beta)]
    if (any(is.na(v.c)) || any(v.c <= 0))
        stop("x.scale must give a positive public scale for every coefficient")
    m.X.std <- sweep(m.X, 2, v.c, "/")
    fam   <- stats::family(mo)
    v.mu  <- stats::fitted(mo)
    v.eta <- fam[["linkfun"]](v.mu)
    v.dmu.deta <- fam[["mu.eta"]](v.eta)
    v.V   <- fam[["variance"]](v.mu)
    v.w   <- v.dmu.deta^2 / v.V
    m.XtWX <- crossprod(sweep(m.X.std, 1, sqrt(v.w), "*"))
    m.dfbeta.std <- sweep(stats::dfbeta(mo), 2, v.c, "*")
    S.local <- 2 * max(sqrt(rowSums(m.dfbeta.std^2)))
    if (is.null(B.x)) B.x <- sqrt(max(rowSums(m.X.std^2)))
    v.y <- stats::model.response(stats::model.frame(mo))
    v.g <- (v.dmu.deta / v.V) * (v.y - v.mu)
    if (is.null(B.r)) B.r <- max(abs(v.g))
    lambda.min <- tryCatch(
        min(eigen(m.XtWX, symmetric = TRUE, only.values = TRUE)[["values"]]),
        error = function(e) NA_real_)
    if (is.na(lambda.min) || lambda.min <= 0)
        return(list(S.smooth.std = NA_real_, S.local.std = S.local,
                    dstep = NA_real_, ratio = NA_real_, lambda.min = lambda.min,
                    sigma.native = NA_real_))
    dstep  <- 2 * B.x * B.r / lambda.min
    beta.p <- eps / (2 * log(2 / delta))
    k.star <- 1 / beta.p - S.local / dstep
    S.smooth <- if (k.star > 0)
        (dstep / (beta.p * exp(1))) * exp(beta.p * S.local / dstep) else S.local
    s.gm <- .dp.gm.multiplier(eps, delta)
    v.sigma.native <- `names<-`((S.smooth * s.gm) / v.c, names(v.beta))
    list(S.smooth.std = S.smooth, S.local.std = S.local, dstep = dstep,
         ratio = S.smooth / S.local, lambda.min = lambda.min,
         sigma.native = v.sigma.native, B.x = B.x, B.r = B.r)
}

# +------------------------------------------------------------------------+
# | gen.syn.dp.projected.subagg -- public API                              |
# +------------------------------------------------------------------------+

#' Worst-case DP generator (subsample-and-aggregate)
#'
#' Generates synthetic data whose refitted GLM coefficient equals a worst-case
#' \eqn{(\varepsilon, \delta)}-DP estimate obtained by subsample-and-aggregate
#' (Nissim-Raskhodnikova-Smith 2007).  The rows are partitioned into \code{m}
#' data-independent blocks, the GLM is fit per block, each coordinate is
#' aggregated by a Widened Winsorized Mean over a public clip box, and Gaussian
#' noise is calibrated to the worst-case sensitivity
#' \eqn{(\mathrm{hi}_j - \mathrm{lo}_j)/m}.  The response is then \strong{drawn
#' from the model's own family} at the mean implied by
#' \eqn{\hat\beta_\mathrm{DP}} and at the declared dispersion.
#'
#' This is the default sound route.  Its \emph{sensitivity} comes from clipping
#' the per-block coefficients into a public box, so the noise calibration needs
#' \strong{neither} a public design \strong{nor} a bounded response: it is sound
#' under full-row replacement, and it holds however large or unbounded \eqn{y} is.
#' Where \code{\link{gen.syn.dp.ols.public}} is unavailable for want of a public \eqn{B_y},
#' this route lives.
#'
#' Its \emph{release}, however, publishes the covariates \strong{as they are}, so
#' the route does assume \eqn{X} is public.  That is a statement about the
#' output, not about the calibration, and the two must not be conflated.  Where
#' the covariates are themselves sensitive, use \code{\link{gen.syn.dp.full}} or
#' \code{\link{gen.syn.dp.logistic}}, which add a DP release of \eqn{X} and
#' compose it into the same mu quadrature.  Passing \eqn{X} through a copula
#' here would produce synthetic-looking covariates and so give the misleading
#' impression that \eqn{X} enjoyed some protection.  It did not.  Releasing
#' \eqn{X} plainly makes the assumption visible in the output instead of hiding
#' it.
#'
#' @section The release depends on the data only through beta_DP:
#' Before version 0.2.0 the response was built by \emph{virtual-response
#' injection}, \eqn{y_\mathrm{virt} = y + \mu(\hat\beta_\mathrm{DP}) -
#' \mu(\hat\beta)}, followed by a projection.  That construction \strong{retains
#' the real y}, so the release reproduced the real data's residuals, in scale
#' and in shape, and therefore did \emph{not} depend on the data only through
#' \eqn{\hat\beta_\mathrm{DP}}, which is the hypothesis the
#' post-processing clause of the subsample-and-aggregate proposition needs.  A
#' canary distinguishing game scored the resulting release at AUC 1.000 against
#' a mu-GDP bound of 0.566 at \eqn{\varepsilon = 1}: a perfect distinguisher, and
#' flat in \eqn{\varepsilon}, because the leak did not go through the noise.
#'
#' The response is now drawn from the family, so the release is a function of
#' \eqn{(\hat\beta_\mathrm{DP}, \mathrm{dispersion}, X, \mathrm{seed})} and
#' nothing else.  This costs the release its \emph{distributional} fidelity to
#' the empirical residuals and keeps its \emph{inferential} fidelity: for the
#' Gaussian family the analyst's refit still recovers \eqn{\hat\beta_\mathrm{DP}}
#' and \eqn{\hat\sigma^2_\mathrm{DP}} exactly.
#'
#' @param d Data frame with the response and predictors.
#' @param formula Model formula.
#' @param epsilon Worst-case privacy budget \eqn{\varepsilon > 0}.
#' @param delta Privacy parameter \eqn{\delta \in (0, 1)}.  Default 1e-6.
#' @param family A GLM family object (default \code{Gamma(link = "log")}).
#' @param m Block count.  Default \code{NULL} applies the public rule
#'   \code{\link{dp.subagg.blocks}}, \eqn{\lfloor n / (5 p) \rfloor}, the
#'   largest count leaving five rows per released coefficient per block.  An
#'   explicit value overrides.
#' @param clip.lo,clip.hi public per-coefficient clip bounds (named numeric
#'   vectors aligned to the coefficient names, or scalars recycled).
#'   Required: the one public input the worst-case guarantee needs.  The box
#'   \emph{width} sets the noise scale, so a box read off the sample makes the
#'   noise a function of the private data and voids the guarantee.
#' @param dispersion A \code{dp.dispersion} object, \strong{required for the
#'   families that have a free dispersion} (gaussian, gamma, quasi-gamma).  A
#'   GLM has two parameters and DP must cover both.  Use
#'   \code{\link{dp.dispersion.public}} when \eqn{\sigma^2} is external or
#'   already published (no budget spent, and \eqn{\hat\beta_\mathrm{DP}} stays
#'   exactly as noisy as before), or \code{\link{dp.dispersion.private}} to pay
#'   for it out of the same subagg channel (costs about 8\% more coefficient
#'   noise at \eqn{p = 6}).  Poisson and binomial have no free dispersion and
#'   must not be given one.
#' @param seed Optional integer RNG seed.
#' @return An object of class \code{synthexpln}: a list with \code{$syn} (the DP
#'   release), \code{$beta.dp}, \code{$sigma2.dp}, \code{$sigma}
#'   (per-coefficient noise sd), \code{$mu.total}, \code{$m}, \code{$epsilon},
#'   \code{$delta}, \code{$dispersion}, \code{$variant = "D2"} and
#'   \code{$scope = "analysis-level"}.
#'
#'   \code{$beta.hat}, \code{$beta.agg} and \code{$n.failed} are
#'   \strong{NON-PRIVATE intermediates} carried for evaluation only
#'   (\code{n.failed} is data-dependent).  They must not be treated as part of
#'   the release or composed downstream.
#' @seealso \code{\link{dp.dispersion.public}},
#'   \code{\link{dp.dispersion.private}}.
#' @export
#' @examples
#' set.seed(1)
#' n  <- 1000
#' x1 <- rnorm(n); x2 <- rnorm(n)
#' y  <- rgamma(n, 2, scale = exp(0.5 + 0.4 * x1 - 0.2 * x2) / 2)
#' d  <- data.frame(y = y, x1 = x1, x2 = x2)
#' # clip.lo/clip.hi are public per-coefficient bounds; a generous budget
#' # keeps the worst-case log-link refit numerically stable.
#' syn <- gen.syn.dp.projected.subagg(
#'          d, y ~ x1 + x2, epsilon = 10, family = Gamma(link = "log"),
#'          clip.lo = -1.5, clip.hi = 1.5,
#'          dispersion = dp.dispersion.public(
#'            0.5, source = "shape 2 gamma, so the dispersion is 1/2 by design"),
#'          seed = 2)
#' coef(glm(y ~ x1 + x2, syn$syn, family = Gamma(link = "log")))
gen.syn.dp.projected.subagg <- function(d, formula, epsilon, delta = 1e-6,
                                        family = stats::Gamma(link = "log"),
                                        m = NULL, clip.lo, clip.hi,
                                        dispersion = NULL,
                                        seed = NULL) {
    if (missing(clip.lo) || missing(clip.hi))
        stop("clip.lo and clip.hi (public per-coefficient bounds) are required")
    if (!is.null(seed)) set.seed(seed)
    seed.use <- if (is.null(seed)) 6489L else seed

    fam.nm.raw <- family[["family"]]
    fam.nm <- switch(fam.nm.raw,
        Gamma = "gamma", poisson = "poisson", tolower(fam.nm.raw))
    if (grepl("quasi", fam.nm)) fam.nm <- "quasi-gamma"
    link.nm <- family[["link"]]

    # ---- the dispersion must be declared --------------------------------------
    # A GLM has two parameters and DP must cover both.  Until 0.2.0 only the
    # coefficient passed through a DP channel and the dispersion rode in on the
    # real y, through the virtual-response construction below.
    b.needs <- .dp.has.dispersion(fam.nm)
    if (b.needs && is.null(dispersion))
        stop("dispersion is required for the ", fam.nm, " family: it has a free ",
             "dispersion parameter, and DP must cover it.\n",
             "  Supply dp.dispersion.public(sigma2, source) if sigma^2 is\n",
             "  external or already published (this spends no budget), or\n",
             "  dp.dispersion.private(lo, hi, source) to spend budget on it.",
             call. = FALSE)
    if (!is.null(dispersion) && !inherits(dispersion, "dp.dispersion"))
        stop("dispersion must be built by dp.dispersion.public() or ",
             "dp.dispersion.private()", call. = FALSE)
    if (!b.needs && !is.null(dispersion))
        stop("the ", fam.nm, " family has no free dispersion (its variance is ",
             "fixed by the mean), so do not pass one", call. = FALSE)

    # A private dispersion rides along as a (p+1)-th aggregated coordinate, so it
    # is paid for out of the same subagg channel and the L2 sensitivity becomes
    # sqrt(p + 1).  A public one costs nothing and is reused.
    b.disp.paid <- !is.null(dispersion) &&
        identical(dispersion[["mode"]], "private")

    l.sa <- .subsample.aggregate.beta.loglink(
        d, formula, family, fam.nm, link.nm, epsilon, delta, m = m,
        clip.lo = clip.lo, clip.hi = clip.hi, seed = seed.use,
        disp    = b.disp.paid,
        disp.lo = if (b.disp.paid) dispersion[["lo"]] else NULL,
        disp.hi = if (b.disp.paid) dispersion[["hi"]] else NULL)

    beta.hat <- l.sa[["beta.hat"]]
    beta.dp  <- l.sa[["beta.dp"]]

    sigma2.use <- if (!b.needs) NA_real_
                  else if (b.disp.paid) l.sa[["sigma2.dp"]]
                  else dispersion[["sigma2"]]

    # ---- build the release from the DP output only ---------------------------
    # This replaces .inject.virtual.y + projection(), which set
    # y_virt = y + mu(beta_DP) - mu(beta_hat) and so RETAINED the real y.  See
    # the "release depends on the data only through beta_DP" section above.
    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)
    d.syn <- .synth.response(d, y.nm, x.nms, fam.nm, link.nm,
                             beta.dp, sigma2.use, seed.use)

    structure(list(
        syn        = d.syn,
        beta.hat   = beta.hat,        # NOT PRIVATE: evaluation only
        beta.dp    = beta.dp,
        beta.agg   = l.sa[["beta.agg"]],   # NOT PRIVATE
        sigma2.dp  = sigma2.use,
        sigma      = l.sa[["sigma"]],
        # subagg's noise is already per-coefficient and on the raw scale
        # (s_j * sigma.Q), so this is it.  Used by calibrated.se().
        sigma.beta.raw = l.sa[["sigma"]][names(beta.dp)],
        sigma.Q    = l.sa[["sigma.Q"]],
        s.j        = l.sa[["s.j"]],
        mu.total   = l.sa[["mu.total"]],
        m          = l.sa[["m"]],
        n.failed   = l.sa[["n.failed"]],   # NOT PRIVATE (data-dependent)
        epsilon    = epsilon,
        delta      = delta,
        dispersion = dispersion,
        variant    = "D2",
        scope      = "analysis-level"
    ), class = "synthexpln")
}
