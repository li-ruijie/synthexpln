# dp-logistic.R -- V5 binary-outcome DP via Bernoulli sampling.
# Public function: gen.syn.dp.logistic()
# Ported from lib/53-dp-logistic.r lines 45-106.
#
# Note: supports both public X (x.public = TRUE) and sensitive X via A's DP
# copula marginals (x.public = FALSE, requires x.spec).  See sub-project B.

#' V5 binary-outcome DP synthetic data generator (Bernoulli sampling)
#'
#' Generates synthetic binary-response data.  The coefficient is released by
#' subsample-and-aggregate over a public clip box, and the response is then drawn
#' as
#' \eqn{Y_i \sim \mathrm{Bernoulli}(\mathrm{expit}(x_i^\top \hat\beta_\mathrm{DP}))},
#' which is post-processing of the DP output and costs no budget.  Covariates may
#' be public (\code{x.public = TRUE}) or released under DP from \code{x.spec}.
#'
#' @details
#' The algorithm proceeds in four steps.
#' \enumerate{
#'   \item \strong{Logistic regression.}  Fit \code{glm(formula, family =
#'     binomial)} on the original data.  This is a non-private reference, kept
#'     for evaluation only.
#'   \item \strong{DP beta by subsample-and-aggregate.}  Rows are partitioned
#'     into \code{m} data-independent blocks, the logistic model is fit per
#'     block, and each coordinate is aggregated by a Widened Winsorized Mean over
#'     a public box.  A one-row change perturbs exactly one block, whose clipped
#'     coefficient moves by at most \eqn{\mathrm{hi}_j - \mathrm{lo}_j}, so the
#'     mean moves by at most \eqn{s_j = (\mathrm{hi}_j - \mathrm{lo}_j)/m}.  The
#'     binomial has \strong{no free dispersion} (its variance is fixed by the
#'     mean), so the aggregate stays \eqn{p}-dimensional and there is no
#'     dispersion channel to leak.
#'   \item \strong{X_syn.}  With \code{x.public = TRUE} the original covariate
#'     matrix is reused.  With \code{x.public = FALSE} the covariates are
#'     released under DP from \code{x.spec} (independent DP marginals by default,
#'     or a DP Gaussian copula when \code{use.copula = TRUE}), so
#'     \eqn{(X_\mathrm{syn}, Y_\mathrm{syn})} is jointly DP and every channel
#'     composes in the same mu quadrature.
#'   \item \strong{Bernoulli sampling.}  Draw
#'     \eqn{Y_i \sim \mathrm{Bernoulli}(\mathrm{expit}(x_i^\top
#'     \hat\beta_\mathrm{DP}))}.  A function of the DP output and public
#'     information alone, so it is post-processing and spends nothing.
#' }
#'
#' @section What changed, and why the old mechanism had no theorem:
#' Before version 0.2.0 the coefficient channel was calibrated to the DFBETA
#' \strong{local} sensitivity, on a column-standardised scale and with a 3-sigma
#' truncation.  Local sensitivity is a function of the realised data, so the
#' release did not satisfy \eqn{(\varepsilon, \delta)}-DP at the stated
#' \eqn{\varepsilon} or at any other.  The Bernoulli draw was always sound
#' post-processing; it was the coefficient channel feeding it, and only that,
#' which was at fault.
#'
#' The legacy path remains reachable as \code{sens.method = "local"} so that
#' pre-0.2.0 runs can be reproduced.  It warns on every call.
#'
#' @param d Data frame containing the response and predictors.
#' @param formula Model formula with a binary (0/1 or TRUE/FALSE) response
#'   (e.g. \code{y ~ x1 + x2}).
#' @param epsilon Total privacy budget \eqn{\varepsilon > 0}.
#' @param delta Privacy parameter \eqn{\delta \in (0, 1)}.  Default
#'   \code{1e-6}.
#' @param x.public Logical flag. If \code{TRUE} (default) the covariate matrix
#'   is public and reused unchanged.  If \code{FALSE}, covariates are privatised
#'   from \code{x.spec} (independent DP marginals, or the copula via
#'   \code{use.copula}) and the budget splits per \code{eps.split}.
#' @param x.spec Named list (one entry per predictor) of \code{list(type =
#'   "continuous", breaks = ...)} or \code{list(type = "categorical", levels =
#'   ...)}.  Required when \code{x.public = FALSE}.
#' @param sens.method \code{"subagg"} (default, sound) releases the coefficient
#'   by subsample-and-aggregate over a public box.  \code{"local"} is the retired
#'   DFBETA local-sensitivity path, kept only to reproduce pre-0.2.0 runs; it
#'   carries \strong{no formal DP guarantee} and warns on every call.
#' @param clip.lo,clip.hi public per-coefficient bounds (named numeric vectors
#'   aligned to the coefficient names, or scalars recycled).  Required for
#'   \code{sens.method = "subagg"}.  The box \emph{width} sets the noise scale,
#'   so a box read off the sample makes the noise a function of the private data
#'   and voids the guarantee.
#' @param m Subsample-and-aggregate block count.  Default \code{NULL} applies
#'   the public rule \code{\link{dp.subagg.blocks}},
#'   \eqn{\lfloor n / (5 p) \rfloor}.  An explicit value overrides.
#' @param mu.split Named numeric \code{c(marg, corr, coef)} summing to 1, the
#'   fractions of the \eqn{\mu^2} budget going to the covariate marginals, the
#'   copula correlation, and the coefficient release.  Used only when
#'   \code{x.public = FALSE}; with public covariates the whole budget goes to the
#'   coefficient.  Default \code{c(marg = 0.05, corr = 0.05, coef = 0.90)}.
#' @param eps.split Legacy epsilon allocation \code{c(X, beta)} summing to 1.
#'   Used only by \code{sens.method = "local"} with \code{x.public = FALSE}.
#' @param eps.corr.frac Fraction of \code{eps.X} spent on the copula
#'   correlation release (the rest goes to the marginals).  Default 0.5.
#' @param use.copula Logical flag to model the joint covariate dependence with
#'   the DP Gaussian copula (requires >= 2 continuous covariates).  Default
#'   \code{FALSE} (independent DP marginals).  Set \code{TRUE} for descriptive
#'   X-fidelity, at some inferential cost.  Since 0.7.0 the copula is
#'   mixed-margin, exactly as in \code{\link{gen.syn.dp.full}}: the categorical
#'   covariates join it through bounded indicator scores and the
#'   \code{back.map}, under assumptions A1 to A6 documented there.
#' @param score.method Copula score transform when \code{use.copula = TRUE}:
#'   \code{"adaptive"} (default, PIT through each covariate's released DP
#'   histogram, best for skewed marginals), \code{"normal"} (public normal-score
#'   PIT, recovers correlation at \eqn{\varepsilon \ge 5} on location-scale
#'   margins), or \code{"box"} (legacy, noise-dominated).
#' @param back.map Back-map from a categorical's copula latents to a level
#'   under \code{use.copula = TRUE}: \code{"probit"} (default, recommended) or
#'   \code{"logit"}, which attenuates the recovered association.  See the
#'   \emph{Mixed-margin copula} section of \code{\link{gen.syn.dp.full}}.
#' @param seed Optional integer RNG seed.  If \code{NULL} the global RNG
#'   state is used.
#' @return An object of class \code{synthexpln}: a named list with
#'   \describe{
#'     \item{\code{$syn}}{Synthetic data frame with DP binary response.  This is
#'       the DP release.  Everything below is evaluation-only.}
#'     \item{\code{$beta.hat}}{Original coefficient vector (NOT private).}
#'     \item{\code{$beta.dp}}{DP coefficient vector.}
#'     \item{\code{$sigma.beta}}{Per-coefficient noise sd.}
#'     \item{\code{$sens.method}}{Which calibration ran.}
#'     \item{\code{$epsilon}, \code{$delta}}{As supplied.}
#'     \item{\code{$mu.total}}{Total mu-GDP level (subagg only).}
#'     \item{\code{$mu.channels}}{Per-channel mu \code{c(marg, corr, coef)}.}
#'     \item{\code{$mu.check}}{The channels re-composed in quadrature.  Equals
#'       \code{$mu.total}.}
#'     \item{\code{$subagg}}{The subagg intermediates.  \code{$beta.agg} and
#'       \code{$n.failed} inside it are \strong{NON-PRIVATE}.}
#'     \item{\code{$variant}}{\code{"V5"}.}
#'     \item{\code{$scope}}{\code{"joint-release"}.}
#'   }
#' @export
#' @examples
#' set.seed(42)
#' n <- 300
#' x1 <- rnorm(n)
#' x2 <- factor(rbinom(n, 1, 0.5))
#' y  <- rbinom(n, 1, plogis(-0.5 + 0.8 * x1))
#' d  <- data.frame(y = y, x1 = x1, x2 = x2)
#' # clip.lo/clip.hi are a public box on the log-odds coefficients: effects of
#' # this size are what the design anticipates, fixed before the data are seen.
#' syn <- gen.syn.dp.logistic(d, y ~ x1 + x2, epsilon = 5, delta = 1e-6,
#'                             clip.lo = -2, clip.hi = 2, seed = 1)
#' table(syn$syn$y)
gen.syn.dp.logistic <- function(d, formula, epsilon, delta = 1e-6,
                                 x.public = TRUE, x.spec = NULL,
                                 sens.method = c("subagg", "local"),
                                 clip.lo = NULL, clip.hi = NULL, m = NULL,
                                 mu.split = c(marg = 0.05, corr = 0.05,
                                              coef = 0.90),
                                 eps.split = c(X = 0.2, beta = 0.8),
                                 eps.corr.frac = 0.5, use.copula = FALSE,
                                 score.method = "adaptive",
                                 back.map = c("probit", "logit"), seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    back.map <- match.arg(back.map)
    sens.method  <- match.arg(sens.method)
    score.method <- match.arg(score.method, c("box", "normal", "adaptive"))
    if (!x.public && is.null(x.spec))
        stop("x.public = FALSE requires x.spec (per-column type/breaks/levels)")

    y.nm <- all.vars(formula)[1]
    if (!all(d[[y.nm]] %in% c(0, 1, TRUE, FALSE)))
        stop("y must be binary (0/1 or TRUE/FALSE)")
    d[[y.nm]] <- as.integer(d[[y.nm]])
    x.nms <- all.vars(formula)[-1]
    # Sensitive-X: code categorical predictors as factors before the fit so
    # beta_hat's dummy names match the sampled X_syn model matrix.
    if (!x.public)
        d <- Reduce(function(acc, nm) {
            spec <- x.spec[[nm]]
            if (spec[["type"]] == "categorical")
                acc[[nm]] <- factor(as.character(acc[[nm]]), levels = spec[["levels"]])
            acc
        }, x.nms, init = d)

    # copula eligibility, hoisted because the budget split needs it.
    v.cont.nms <- if (x.public) character(0)
                  else x.nms[vapply(x.nms, function(nm)
                      x.spec[[nm]][["type"]] == "continuous", logical(1))]
    do.copula  <- !x.public && isTRUE(use.copula) && length(v.cont.nms) >= 2L
    n          <- nrow(d)

    # ---- allocate the budget -------------------------------------------------
    # subagg: every channel is a Gaussian mechanism, so the release composes in
    # one mu-GDP quadrature.  With public X there is no X channel and the whole
    # mu goes to the coefficient.  The Bernoulli draw from expit(x' beta_DP) is
    # post-processing and costs nothing.
    if (sens.method == "subagg") {
        if (is.null(clip.lo) || is.null(clip.hi))
            stop("sens.method = \"subagg\" requires the public per-coefficient ",
                 "bounds clip.lo and clip.hi", call. = FALSE)
        if (abs(sum(mu.split) - 1) > 1e-8) stop("mu.split must sum to 1")
        mu.total <- .gdp.mu.from.eps.delta(epsilon, delta)
        mu.coef  <- if (x.public) mu.total
                    else mu.total * sqrt(mu.split[["coef"]])
        mu.corr  <- if (do.copula) mu.total * sqrt(mu.split[["corr"]]) else 0
        mu.marg  <- if (x.public) 0
                    else if (do.copula) mu.total * sqrt(mu.split[["marg"]])
                    else .gdp.compose(mu.total *
                             sqrt(c(mu.split[["marg"]], mu.split[["corr"]])))
        mu.each  <- if (x.public) 0 else mu.marg / sqrt(length(x.nms))
        eps.X <- NA_real_; eps.beta <- NA_real_
        eps.corr <- NA_real_; eps.marg <- NA_real_; eps.each <- NULL
    } else {
        warning("sens.method = \"local\" calibrates beta to the DFBETA LOCAL ",
                "sensitivity, which is a function of the realised data. The ",
                "release therefore carries NO formal (eps, delta)-DP guarantee ",
                "(retired 2026-07-07). Retained only to reproduce pre-0.2.0 ",
                "runs. Use sens.method = \"subagg\".", call. = FALSE)
        if (abs(sum(eps.split) - 1) > 1e-8) stop("eps.split must sum to 1")
        # Public X spends the whole budget on beta; sensitive X splits it.
        if (x.public) { eps.X <- 0; eps.beta <- epsilon }
        else { eps.X <- epsilon * eps.split[["X"]]
               eps.beta <- epsilon * eps.split[["beta"]] }
        eps.corr <- if (do.copula) eps.X * eps.corr.frac else 0
        eps.marg <- eps.X - eps.corr
        eps.each <- eps.marg / length(x.nms)
        mu.total <- NA_real_; mu.coef <- NA_real_
        mu.corr  <- NA_real_; mu.marg <- NA_real_; mu.each <- NULL
    }

    dp.marg.cont <- function(x, breaks)
        if (sens.method == "subagg") .dp.hist.continuous(x, breaks, mu = mu.each)
        else .dp.hist.continuous(x, breaks, eps = eps.each)
    dp.marg.cat <- function(x, levels)
        if (sens.method == "subagg") .dp.cat(x, levels, mu = mu.each)
        else .dp.cat(x, levels, eps = eps.each)
    dp.marg.cat.probs <- function(x, levels)
        if (sens.method == "subagg") .dp.cat.probs(x, levels, mu = mu.each)
        else .dp.cat.probs(x, levels, eps = eps.each)

    # Step 1: fit logistic regression (non-private reference, evaluation only)
    mo       <- stats::glm(formula, d, family = stats::binomial(link = "logit"))
    beta.hat <- stats::coef(mo)
    m.X      <- stats::model.matrix(formula, d)
    p        <- ncol(m.X)

    # Step 2: DP beta
    if (sens.method == "subagg") {
        # Subsample-and-aggregate, binomial/logit.  Worst-case per-coordinate
        # sensitivity s_j = (hi_j - lo_j)/m over a data-independent partition, so
        # the release is mu_coef-GDP and is sound under full-row replacement,
        # which is the relation that matters once X is sensitive.  The binomial
        # has no free dispersion (its variance is fixed by the mean), so the
        # aggregate stays p-dimensional and there is no dispersion channel to
        # leak.
        l.sa <- .subsample.aggregate.beta.loglink(
            d, formula, stats::binomial(link = "logit"), "binomial", "logit",
            eps = NA_real_, delta = delta, m = m,
            clip.lo = clip.lo, clip.hi = clip.hi,
            seed = if (is.null(seed)) 6489L else seed + 10L,
            mu = mu.coef)
        beta.dp <- l.sa[["beta.dp"]]
        sigma   <- l.sa[["sigma"]]
        v.s     <- NULL
    } else {
        l.sa <- NULL
        # retired: DFBETA local sensitivity.  Data-dependent, so no formal
        # (eps, delta)-DP.  Kept only to reproduce pre-0.2.0 runs.
        v.sd <- apply(m.X, 2, stats::sd)
        v.s  <- ifelse(v.sd > 1e-8, v.sd, 1)
        m.dfb    <- stats::dfbeta(mo)
        LS.tilde <- 2 * max(sqrt(rowSums(sweep(m.dfb, 2, v.s, "*")^2)))
        sigma    <- LS.tilde * .dp.gm.multiplier(eps.beta, delta)

        eta.tilde <- stats::rnorm(p, 0, sigma)
        # 3-sigma L2 truncation (post-processing)
        s.norm <- sqrt(sum(eta.tilde^2))
        s.max  <- 3 * sqrt(p) * sigma
        if (s.norm > s.max) eta.tilde <- eta.tilde * (s.max / s.norm)

        beta.tilde <- beta.hat * v.s
        beta.dp    <- `names<-`((beta.tilde + eta.tilde) / v.s, names(beta.hat))
    }

    # Step 3: construct X_syn
    if (x.public) {
        d.syn      <- d
        m.X.syn    <- m.X
        m.R.dp     <- NULL
    } else {
        if (do.copula) {
            # Mixed-margin copula (0.7.0): the categoricals join the copula
            # through their K-1 indicator scores, exactly as gen.syn.dp.full's
            # Stage 1.  The released vech(S) covers p_total columns at the same
            # M^2 p / n sensitivity form, and sampling runs through the
            # calibrated back-map with the biserial de-attenuation (A6).
            # Assumptions A1-A6 on gen.syn.dp.full apply verbatim.
            v.cat.nms <- setdiff(x.nms, v.cont.nms)
            l.cont.marg <- Map(function(nm) {
                breaks <- x.spec[[nm]][["breaks"]]
                list(counts = dp.marg.cont(d[[nm]], breaks), breaks = breaks)
            }, v.cont.nms)
            l.cat.marg <- Map(function(nm)
                list(levels = x.spec[[nm]][["levels"]],
                     probs  = dp.marg.cat.probs(d[[nm]],
                                                x.spec[[nm]][["levels"]])),
                v.cat.nms)
            s.tau <- if (score.method == "box") 1 / (2 * n) else 0.005
            s.M   <- stats::qnorm(1 - s.tau)
            m.W.cont <- vapply(v.cont.nms, function(nm) {
                v.bx <- range(x.spec[[nm]][["breaks"]])
                .dp.copula.scores(d[[nm]], v.bx[1], v.bx[2], s.tau,
                                  method = score.method,
                                  counts = l.cont.marg[[nm]][["counts"]],
                                  breaks = l.cont.marg[[nm]][["breaks"]])
            }, numeric(n))
            m.W <- cbind(m.W.cont,
                         do.call(cbind, Map(function(nm) .cat.indicator.scores(
                             d[[nm]], l.cat.marg[[nm]][["levels"]],
                             l.cat.marg[[nm]][["probs"]], s.M), v.cat.nms)))
            m.R.dp <- if (sens.method == "subagg")
                          .dp.copula.corr(m.W, NULL, delta, s.M, mu = mu.corr)
                      else .dp.copula.corr(m.W, eps.corr, delta, s.M)
            d.syn <- .sample.from.dp.mixed.copula(m.R.dp,
                         c(l.cont.marg, l.cat.marg), v.cont.nms, v.cat.nms,
                         back.map, n)
            d.syn <- d.syn[x.nms]
        } else {
            m.R.dp <- NULL
            d.syn  <- Reduce(function(acc, nm) {
                spec <- x.spec[[nm]]
                if (spec[["type"]] == "continuous") {
                    counts    <- dp.marg.cont(d[[nm]], spec[["breaks"]])
                    acc[[nm]] <- .sample.from.dp.hist(counts, spec[["breaks"]], n)
                } else {
                    acc[[nm]] <- factor(dp.marg.cat(d[[nm]], spec[["levels"]]),
                                        levels = spec[["levels"]])
                }
                acc
            }, x.nms, init = data.frame(row.names = seq_len(n)))
        }
        m.X.syn <- stats::model.matrix(stats::reformulate(x.nms), d.syn)
    }

    # Align beta.dp to columns of m.X.syn (handles missing factor levels)
    v.common     <- intersect(names(beta.dp), colnames(m.X.syn))
    beta.aligned <- `names<-`(rep(0, ncol(m.X.syn)), colnames(m.X.syn))
    beta.aligned[v.common] <- beta.dp[v.common]

    # Step 4: Bernoulli sampling from DP linear predictor
    lin.pred      <- as.numeric(m.X.syn %*% beta.aligned)
    p.syn         <- stats::plogis(lin.pred)
    d.syn[[y.nm]] <- stats::rbinom(nrow(m.X.syn), 1, p.syn)

    # NOT private, evaluation only: beta.hat is the full-data fit, and subagg
    # carries its own non-private intermediates.  The DP release is $syn.
    structure(list(
        syn         = d.syn,
        beta.hat    = beta.hat,       # NOT PRIVATE
        beta.dp     = beta.dp,
        sigma.beta  = sigma,
        # Raw-scale noise sd per coefficient, subset by name.  See
        # gen.syn.dp.full.  (The binomial has no dispersion, so there is no
        # (p+1)-th coordinate here, but index by name anyway: positional renaming
        # is the class of defect that stays correct until it silently does not.)
        sigma.beta.raw = if (sens.method == "subagg")
                             sigma[names(beta.dp)]
                         else `names<-`(as.numeric(sigma) / v.s,
                                        names(beta.hat)),
        sens.method = sens.method,
        epsilon     = epsilon,
        delta       = delta,
        x.public    = x.public,
        use.copula  = do.copula,
        dp.corr     = m.R.dp,
        # mu-GDP accounting (subagg).  mu.check re-composes the channels.
        mu.total    = mu.total,
        mu.channels = c(marg = mu.marg, corr = mu.corr, coef = mu.coef),
        mu.check    = if (sens.method == "subagg")
                          .gdp.compose(c(mu.marg, mu.corr, mu.coef))
                      else NA_real_,
        subagg      = l.sa,
        # legacy (sens.method = "local") diagnostics
        s.j         = v.s,
        eps.allocation = if (sens.method == "local")
                             c(X = eps.X, beta = eps.beta, marg = eps.marg,
                               corr = eps.corr)
                         else NULL,
        delta.total = if (sens.method == "subagg") delta
                      else if (do.copula) 2 * delta else delta,
        variant     = "V5",
        scope       = "joint-release"
    ), class = "synthexpln")
}
