# dp-full.R -- V4 sensitive-X full DP synthesiser.
# Public function: gen.syn.dp.full()
# Ported from lib/35-dp-full.r lines 107-209.

#' V4 full-DP OLS synthetic data generator (sensitive covariates)
#'
#' Generates synthetic data when the covariate matrix \eqn{X} is also sensitive.
#' Every channel is a Gaussian mechanism, so the whole release composes in one
#' mu-GDP quadrature,
#' \eqn{\mu^2 = \mu_\mathrm{marg}^2 + \mu_\mathrm{corr}^2 + \mu_\mathrm{coef}^2},
#' converted to \eqn{(\varepsilon, \delta)} by the exact Gaussian-DP dual.
#'
#' @details
#' The algorithm proceeds in four stages.
#' \enumerate{
#'   \item \strong{DP X marginals.}  For each covariate named in \code{x.spec},
#'     a DP histogram (continuous) or DP categorical release (categorical).  A
#'     replace-one row swap decreases one count and increases another, so
#'     the count vector's L2 sensitivity is \eqn{\sqrt 2} and the noise sd is
#'     \eqn{\sqrt 2 / \mu}.  The \eqn{k} marginals are \eqn{k} functions of the
#'     same row, so they compose in quadrature and each gets
#'     \eqn{\mu_\mathrm{marg} / \sqrt k}.
#'   \item \strong{DP covariate dependence} (optional).  A Gaussian copula over
#'     the whole covariate block.  Each continuous column enters by a bounded
#'     normal-score PIT; since 0.6.0 each categorical column joins the same
#'     copula by its \eqn{K-1} bounded indicator scores computed from its
#'     released DP proportions (the mixed-margin copula), so the
#'     categorical-to-continuous joint structure is preserved rather than
#'     discarded.  The second-moment matrix of the extended score matrix is
#'     released at \eqn{\Delta_2 = M^2 p_\mathrm{total} / n} (the L2
#'     sensitivity of the released \eqn{vech(S)}, with \eqn{p_\mathrm{total}}
#'     counting the indicator columns), then projected to the nearest PSD
#'     correlation matrix as post-processing.  Categorical columns are
#'     reconstructed from their latents through a calibrated back-map (see
#'     \code{back.map}); continuous columns through their DP-histogram
#'     quantiles.  See the \emph{Mixed-margin copula} section for the
#'     assumptions this carries.
#'   \item \strong{DP beta and sigma\eqn{^2}, jointly, by
#'     subsample-and-aggregate.}  Rows are partitioned into \code{m}
#'     data-independent blocks, OLS is fit per block, and each coordinate,
#'     including \eqn{\sigma^2} as a \eqn{(p+1)}-th coordinate, is aggregated by
#'     a Widened Winsorized Mean over a public box.  One row changes exactly one
#'     block, and that block's clipped coefficient moves by at most
#'     \eqn{\mathrm{hi}_j - \mathrm{lo}_j}, so the mean moves by at most
#'     \eqn{s_j = (\mathrm{hi}_j - \mathrm{lo}_j)/m}.  The scaled aggregate has
#'     L2 sensitivity \eqn{\sqrt{p+1}}.
#'
#'     This is sound under \strong{full-row} replacement, which is the relation
#'     that matters once \eqn{X} is sensitive: no design constant and no
#'     data-dependent quantity enters, so it holds whether the neighbour changed
#'     \eqn{x}, \eqn{y}, or both.
#'   \item \strong{Synthetic y.}  Constructed as
#'     \eqn{X_\mathrm{syn}\hat\beta_\mathrm{DP} + e_\mathrm{proj}}, with
#'     \eqn{e_\mathrm{proj}} freshly drawn, projected onto the null space of
#'     \eqn{X_\mathrm{syn}} and rescaled so that
#'     \eqn{\|e_\mathrm{proj}\|^2 = (n-p)\hat\sigma^2_\mathrm{DP}}.  The real
#'     \eqn{y} appears nowhere: this stage is post-processing of the DP output.
#' }
#'
#' @section What changed, and why the old mechanism had no theorem:
#' Before version 0.2.0 the coefficient channel was calibrated to the DFBETA
#' \strong{local} sensitivity, which is a function of the realised data, so the
#' release did not satisfy \eqn{(\varepsilon, \delta)}-DP at the stated
#' \eqn{\varepsilon} or at any other: neighbouring datasets received different
#' noise scales, no bound held across the pair, and there was no theorem to
#' instantiate.  Supplying a public \code{B.y} was necessary but \strong{not
#' sufficient}, since \code{B.y} reached only the variance channel.
#'
#' The variance channel was itself scoped wrongly.  It released the residual norm
#' under \emph{response-row} replacement holding the design fixed, which is not
#' the relation that applies when the design is the object being protected.  It
#' was also Laplace, hence pure eps-DP, and so could not enter the mu-GDP
#' quadrature the covariate channels needed.  Subsample-and-aggregate replaces
#' both channels at once and needs no \code{B.y} at all.
#'
#' The legacy path remains reachable as \code{sens.method = "local"} so that
#' pre-0.2.0 runs can be reproduced.  It warns on every call.
#'
#' @section Mixed-margin copula, assumptions in application (A1 to A5):
#' With \code{use.copula = TRUE} the categorical covariates join the same
#' Gaussian copula as the continuous block, each encoded by \eqn{K-1} bounded
#' normal-scored level indicators computed from its released DP proportions
#' (public, so no extra privacy cost).  This is DESCRIPTIVE fidelity only: the
#' analyst's refit recovers \eqn{\hat\beta_\mathrm{DP}} exactly whatever joint
#' structure \eqn{X_\mathrm{syn}} carries, so inferential fidelity is
#' untouched either way.  The recovered joint carries the following
#' assumptions, documented identically in the paper's supplement (S6) and in
#' the research implementation \code{35-dp-full.r}:
#' \describe{
#'   \item{A1}{Gaussian-copula fit: only the pairwise, latent-monotone
#'     association is recovered; higher-order or non-elliptical dependence is
#'     not.}
#'   \item{A2}{Back-map: \code{"probit"} (argmax of the latents) is
#'     deterministic and transfers the full released correlation.
#'     \code{"logit"} (softmax, then one categorical draw per row) injects
#'     independent choice noise, so it \strong{attenuates the very
#'     association the mixed copula exists to preserve}; use it only when a
#'     stochastic assignment is specifically wanted, and expect the recovered
#'     association to land well below the probit one at the same budget.}
#'   \item{A3}{The extended sensitivity \eqn{M^2 p_\mathrm{total} / n} is a
#'     sound upper bound, attained only in the continuous sub-block: a
#'     categorical indicator column takes one of K fixed winsorised values
#'     strictly inside the clamp, so its channel is over-noised in the safe
#'     direction (delivered mu no larger than claimed).}
#'   \item{A4}{\eqn{p_\mathrm{total}} grows by \eqn{K-1} per nominal
#'     covariate, raising the correlation-channel noise.  Acceptable at the
#'     default split, where that channel is budget-starved and the
#'     coefficient channel dominates.}
#'   \item{A5}{Reference-level encoding: recovered associations are relative
#'     to the first declared level of each categorical.}
#'   \item{A6}{De-attenuation (0.7.0): the categorical-continuous cross
#'     entries of the released correlation are divided by the public biserial
#'     factor \eqn{c_k = \phi(\Phi^{-1}(p_k)) / \sqrt{p_k (1 - p_k)}} before
#'     sampling, correcting the point-biserial fit to the latent scale the
#'     back-map needs.  \eqn{c_k} reads the released proportions only, so the
#'     step is public post-processing and the release, its sensitivity, and
#'     the privacy accounting are untouched.  Cross-categorical blocks are
#'     not corrected.}
#' }
#'
#' @section x.spec format:
#' A named list with one entry per predictor (names must match the RHS of
#' \code{formula}).  Each entry is a list with at least:
#' \describe{
#'   \item{\code{type}}{Either \code{"continuous"} or \code{"categorical"}.}
#'   \item{\code{breaks}}{(continuous only) Numeric break-point vector of
#'     length \eqn{B+1} for a \eqn{B}-bin histogram.  Alternative: supply
#'     \code{bounds = c(lo, hi)} to derive 10 equal-width bins automatically
#'     via \code{seq(lo, hi, length.out = 11)}.}
#'   \item{\code{levels}}{(categorical only) Character vector of all possible
#'     levels.}
#' }
#'
#' @param d Data frame containing response and predictors.
#' @param formula Model formula (e.g. \code{y ~ age + sex}).
#' @param x.spec Named list specifying the type and public parameters of each
#'   predictor.  See the \emph{x.spec format} section.
#' @param epsilon Total privacy budget \eqn{\varepsilon > 0}.
#' @param delta Privacy parameter \eqn{\delta \in (0, 1)}.  Default
#'   \code{1e-6}.
#' @param sens.method \code{"subagg"} (default, sound) calibrates the
#'   coefficient and dispersion channels by subsample-and-aggregate over a
#'   public box.  \code{"local"} is the retired DFBETA local-sensitivity path,
#'   kept only to reproduce pre-0.2.0 runs; it carries \strong{no formal DP
#'   guarantee} and warns on every call.
#' @param clip.lo,clip.hi public per-coefficient bounds (named numeric vectors
#'   aligned to the coefficient names, or scalars recycled).  Required for
#'   \code{sens.method = "subagg"}.  The box \emph{width} sets the noise scale,
#'   so a box read off the sample makes the noise a function of the private data
#'   and voids the guarantee.  This is the unavoidable cost of a sensitive
#'   design, the old API having needed no bounds at all, which is how the
#'   defect went unnoticed.
#' @param disp.lo,disp.hi public bounds on \eqn{\sigma^2}.  Required for
#'   \code{sens.method = "subagg"}, where the dispersion rides along as a
#'   \eqn{(p+1)}-th aggregated coordinate.
#' @param m Subsample-and-aggregate block count.  Default \code{NULL} applies
#'   the public rule \code{\link{dp.subagg.blocks}},
#'   \eqn{\lfloor n / (5 p) \rfloor}.  An explicit value overrides.
#' @param mu.split Named numeric \code{c(marg, corr, coef)} summing to 1, the
#'   fractions of the \eqn{\mu^2} budget going to the covariate marginals, the
#'   copula correlation, and the coefficient release.  Default
#'   \code{c(marg = 0.05, corr = 0.05, coef = 0.90)}: the analyst's refit
#'   recovers \eqn{\hat\beta_\mathrm{DP}} exactly, so the joint structure of
#'   \eqn{X_\mathrm{syn}} does not touch inferential fidelity and the marginal
#'   channels are starved in its favour.  Shift it towards \code{marg}/
#'   \code{corr} for releases that also need descriptive fidelity of \eqn{X}.
#' @param eps.split Legacy epsilon allocation \code{c(X, beta, sigma)} summing
#'   to 1.  Used only by \code{sens.method = "local"}.
#' @param B.y Public upper bound on \eqn{|y|}.  Used only by
#'   \code{sens.method = "local"}, where it sets the variance channel's noise
#'   scale, and \strong{required} there.  The sound \code{"subagg"} path needs no
#'   bound on \eqn{|y|} at all: its sensitivity comes from clipping the per-block
#'   coefficients into a public box, so it holds however large or unbounded
#'   \eqn{y} is.  That is the reason the mechanism exists.
#' @param use.copula Logical; if \code{TRUE} (default) and at least two
#'   covariates are continuous, the covariates' joint dependence is released
#'   under DP via a Gaussian copula (second-moment Gaussian mechanism on
#'   bounded scores, nearest-PSD post-processing) rather than sampling each
#'   marginal independently.  Since 0.6.0 the copula is mixed-margin: the
#'   categorical covariates join it through bounded indicator scores and the
#'   \code{back.map}, so their joint structure with the continuous block is
#'   preserved too (see the \emph{Mixed-margin copula} section).  With fewer
#'   than two continuous covariates the copula is a no-op, every marginal is
#'   sampled independently, and this is reported as \code{FALSE} in the
#'   result.
#' @param eps.corr.frac Fraction of \code{eps.X} spent on the DP correlation
#'   release when the copula is active; the rest funds the marginals.
#'   Default \code{0.5}.
#' @param score.method Copula score transform: \code{"adaptive"} (default, PIT
#'   through each covariate's released DP histogram, handling skewed or bounded
#'   marginals), \code{"normal"} (public normal-score PIT that recovers the
#'   covariate correlation at \eqn{\varepsilon \ge 5} on location-scale margins),
#'   or \code{"box"} (legacy, noise-dominated).  Affects only the
#'   \code{use.copula = TRUE} X-fidelity, not coefficient recovery.
#' @param back.map Back-map from a categorical's \eqn{K-1} copula latents to a
#'   level, under \code{use.copula = TRUE}: \code{"probit"} (default,
#'   recommended) assigns \code{argmax(0, Z_k + shift_k)}, deterministic, so
#'   it transfers the full released latent correlation; \code{"logit"} draws
#'   one level per row from \code{softmax(0, Z_k + shift_k)}, whose
#'   independent choice noise attenuates the recovered association
#'   (assumption A2 in the \emph{Mixed-margin copula} section).  Under either
#'   map the shifts are calibrated so the back-mapped proportions reproduce
#'   the released DP marginal.  Both maps are post-processing of the released
#'   noisy correlation and spend no privacy budget.
#' @param seed Optional integer RNG seed.  If \code{NULL} the global RNG
#'   state is used.
#' @return An object of class \code{synthexpln}: a list with
#'   \describe{
#'     \item{\code{$syn}}{Synthetic data frame (DP covariates + DP response).
#'       This is the DP release.  Everything below is evaluation-only.}
#'     \item{\code{$beta.hat}}{Original OLS coefficient vector (NOT private).}
#'     \item{\code{$beta.dp}}{DP coefficient vector.}
#'     \item{\code{$sigma.beta}}{Per-coefficient noise sd.}
#'     \item{\code{$sigma2.hat}}{Original OLS residual variance (NOT private).}
#'     \item{\code{$sigma2.dp}}{DP residual variance.}
#'     \item{\code{$sens.method}}{Which calibration ran.}
#'     \item{\code{$epsilon}, \code{$delta}}{As supplied.}
#'     \item{\code{$mu.total}}{Total mu-GDP level (subagg only).}
#'     \item{\code{$mu.channels}}{Per-channel mu \code{c(marg, corr, coef)}.}
#'     \item{\code{$mu.check}}{The three channels re-composed in quadrature.
#'       Equals \code{$mu.total}.}
#'     \item{\code{$subagg}}{The subagg intermediates.  \code{$beta.agg} and
#'       \code{$n.failed} inside it are \strong{NON-PRIVATE} and must not be
#'       composed downstream.}
#'     \item{\code{$dp.corr}}{DP correlation matrix used by the copula
#'       (dimension \eqn{p_\mathrm{total}}: the continuous columns followed by
#'       each categorical's \eqn{K-1} indicator columns), or \code{NULL} when
#'       the copula is inactive.}
#'     \item{\code{$delta.total}}{Composed delta.  Under \code{"subagg"} the
#'       whole release is calibrated to one \eqn{(\varepsilon, \delta)} through
#'       \eqn{\mu}, so delta no longer accumulates across channels: it was
#'       \code{2 * delta} under the copula on the legacy path.}
#'     \item{\code{$use.copula}}{Whether the Gaussian copula was applied.}
#'     \item{\code{$variant}}{\code{"V4"}.}
#'     \item{\code{$scope}}{\code{"joint-release"}.}
#'   }
#'   \code{$resid.norm.dp}, \code{$s.j} and \code{$epsilon.split} are legacy
#'   diagnostics and are \code{NA}/\code{NULL} under \code{"subagg"}.
#' @seealso \code{\link{gen.syn.dp.projected.subagg}} for the public-design
#'   route, \code{\link{dp.dispersion.public}} for the dispersion modes.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(200), age = rnorm(200, 50, 10),
#'                 sex = factor(sample(c("M","F"), 200, replace = TRUE)))
#' x.spec <- list(age = list(type = "continuous", bounds = c(0, 100)),
#'                sex = list(type = "categorical", levels = c("M", "F")))
#' # clip.lo/clip.hi and disp.lo/disp.hi are public boxes.  Here y is standard
#' # normal and independent of the covariates by construction, so the population
#' # coefficients are 0 and sigma^2 is 1: these are facts about the example's own
#' # DGP, not quantities read off d.
#' syn <- gen.syn.dp.full(d, y ~ age + sex, x.spec = x.spec,
#'                         epsilon = 5, delta = 1e-6,
#'                         clip.lo = -1, clip.hi = 1,
#'                         disp.lo = 0, disp.hi = 4, seed = 1)
#' nrow(syn$syn)
#' syn$mu.check   # == syn$mu.total, the three channels recomposed
gen.syn.dp.full <- function(d, formula, x.spec,
                             epsilon, delta = 1e-6,
                             sens.method = c("subagg", "local"),
                             clip.lo = NULL, clip.hi = NULL,
                             disp.lo = NULL, disp.hi = NULL, m = NULL,
                             mu.split = c(marg = 0.05, corr = 0.05, coef = 0.90),
                             eps.split = c(X = 0.1, beta = 0.8, sigma = 0.1),
                             B.y = NULL, use.copula = TRUE,
                             eps.corr.frac = 0.5, score.method = "adaptive",
                             back.map = c("probit", "logit"),
                             seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    sens.method  <- match.arg(sens.method)
    score.method <- match.arg(score.method, c("box", "normal", "adaptive"))
    back.map     <- match.arg(back.map)

    # ---- parse formula -------------------------------------------------------
    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)
    stopifnot(setequal(names(x.spec), x.nms))

    # ---- fit OLS on real data (non-private reference, evaluation only) --------
    mo         <- stats::lm(formula, d)
    beta.hat   <- stats::coef(mo)
    sigma2.hat <- stats::summary.lm(mo)[["sigma"]]^2
    m.X        <- stats::model.matrix(formula, d)
    n          <- nrow(m.X)
    p          <- ncol(m.X)

    # ---- decide copula eligibility (>= 2 continuous covariates) --------------
    v.is.cont  <- vapply(x.nms,
                         function(nm) x.spec[[nm]][["type"]] == "continuous",
                         logical(1))
    v.cont.nms <- x.nms[v.is.cont]
    do.copula  <- isTRUE(use.copula) && length(v.cont.nms) >= 2L

    # ---- allocate the budget -------------------------------------------------
    # subagg: every channel is a Gaussian mechanism, so the release composes in
    # one mu-GDP quadrature, mu^2 = mu_marg^2 + mu_corr^2 + mu_coef^2, converted
    # to (eps, delta) by the exact dual.  mu.split holds fractions of mu^2.  The
    # coefficient channel dominates by default: the analyst's refit recovers
    # beta_DP exactly, so X_syn's joint structure does not touch inferential
    # fidelity and the marginal channels are starved in its favour.  Override
    # towards marg/corr for releases that also need descriptive fidelity of X.
    #
    # When the copula is off, the correlation share is folded back into the
    # marginals in quadrature, so the total mu is preserved either way.  The k
    # marginal releases are k functions of the same row, so they too compose in
    # quadrature: each column gets mu_marg / sqrt(k).
    if (sens.method == "subagg") {
        if (is.null(clip.lo) || is.null(clip.hi))
            stop("sens.method = \"subagg\" requires the public per-coefficient ",
                 "bounds clip.lo and clip.hi", call. = FALSE)
        if (is.null(disp.lo) || is.null(disp.hi))
            stop("sens.method = \"subagg\" requires the public sigma^2 bounds ",
                 "disp.lo and disp.hi", call. = FALSE)
        if (abs(sum(mu.split) - 1) > 1e-8) stop("mu.split must sum to 1")
        mu.total <- .gdp.mu.from.eps.delta(epsilon, delta)
        mu.coef  <- mu.total * sqrt(mu.split[["coef"]])
        mu.corr  <- if (do.copula) mu.total * sqrt(mu.split[["corr"]]) else 0
        mu.marg  <- if (do.copula) mu.total * sqrt(mu.split[["marg"]])
                    else .gdp.compose(mu.total *
                             sqrt(c(mu.split[["marg"]], mu.split[["corr"]])))
        mu.X.each <- mu.marg / sqrt(length(x.spec))
        eps.X <- NA_real_; eps.beta <- NA_real_; eps.sigma <- NA_real_
        eps.corr <- NA_real_; eps.marg <- NA_real_; eps.marg.each <- NULL
    } else {
        warning("sens.method = \"local\" calibrates beta to the DFBETA LOCAL ",
                "sensitivity, which is a function of the realised data. The ",
                "release therefore carries NO formal (eps, delta)-DP guarantee ",
                "(retired 2026-07-07). Retained only to reproduce pre-0.2.0 ",
                "runs. Use sens.method = \"subagg\".", call. = FALSE)
        if (is.null(B.y))
            stop("B.y is required by sens.method = \"local\" and must be ",
                 "PUBLIC.\n",
                 "  It sets the variance channel's noise scale, so a bound taken\n",
                 "  from the sample makes the noise a function of the private\n",
                 "  data.  Supply a bound from the DOMAIN, not from the data.",
                 call. = FALSE)
        if (abs(sum(eps.split) - 1) > 1e-8) stop("eps.split must sum to 1")
        eps.X     <- epsilon * eps.split[["X"]]
        eps.beta  <- epsilon * eps.split[["beta"]]
        eps.sigma <- epsilon * eps.split[["sigma"]]
        eps.corr  <- if (do.copula) eps.X * eps.corr.frac else 0
        eps.marg  <- eps.X - eps.corr
        # Laplace(2/eps): the replace-one count sensitivity is 2.  Pure eps-DP,
        # so it does NOT enter a mu-GDP quadrature, which is one reason this
        # path is retired.
        eps.marg.each <- eps.marg / length(x.spec)
        mu.total <- NA_real_; mu.coef <- NA_real_
        mu.corr  <- NA_real_; mu.marg <- NA_real_; mu.X.each <- NULL
    }

    # channel dispatch: the marginal mechanisms take exactly one of eps or mu.
    dp.marg.cont <- function(x, breaks)
        if (sens.method == "subagg")
            .dp.hist.continuous(x, breaks, mu = mu.X.each)
        else .dp.hist.continuous(x, breaks, eps = eps.marg.each)
    dp.marg.cat <- function(x, levels)
        if (sens.method == "subagg") .dp.cat(x, levels, mu = mu.X.each)
        else .dp.cat(x, levels, eps = eps.marg.each)
    dp.marg.cat.probs <- function(x, levels)
        if (sens.method == "subagg") .dp.cat.probs(x, levels, mu = mu.X.each)
        else .dp.cat.probs(x, levels, eps = eps.marg.each)

    # ---- Stage 1: DP X marginals (+ optional mixed-margin Gaussian copula) ---
    if (do.copula) {
        v.cat.nms <- setdiff(x.nms, v.cont.nms)
        # DP marginals, all released: histograms (counts + breaks) for the
        # continuous columns, probability vectors for the categorical ones.
        # Since 0.6.0 the categoricals join the copula (mixed-margin)
        # rather than being sampled independently, so their released probs
        # feed both the indicator scores and the back-map calibration.
        l.cont.marg <- Map(function(nm) {
            spec   <- x.spec[[nm]]
            breaks <- if (!is.null(spec[["breaks"]])) spec[["breaks"]]
                      else seq(spec[["bounds"]][1], spec[["bounds"]][2],
                               length.out = 11)
            list(counts = dp.marg.cont(d[[nm]], breaks), breaks = breaks)
        }, v.cont.nms)
        l.cat.marg <- Map(function(nm)
            list(levels = x.spec[[nm]][["levels"]],
                 probs  = dp.marg.cat.probs(d[[nm]],
                                            x.spec[[nm]][["levels"]])),
            v.cat.nms)

        # Extended score matrix: continuous PIT scores, then K-1 bounded
        # indicator scores per categorical (from the released probs, so no
        # extra privacy cost).  Every column lands in [-s.M, s.M], so
        # .dp.copula.corr reads p = ncol(m.W) = p_total and the released
        # vech(S) sensitivity is M^2 p_total / n unchanged in form.  attained
        # only in the continuous sub-block (assumption A3).
        s.tau <- if (score.method == "box") 1 / (2 * n) else 0.005
        s.M   <- stats::qnorm(1 - s.tau)
        m.W.cont <- vapply(v.cont.nms, function(nm) {
            spec <- x.spec[[nm]]
            v.bx <- if (!is.null(spec[["breaks"]])) range(spec[["breaks"]])
                    else spec[["bounds"]]
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
        d.syn <- d.syn[x.nms]                       # restore declared order
    } else {
        m.R.dp <- NULL
        d.syn <- Reduce(function(acc, nm) {
            spec <- x.spec[[nm]]
            if (spec[["type"]] == "continuous") {
                breaks <- if (!is.null(spec[["breaks"]])) spec[["breaks"]]
                          else seq(spec[["bounds"]][1], spec[["bounds"]][2], length.out = 11)
                counts    <- dp.marg.cont(d[[nm]], breaks)
                col.value <- .sample.from.dp.hist(counts, breaks, n)
            } else {  # categorical
                col.value <- factor(dp.marg.cat(d[[nm]], spec[["levels"]]),
                                    levels = spec[["levels"]])
            }
            acc[[nm]] <- col.value
            acc
        }, x.nms, init = data.frame(row.names = seq_len(n)))
    }

    # ---- Stages 2 and 3: coefficient and dispersion --------------------------
    if (sens.method == "subagg") {
        # Subsample-and-aggregate on (beta, sigma^2) JOINTLY.  Sound under
        # full-row replacement, which is the relation that matters once X is
        # sensitive: the block partition is data-independent, and every
        # coordinate is clipped into a public box, so a one-row change moves the
        # mean by at most s_j = (hi_j - lo_j)/m whether it altered x, y or both.
        # The scaled aggregate has L2 sensitivity sqrt(p + 1), so the release is
        # mu_coef-GDP and composes in quadrature with the two X channels.
        #
        # This replaces both retired channels: the DFBETA local-sensitivity beta
        # release (data-dependent, no formal guarantee) and the Laplace
        # residual-norm sigma^2 release (valid only under response-row
        # replacement holding X fixed, so it did not cover the
        # sensitive-X case it shipped in).
        l.sa <- .subsample.aggregate.beta.loglink(
            d, formula, stats::gaussian(), "gaussian", "identity",
            eps = NA_real_, delta = delta, m = m,
            clip.lo = clip.lo, clip.hi = clip.hi,
            seed = if (is.null(seed)) 6489L else seed + 10L,
            mu = mu.coef, disp = TRUE, disp.lo = disp.lo, disp.hi = disp.hi)
        beta.dp       <- l.sa[["beta.dp"]]
        sigma2.dp     <- l.sa[["sigma2.dp"]]
        sigma.beta    <- l.sa[["sigma"]]
        resid.norm.dp <- NA_real_
        v.s           <- NULL
    } else {
        l.sa <- NULL
        # retired: DFBETA local sensitivity.  Data-dependent, so it delivers no
        # formal (eps, delta)-DP.  Reachable only through sens.method = "local",
        # which warns loudly above.  Kept so pre-0.2.0 runs stay reproducible.
        v.sd <- apply(m.X, 2, stats::sd)
        v.s  <- ifelse(v.sd > 1e-8, v.sd, 1)
        m.dfb    <- stats::dfbeta(mo)
        LS.tilde <- 2 * max(sqrt(rowSums(sweep(m.dfb, 2, v.s, "*")^2)))
        sigma.beta <- LS.tilde * .dp.gm.multiplier(eps.beta, delta)
        eta.tilde  <- stats::rnorm(p, 0, sigma.beta)
        # 3-sigma L2 truncation (post-processing)
        s.norm <- sqrt(sum(eta.tilde^2))
        s.max  <- 3 * sqrt(p) * sigma.beta
        if (s.norm > s.max) eta.tilde <- eta.tilde * (s.max / s.norm)
        beta.tilde <- beta.hat * v.s
        beta.dp    <- `names<-`((beta.tilde + eta.tilde) / v.s, names(beta.hat))

        # Legacy sigma^2: Laplace on the residual norm g = ||(I - H)y||_2, whose
        # one-row y sensitivity is 2 B.y for any fixed design.  sigma^2 is the
        # debiased square, E[(g + L)^2] = g^2 + 2 b^2, so subtract 2 b^2.
        b.sigma       <- 2 * B.y / eps.sigma
        resid.norm.dp <- sqrt(sum(stats::residuals(mo)^2)) +
            .rlaplace(1, 0, b.sigma)
        sigma2.dp     <- max((resid.norm.dp^2 - 2 * b.sigma^2) / (n - p), 1e-6)
    }

    # ---- Stage 4: build y using X_syn (not X_orig) --------------------------
    fml.rhs  <- stats::reformulate(x.nms)
    m.X.syn  <- stats::model.matrix(fml.rhs, d.syn)
    # Align beta.dp to columns of m.X.syn (missing factor levels set to 0)
    v.common      <- intersect(names(beta.dp), colnames(m.X.syn))
    beta.aligned  <- stats::setNames(rep(0, ncol(m.X.syn)), colnames(m.X.syn))
    beta.aligned[v.common] <- beta.dp[v.common]
    # Project fresh residuals onto null space of X_syn
    H.syn  <- m.X.syn %*%
        solve(crossprod(m.X.syn) + 1e-10 * diag(ncol(m.X.syn))) %*%
        t(m.X.syn)
    e.iid  <- stats::rnorm(n, 0, sqrt(sigma2.dp))
    e.proj <- as.numeric((diag(n) - H.syn) %*% e.iid)
    if (sum(e.proj^2) > 0)
        e.proj <- e.proj * sqrt((n - p) * sigma2.dp / sum(e.proj^2))
    d.syn[[y.nm]] <- as.numeric(m.X.syn %*% beta.aligned) + e.proj

    # ---- Return synthexpln object ----------------------------------------------
    # NOT private, evaluation only: beta.hat and sigma2.hat are the full-data
    # fits, and subagg carries its own non-private intermediates (beta.agg,
    # n.failed).  The DP release is $syn.  Do not compose the rest downstream.
    structure(list(
        syn           = d.syn,
        beta.hat      = beta.hat,     # NOT PRIVATE
        beta.dp       = beta.dp,
        sigma.beta    = sigma.beta,
        # The noise sd on the raw scale of each coefficient.  Under subagg the
        # noise is already raw and per-coefficient (s_j * sigma.Q), but it is
        # indexed over the AGGREGATED coordinates, which with disp = TRUE include
        # a (p+1)-th entry for sigma^2.  subset by name: renaming positionally
        # silently pairs the dispersion's noise with a coefficient.  The legacy
        # path perturbs the STANDARDISED beta, so its raw sd is sigma.beta / s_j.
        sigma.beta.raw = if (sens.method == "subagg")
                             sigma.beta[names(beta.dp)]
                         else `names<-`(as.numeric(sigma.beta) / v.s,
                                        names(beta.hat)),
        sigma2.hat    = sigma2.hat,   # NOT PRIVATE
        sigma2.dp     = sigma2.dp,
        sens.method   = sens.method,
        epsilon       = epsilon,
        delta         = delta,
        # mu-GDP accounting (subagg).  mu.check re-composes the three channels in
        # quadrature and must reproduce mu.total.
        mu.total      = mu.total,
        mu.channels   = c(marg = mu.marg, corr = mu.corr, coef = mu.coef),
        mu.check      = if (sens.method == "subagg")
                            .gdp.compose(c(mu.marg, mu.corr, mu.coef))
                        else NA_real_,
        subagg        = l.sa,
        # legacy (sens.method = "local") diagnostics
        resid.norm.dp = resid.norm.dp,
        s.j           = v.s,
        epsilon.split = if (sens.method == "local")
                            c(X = eps.X, beta = eps.beta, sigma = eps.sigma,
                              marg = eps.marg, corr = eps.corr)
                        else NULL,
        dp.corr       = m.R.dp,
        # Under subagg the whole release is calibrated to one (eps, delta)
        # through mu, so delta does not accumulate across channels as it did.
        delta.total   = if (sens.method == "subagg") delta
                        else if (do.copula) 2 * delta else delta,
        use.copula    = do.copula,
        B.y           = B.y,
        variant       = "V4",
        scope         = "joint-release"
    ), class = "synthexpln")
}

# ---------------------------------------------------------------------------
# Private helper: sample from DP histogram
# ---------------------------------------------------------------------------

#' Sample n values from a DP-noised histogram
#'
#' @param counts Non-negative numeric vector of DP-noised bin counts
#'   (as returned by \code{.dp.hist.continuous}).
#' @param breaks Numeric vector of bin boundaries (length B + 1).
#' @param n Integer number of draws.
#' @return Numeric vector of length \code{n} with values drawn uniformly
#'   within bins, weighted by \code{counts}.
#' @noRd
.sample.from.dp.hist <- function(counts, breaks, n) {
    B     <- length(counts)
    probs <- if (sum(counts) > 0) counts / sum(counts) else rep(1 / B, B)
    bins  <- sample.int(B, n, replace = TRUE, prob = probs)
    stats::runif(n, breaks[bins], breaks[bins + 1L])
}
