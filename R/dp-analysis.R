# dp-analysis.R -- V1 analysis-level DP generator for synthexpln.
# Public function:  gen.syn.dp.projected()
# Internal helper:  .inject.virtual.y()
# Ported (with simplifications) from lib/28-dp-projection.r lines 218-442.

# +------------------------------------------------------------------------+
# | .inject.virtual.y -- internal helper                                   |
# +------------------------------------------------------------------------+
# Constructs a virtual dataset d.virt whose MLE equals beta.dp by
# shifting the response by (mu(beta.dp) - mu(beta.hat)) at each row.
#
# Ported from inject.virtual.y() in lib/28-dp-projection.r lines 218-254.
# Factor-level alignment follows lines 224-232 of that source.
#
# Arguments:
#   d        Original data frame.
#   y.nm     Response column name (character).
#   beta.hat Original MLE coefficient vector (named).
#   beta.dp  DP-perturbed coefficient vector (named, same length).
#   family   Lowercase family string: "gaussian", "gamma", "poisson", etc.
#   link     Link string: "identity", "log", etc.
# Returns: modified data frame with updated response column.

#' Virtual response injection (internal)
#'
#' @section It retains the real y, so it must not build a DP release:
#' \code{y_virt = y + mu(beta_DP) - mu(beta_hat)} keeps the real response in the
#' output.  Anything built from it therefore depends on the real data by more
#' than \code{beta_DP}, and cannot claim the post-processing clause of the
#' subsample-and-aggregate proposition, whose hypothesis is a release depending
#' on the data "only through beta_DP".  Measured: the released residuals match
#' the real ones in scale to four significant figures and in shape too, and a
#' canary distinguishing game scores such a release at AUC 1.000 against a mu-GDP
#' bound of 0.566, flat in epsilon, because the leak did not go through the
#' noise.
#'
#' This is appropriate for \code{\link{gen.syn.dp.projected}}, the analysis-scope
#' generator, which makes no formal DP claim.  The sound routes use
#' \code{.synth.response} instead, which draws the response from the family at
#' \code{beta_DP} and so does not touch the real y.
#' @noRd
.inject.virtual.y <- function(d, y.nm, beta.hat, beta.dp, family, link) {
    v.x.nms <- setdiff(names(d), y.nm)
    fml <- stats::reformulate(v.x.nms, response = y.nm)
    m.X <- stats::model.matrix(fml, d)

    # Align beta.hat and beta.dp to columns of m.X.  Factor levels in
    # beta.hat/beta.dp may name columns that map to model.matrix columns;
    # any column not present in either vector is set to zero (conservative).
    v.beta.h.aligned <- stats::setNames(rep(0, ncol(m.X)), colnames(m.X))
    v.common.h <- intersect(names(beta.hat), colnames(m.X))
    v.beta.h.aligned[v.common.h] <- beta.hat[v.common.h]

    v.beta.d.aligned <- stats::setNames(rep(0, ncol(m.X)), colnames(m.X))
    v.common.d <- intersect(names(beta.dp), colnames(m.X))
    v.beta.d.aligned[v.common.d] <- beta.dp[v.common.d]

    v.eta.h  <- as.numeric(m.X %*% v.beta.h.aligned)
    v.eta.dp <- as.numeric(m.X %*% v.beta.d.aligned)

    # Apply inverse link to obtain mu on the response scale.
    # Only identity and log are needed for this Task 6 scope (Gaussian/OLS
    # default and canonical GLMs); further links can be added in later tasks.
    linkinv <- switch(link,
        identity = function(eta) eta,
        log      = exp,
        # fallback: treat as identity (rare for the supported families)
        function(eta) eta)
    v.mu.h  <- linkinv(v.eta.h)
    v.mu.dp <- linkinv(v.eta.dp)

    d.virt <- d
    d.virt[[y.nm]] <- d[[y.nm]] + v.mu.dp - v.mu.h

    # Poisson: round to integer lattice and clip to >= 0.
    if (family == "poisson")
        d.virt[[y.nm]] <- pmax(0L, round(d.virt[[y.nm]]))

    # Gamma / quasi-gamma: clip to small positive to avoid log(0).
    if (family %in% c("gamma", "quasi-gamma"))
        d.virt[[y.nm]] <- pmax(1e-8, d.virt[[y.nm]])

    d.virt
}

# +------------------------------------------------------------------------+
# | gen.syn.dp.projected -- public API                                     |
# +------------------------------------------------------------------------+

#' V1 analysis-level DP synthetic data generator
#'
#' Generates synthetic data whose refitted regression coefficients recover
#' the DP-perturbed estimate \eqn{\hat\beta_\mathrm{DP}}.  Recovery is to
#' IEEE 754 double-precision machine epsilon (\eqn{2^{-52} \approx 2.2 \times
#' 10^{-16}}, measured within 1 ulp) for Gaussian with the identity link,
#' where the response can be bent to any target without error.  For other
#' GLM families it is
#' approximate, and the error scales with the size of the DP perturbation
#' rather than sitting at the Newton tolerance (\eqn{\sim 10^{-7}}) that
#' \code{\link{projection}} reaches on the same families undisplaced.  It
#' falls as \eqn{\varepsilon} rises and the perturbation shrinks.  Budget for
#' this route on that basis, not on the undisplaced figure.  The DP
#' guarantee is analysis-level (Variant 1 in the plan): only the
#' regression coefficients are privatised; the covariate distribution is
#' released by the non-private projection step.
#'
#' @details
#' The algorithm is:
#' \enumerate{
#'   \item Fit the original model and obtain \eqn{\hat\beta}.
#'   \item Compute the local (DFBETA-based) sensitivity
#'     \eqn{LS = 2 \max_i \|\mathrm{DFBETA}_i\|_2}.
#'   \item Compute the Gaussian mechanism scale \eqn{\sigma = LS \cdot m},
#'     where \eqn{m} is the calibrated noise multiplier: the classical
#'     \eqn{\sqrt{2 \ln(1.25/\delta)} / \varepsilon} bound (Dwork-Roth
#'     2014, Thm A.1) for \eqn{\varepsilon < 1}, and the exact analytic
#'     Gaussian mechanism (Balle-Wang 2018) for \eqn{\varepsilon \ge 1},
#'     where the classical bound under-delivers \eqn{\delta}.
#'   \item Sample i.i.d. \eqn{N(0, \sigma^2)} noise and add it to
#'     \eqn{\hat\beta} to obtain \eqn{\hat\beta_\mathrm{DP}}.
#'   \item Construct a virtual dataset \eqn{d_\mathrm{virt}} by shifting
#'     the response by \eqn{\mu(\hat\beta_\mathrm{DP}) - \mu(\hat\beta)}.
#'   \item Run the non-DP projection on \eqn{d_\mathrm{virt}} so that the
#'     refitted synthetic coefficients equal \eqn{\hat\beta_\mathrm{DP}}.
#' }
#'
#' @section Simplifications:
#' Compared to the full research implementation, \code{28-dp-projection.r}
#' in the replication repository the paper's Data availability section
#' cites, this
#' implementation omits: (a) column standardisation of the design matrix
#' (so noise variance is not balanced across coefficients), (b) L2 norm
#' clipping of the noise vector, and (c) smooth sensitivity.  These
#' omissions leave the \eqn{(\varepsilon, \delta)}-DP guarantee intact but
#' may produce less efficient privatisation for log-link GLMs or highly
#' imbalanced covariates.  Standardisation and norm clipping are deferred
#' to a later release.
#'
#' @param d Data frame containing the response and predictors.
#' @param formula Model formula (e.g. \code{y ~ x1 + x2}).
#' @param epsilon Privacy budget \eqn{\varepsilon > 0}.
#' @param delta Privacy parameter \eqn{\delta \in (0, 1)}.
#'   Default \code{1e-6}.
#' @param family A GLM family object (default \code{gaussian()}).
#' @param seed Optional integer RNG seed passed to
#'   \code{\link{projection}}.  If \code{NULL} the global RNG state is
#'   used (reproducibility is the caller's responsibility).
#' @return An object of class \code{synthexpln}: a list with
#'   \describe{
#'     \item{\code{$syn}}{Synthetic data frame.}
#'     \item{\code{$beta.dp}}{DP-perturbed coefficient vector.}
#'     \item{\code{$beta.hat}}{Original MLE coefficient vector.}
#'     \item{\code{$sigma.beta}}{Per-coefficient Gaussian noise sd
#'       (scalar, since no column standardisation is applied).}
#'     \item{\code{$s.j}}{Column scales (all 1, no standardisation).}
#'     \item{\code{$epsilon}}{Epsilon used.}
#'     \item{\code{$delta}}{Delta used.}
#'     \item{\code{$variant}}{\code{"V1"}.}
#'     \item{\code{$scope}}{\code{"analysis-level"}.}
#'   }
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(200), x = rnorm(200))
#' syn <- gen.syn.dp.projected(d, y ~ x, epsilon = 1, delta = 1e-6)
#' coef(lm(y ~ x, syn$syn))   # == syn$beta.dp (up to floating point)
gen.syn.dp.projected <- function(d, formula,
                                  epsilon, delta = 1e-6,
                                  family = stats::gaussian(),
                                  seed = NULL) {
    # V1 is the retired analysis-scope route.  Its noise is calibrated to the
    # DFBETA local sensitivity, which is a function of the realised data, so the
    # release carries no formal (eps, delta)-DP at the stated epsilon or at any
    # other.  It previously made no statement about this, and dp.project() fell
    # through to it by default, so the package's headline entry point returned a
    # release with no guarantee, silently.  It now warns on every call, and
    # dp.project() no longer dispatches here.
    warning("gen.syn.dp.projected() (V1) calibrates its noise to the DFBETA ",
            "LOCAL sensitivity, which is a function of the realised data. The ",
            "release therefore carries NO formal (eps, delta)-DP guarantee at ",
            "the stated epsilon, or at any other (retired 2026-07-07). It is an ",
            "EMPIRICAL-DP baseline. For a release with a guarantee use one of ",
            "the three sound routes: gen.syn.dp.projected.subagg() (needs only ",
            "a public coefficient box), gen.syn.dp.ols.public() (needs a public ",
            "design and a public bound on |y|), or gen.syn.dp.suffstats() ",
            "(keeps the design private, needs a public whitening matrix, a ",
            "public row bound and a public response clip).", call. = FALSE)

    # ---- bridge formula / family to internal string representations ----
    # Follows the same switch() pattern as projection() in R/projection.R
    # lines 47-53 so that the two functions are consistent.
    y.nm <- all.vars(formula)[1]

    fam.nm.raw <- family[["family"]]
    fam.nm <- switch(fam.nm.raw,
        Gamma               = "gamma",
        `Negative Binomial` = "neg.binomial",
        `inverse.gaussian`  = "inv.gaussian",
        `Inverse Gaussian`  = "inv.gaussian",
        tolower(fam.nm.raw))   # gaussian, poisson, binomial, quasi*
    link.nm <- family[["link"]]

    # ---- Step 1: fit original model, extract beta_hat ------------------
    fit <- if (fam.nm == "gaussian" && link.nm == "identity") {
        stats::lm(formula, d)
    } else {
        stats::glm(formula, family = family, data = d)
    }
    beta.hat <- stats::coef(fit)

    # ---- Step 2: local sensitivity via DFBETA -------------------------
    LS <- sensitivity.local(fit, bound = "dfbeta")

    # ---- Step 3: Gaussian mechanism scale -----------------------------
    # Calibrated multiplier: classical for eps < 1, analytic for eps >= 1
    sigma <- LS * .dp.gm.multiplier(epsilon, delta)

    # ---- Step 4: sample noise and perturb beta ------------------------
    eta     <- stats::rnorm(length(beta.hat), 0, sigma)
    beta.dp <- `names<-`(beta.hat + eta, names(beta.hat))

    # ---- Step 5: inject virtual response ------------------------------
    d.virt <- .inject.virtual.y(d, y.nm, beta.hat, beta.dp,
                                 fam.nm, link.nm)

    # ---- Step 6: non-DP projection on the virtual dataset -------------
    proj <- projection(d.virt, formula, family = family, seed = seed)

    # ---- Return synthexpln object ----------------------------------------
    structure(list(
        syn        = proj[["syn"]],
        beta.dp    = beta.dp,
        beta.hat   = beta.hat,
        sigma.beta = sigma,
        # V1 does not column-standardise, so the noise is already raw-scale and
        # isotropic.  Carried under the same name as every other generator so
        # calibrated.se() does not have to guess which scale it is looking at.
        sigma.beta.raw = `names<-`(rep(sigma, length(beta.hat)),
                                   names(beta.hat)),
        s.j        = `names<-`(rep(1, length(beta.hat)), names(beta.hat)),
        epsilon    = epsilon,
        delta      = delta,
        variant    = "V1",
        scope      = "analysis-level"
    ), class = "synthexpln")
}

# +------------------------------------------------------------------------+
# | gen.syn.dp.projected.zi -- full DP-ZI with a private zero-model        |
# +------------------------------------------------------------------------+

#' Full differentially-private zero-inflated generator (private zero-model)
#'
#' Extends \code{\link{gen.syn.dp.projected}} to zero-inflated log-link
#' families by privately releasing the zero-model as well as the
#' positive-component coefficient.  It fits a logistic zero-model
#' \eqn{Z = I(Y > 0) \sim X} and privatises its coefficient to
#' \eqn{\hat\gamma_\mathrm{DP}} at \eqn{\mu_Z}-GDP, privatises the
#' positive-component coefficient to \eqn{\hat\beta_\mathrm{DP}} at
#' \eqn{\mu_P}-GDP, and composes the two in quadrature to
#' \eqn{\{\mu_Z^2 + \mu_P^2\}^{1/2}}-GDP.  Both coefficient channels are the
#' assumption-free subsample-and-aggregate Gaussian mechanism over a
#' \strong{public} per-coordinate clip box, so the release needs no bound on
#' \eqn{|y|} at all.
#'
#' @details
#' The construction is sound mu-Gaussian DP (Dong, Roth and Su, 2022,
#' \doi{10.1111/rssb.12454}):
#' \enumerate{
#'   \item DP zero-model \eqn{\hat\gamma_\mathrm{DP}} by
#'     subsample-and-aggregate over \code{m} data-independent blocks of the
#'     full design, at \eqn{\mu_Z}.
#'   \item DP positive coefficient \eqn{\hat\beta_\mathrm{DP}} by
#'     subsample-and-aggregate on the positive rows, at \eqn{\mu_P}.  The
#'     partition is taken over the \strong{fixed full \eqn{N}} rows and
#'     restricted to the positive subset, so a one-row replacement that toggles
#'     a row's positivity (and so changes the sensitive count \eqn{|P|}) still
#'     perturbs exactly one block.  Partitioning the positive subset directly
#'     would make the partition a function of \eqn{|P|} and void the guarantee.
#'   \item Private synthetic indicator
#'     \eqn{Z_\mathrm{syn} \sim \mathrm{Bernoulli}(\mathrm{expit}(x^\top
#'     \hat\gamma_\mathrm{DP}))}, a post-processing of \eqn{\hat\gamma_\mathrm{DP}}.
#'   \item Response drawn from the model's own family at
#'     \eqn{(\hat\beta_\mathrm{DP}, \phi)} with \eqn{\phi} the declared public
#'     dispersion, zeroed wherever \eqn{Z_\mathrm{syn} = 0}.  The real response
#'     is not an input to the draw, so the release depends on the data only
#'     through the two DP coefficients and post-processing applies.
#' }
#' The covariates are released unchanged (public design).  The clip boxes set
#' the per-coordinate sensitivity \eqn{s_j = (u_j - \ell_j)/m}, so they must be
#' public: a box read off the sample makes the noise a function of the private
#' data and voids the guarantee.
#'
#' @param d Data frame containing the zero-inflated response and predictors.
#' @param formula Model formula (log-link GLM, e.g. \code{y ~ x1 + x2}).
#' @param epsilon Total privacy budget \eqn{\varepsilon > 0}.
#' @param delta Total privacy parameter \eqn{\delta \in (0, 1)}.
#'   Default \code{1e-6}.
#' @param family A log-link GLM family object: \code{Gamma(link = "log")},
#'   \code{poisson(link = "log")}, or \code{quasi(link = "log", variance =
#'   "mu^2")} (the default, quasi-gamma).
#' @param clip.lo,clip.hi public per-coordinate bounds on the
#'   positive-component coefficients (named vectors aligned to the coefficient
#'   names, or scalars recycled).  Required.  The box width sets the noise
#'   scale, so it must not be derived from \code{d}.
#' @param clip.z.lo,clip.z.hi public per-coordinate bounds on the logistic
#'   zero-model coefficients.  Required.
#' @param dispersion A dispersion declaration from
#'   \code{\link{dp.dispersion.public}}, naming where the positive component's
#'   dispersion comes from.  Required for \code{gamma} and \code{quasi-gamma}
#'   (which carry a free dispersion); a private dispersion channel is not
#'   implemented here, so it must be public.
#' @param mu.split Named numeric \code{c(Z, P)} summing to 1, the fractions of
#'   the \eqn{\mu^2} budget for the zero-model and the positive component.
#'   Default \code{c(Z = 0.5, P = 0.5)}.
#' @param m Subsample-and-aggregate block count.  Default \code{NULL} applies
#'   the public rule \code{\link{dp.subagg.blocks}} per channel,
#'   \eqn{\lfloor n / (5 p) \rfloor} for the zero-model (fitted on every row)
#'   and \eqn{\lfloor n (1 - \pi_0) / (5 p) \rfloor} for the positive
#'   component, and then requires \code{pos.rate} and per-coefficient vector
#'   clip boxes (the rule reads \eqn{p} from the box length).  An explicit
#'   scalar \code{m} covers both channels.
#' @param pos.rate A \code{dp.rate} object from \code{\link{dp.rate.public}}
#'   carrying the \strong{published} zero rate \eqn{\pi_0}.  Required when
#'   \code{m} is left to the rule.  A rate computed from the sample, such as
#'   \code{mean(y == 0)}, makes the block count a function of the private data
#'   and voids the calibration.
#' @param seed Optional integer RNG seed.
#' @return An object of class \code{synthexpln}: a list with \code{$syn} (the
#'   DP release), \code{$beta.dp} and \code{$gamma.dp} (the two DP coefficient
#'   vectors), \code{$beta.hat} and \code{$gamma.hat} (the non-private
#'   full-data fits, for evaluation only), \code{$Z.syn}, \code{$sigma.beta}
#'   and \code{$sigma.gamma} (per-coordinate noise sds), \code{$dispersion},
#'   \code{$mu.total}, \code{$mu.Z}, \code{$mu.P}, \code{$m} (the named
#'   \code{c(z, p)} per-channel block counts), \code{$epsilon},
#'   \code{$delta}, \code{$variant = "ZI-subagg"}, and
#'   \code{$scope = "analysis-level"}.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 800
#' x1 <- rnorm(n); x2 <- rnorm(n)
#' z  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1))
#' mu <- exp(0.5 + 0.4 * x1 - 0.2 * x2)
#' y  <- ifelse(z == 1, rgamma(n, 2, scale = mu / 2), 0)
#' d  <- data.frame(y = y, x1 = x1, x2 = x2)
#' box <- c(`(Intercept)` = 8, x1 = 8, x2 = 8)
#' syn <- gen.syn.dp.projected.zi(
#'     d, y ~ x1 + x2, epsilon = 10,
#'     pos.rate = dp.rate.public(0.5, source = "design zero rate"),
#'     clip.lo = -box, clip.hi = box, clip.z.lo = -box, clip.z.hi = box,
#'     dispersion = dp.dispersion.public(1, source = "design residual variance"),
#'     family = Gamma(link = "log"), seed = 1)
#' mean(syn$syn$y == 0)   # private zero fraction
#' }
gen.syn.dp.projected.zi <- function(d, formula, epsilon, delta = 1e-6,
                                    family = stats::quasi(link = "log",
                                                          variance = "mu^2"),
                                    clip.lo = NULL, clip.hi = NULL,
                                    clip.z.lo = NULL, clip.z.hi = NULL,
                                    dispersion = NULL,
                                    mu.split = c(Z = 0.5, P = 0.5),
                                    m = NULL, pos.rate = NULL, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)

    fam.nm.raw <- family[["family"]]
    fam.nm <- switch(fam.nm.raw,
        Gamma = "gamma", poisson = "poisson", tolower(fam.nm.raw))
    if (grepl("quasi", fam.nm)) fam.nm <- "quasi-gamma"
    link.nm <- family[["link"]]

    if (is.null(clip.lo) || is.null(clip.hi))
        stop("clip.lo and clip.hi (public bounds on the positive-component ",
             "coefficients) are required")
    if (is.null(clip.z.lo) || is.null(clip.z.hi))
        stop("clip.z.lo and clip.z.hi (public bounds on the zero-model ",
             "coefficients) are required")
    if (.dp.has.dispersion(fam.nm) && is.null(dispersion))
        stop("dispersion is required for the ", fam.nm, " family. Supply ",
             "dp.dispersion.public(sigma2, source): the positive component has ",
             "a free dispersion parameter and DP must cover it.")
    if (!is.null(dispersion) && !inherits(dispersion, "dp.dispersion"))
        stop("dispersion must be built by dp.dispersion.public() or ",
             "dp.dispersion.private()")
    if (!is.null(dispersion) && dispersion[["mode"]] == "private")
        stop("dp.dispersion.private() is not implemented for the zero-inflated ",
             "generator; supply dp.dispersion.public().")
    if (link.nm != "log")
        stop("gen.syn.dp.projected.zi handles the log link only")
    if (!(fam.nm %in% c("gamma", "poisson", "quasi-gamma")))
        stop("family must be one of gamma, poisson, quasi-gamma")
    if (abs(sum(mu.split) - 1) > 1e-8) stop("mu.split must sum to 1")
    if (!any(d[[y.nm]] == 0))
        stop("no zeros in response; use gen.syn.dp.projected for non-ZI data")

    # Block-count rule, per channel.  The zero-model fits every row, the
    # positive fit sees only the expected positive share, and both partitions
    # cover the fixed n rows, so both read the public nrow(d).  The positive
    # share must be the published rate: reading it off d would calibrate the
    # mechanism to the data it protects.
    if (is.null(m)) {
        if (is.null(pos.rate))
            stop("with m = NULL (the block-count rule) pos.rate is required.  ",
                 "Supply dp.rate.public(pi0, source) with the PUBLISHED zero ",
                 "rate.", call. = FALSE)
        if (!inherits(pos.rate, "dp.rate"))
            stop("pos.rate must be built by dp.rate.public(), so that its ",
                 "provenance is on the record", call. = FALSE)
        if (length(clip.lo) < 2L || length(clip.z.lo) < 2L)
            stop("with m = NULL the clip boxes must be per-coefficient ",
                 "vectors: the rule reads p from the box length.",
                 call. = FALSE)
        m.z <- dp.subagg.blocks(nrow(d), length(clip.z.lo))
        m.p <- dp.subagg.blocks(nrow(d), length(clip.lo), rate = pos.rate)
    } else {
        m.z <- m
        m.p <- m
    }

    s.sigma2  <- if (is.null(dispersion)) NA_real_ else dispersion[["sigma2"]]
    seed.part <- if (is.null(seed)) 6489L else seed

    # Split the total mu-GDP budget in quadrature: mu_Z^2 + mu_P^2 = mu_total^2.
    s.mu.total <- .gdp.mu.from.eps.delta(epsilon, delta)
    s.mu.Z <- s.mu.total * sqrt(mu.split[["Z"]])
    s.mu.P <- s.mu.total * sqrt(mu.split[["P"]])

    v.pos <- d[[y.nm]] > 0
    d.pos <- d[v.pos, , drop = FALSE]

    # ---- Step 1: DP zero-model gamma_DP via subsample-and-aggregate ----
    # Logistic Z = I(y > 0) ~ X, privatised to mu_Z-GDP over the full design (a
    # fixed row set), the caller's public clip box setting the sensitivity.
    d.z <- d[x.nms]
    d.z[[".Z"]] <- as.integer(v.pos)
    l.z <- .subsample.aggregate.beta.loglink(
        d.z, stats::reformulate(x.nms, response = ".Z"),
        stats::binomial(link = "logit"), "binomial", "logit",
        eps = NA_real_, delta = delta, m = m.z,
        clip.lo = clip.z.lo, clip.hi = clip.z.hi,
        seed = seed.part + 100L, mu = s.mu.Z)
    v.gamma.dp <- l.z[["beta.dp"]]

    # ---- Step 2: DP positive-component beta_DP via subsample-and-aggregate ----
    # The positive block is fit on d.pos, whose row count |P| is sensitive.
    # Partition the fixed full N rows and restrict the assignment to the positive
    # rows (passed as `block`), so a one-row replacement that toggles a row's
    # positivity changes exactly one block and the (hi_j - lo_j)/m sensitivity
    # holds.  Partitioning d.pos directly would make the partition a function of
    # |P| and let a toggle reshuffle it, voiding mu_P-GDP (test-dp-zi.R).
    set.seed(seed.part + 3L)
    v.block.full <- sample(rep_len(seq_len(m.p), nrow(d)))
    l.p <- .subsample.aggregate.beta.loglink(
        d.pos, formula, family, fam.nm, "log",
        eps = NA_real_, delta = delta, m = m.p,
        clip.lo = clip.lo, clip.hi = clip.hi,
        seed = seed.part, mu = s.mu.P, block = v.block.full[v.pos])
    v.beta.hat <- l.p[["beta.hat"]]
    v.beta.dp  <- l.p[["beta.dp"]]

    # ---- Step 3: private zero indicator Z_syn ----
    # Z_syn ~ Bernoulli(expit(X gamma_DP)): post-processing of gamma_DP and the
    # public covariates.  X is released as it stands (public design).
    m.X.syn <- stats::model.matrix(stats::reformulate(x.nms), d)
    v.gamma.aligned <- `names<-`(rep(0, ncol(m.X.syn)), colnames(m.X.syn))
    v.common.z <- intersect(names(v.gamma.dp), colnames(m.X.syn))
    v.gamma.aligned[v.common.z] <- v.gamma.dp[v.common.z]
    v.p.syn <- stats::plogis(as.numeric(m.X.syn %*% v.gamma.aligned))
    set.seed(seed.part + 2L)
    v.Z.syn <- stats::rbinom(nrow(m.X.syn), 1L, v.p.syn)

    # ---- Step 4: draw the positive response from (beta_DP, dispersion) ----
    # Drawn from the model's own family, not from the real y, then zeroed
    # wherever the private zero-model switched the row off.  The release depends
    # on the data only through (beta_DP, gamma_DP), so post-processing applies.
    d.syn <- .synth.response(d, y.nm, x.nms, fam.nm, "log",
                             v.beta.dp, s.sigma2, seed.part)
    d.syn <- `[[<-`(d.syn, y.nm,
                    value = `[<-`(d.syn[[y.nm]], v.Z.syn == 0L, 0))

    # beta.hat and gamma.hat are non-private full-data fits, returned for
    # evaluation only.  They are not part of the release.
    structure(list(
        syn         = d.syn,
        beta.hat    = v.beta.hat,   # NOT PRIVATE
        beta.dp     = v.beta.dp,
        gamma.hat   = l.z[["beta.hat"]],   # NOT PRIVATE
        gamma.dp    = v.gamma.dp,
        Z.syn       = v.Z.syn,
        sigma.beta  = l.p[["sigma"]],
        sigma.gamma = l.z[["sigma"]],
        dispersion  = dispersion,
        mu.total    = .gdp.compose(c(s.mu.Z, s.mu.P)),
        mu.Z        = s.mu.Z,
        mu.P        = s.mu.P,
        m           = c(z = m.z, p = m.p),
        epsilon     = epsilon,
        delta       = delta,
        variant     = "ZI-subagg",
        scope       = "analysis-level"
    ), class = "synthexpln")
}
