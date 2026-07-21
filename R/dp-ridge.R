# dp-ridge.R -- V3 full-DP ridge synthesiser for public covariates.
# Public function: gen.syn.dp.ridge.public()
# Ported (with math correction) from lib/34-dp-ridge-projection.r.
#
# Note on the mathematics.  The source (lines 120-132) projects iid residuals via the ridge
# hat matrix H = X A^{-1} X' (A = X'X + lambda*I) and then sets
# y_syn = X beta_dp + (I - H) e_iid.  This does NOT guarantee that
# re-fitting ridge on y_syn recovers beta_dp exactly (see derivation in
# task description).  The corrected construction below uses the OLS null-space
# projection so that X' e_null = 0 exactly, and adds a lambda-correction to
# the fitted part so that X' y_syn = A beta_dp, giving exact recovery.

# ---------------------------------------------------------------------------
# Retired: ridge DFBETA local sensitivity
# ---------------------------------------------------------------------------
# Why this function is still here.  It is kept to demonstrate that the idea
# does not work.  A reader who wonders whether a local sensitivity could
# calibrate this mechanism can run the retired quantity against the sound one
# and watch it fail, which the test suite does on every run by asserting that
# this sensitivity moves between neighbouring datasets.  The failure is
# therefore demonstrated rather than asserted, and the function cannot be
# reinstated by someone who believes it was sound.  Keeping it also leaves the
# pre-0.2.0 releases reproducible.
#
# It is unexported, and no sound path may call it.
#
# Why it fails.  It is the local sensitivity of beta_ridge, read off the
# realised residuals (v.e <- v.y - v.mu below).  A sensitivity that is a
# function of the data it is meant to protect gives neighbouring datasets
# different noise scales, so no valid bound exists and there is no theorem to
# instantiate.  Not a weaker epsilon, and not any epsilon.
#
# The analysis-level proposition resting on this was retired on 2026-07-07, and
# V2 was migrated to the global sensitivity the same day.  V3 was not, in the
# replication repository until 2026-07-13 and in this package until 0.2.0,
# because the earlier repair migrated only its sigma^2 channel (to
# .gs.sigma2.ridge, sound) and left the beta channel behind.  That is the
# identical half-migration later found in V4 and V5: three generators, one
# mistake, three separate discoveries.
#
# V3 had no excuse.  V4 and V5 reached for DFBETA because their design is
# sensitive, so a public Delta_beta (a function of X) would itself leak the
# design and the global route was genuinely unavailable.  V3's design is
# public.  The ridge global sensitivity was always available and it is two
# lines, see .gs.beta.ridge.public in R/sensitivity.R.

#' Ridge DFBETA local sensitivity (retired, not a DP calibration)
#'
#' Computes the local sensitivity of the ridge coefficient vector at a given
#' dataset via the closed-form ridge DFBETA formula:
#'   beta_ridge(-i) - beta_ridge = -(X'X + lambda I)^{-1} x_i (y_i - x_i' beta) / (1 - h_ii)
#' where h_ii = x_i' (X'X + lambda I)^{-1} x_i.
#'
#' Ported from \code{sens.beta.ridge} in \file{lib/34-dp-ridge-projection.r}.
#' Retired there on 2026-07-13 and here at 0.2.0.  Do not calibrate noise to it.
#'
#' @param m.X Design matrix (already standardised by caller if needed).
#' @param v.y Response vector.
#' @param lambda Ridge penalty (positive).
#' @return Named list: \code{beta} (ridge coefficient vector) and
#'   \code{LS} (local sensitivity = 2 * max row-wise DFBETA L2 norm).
#' @noRd
.sens.beta.ridge <- function(m.X, v.y, lambda) {
    m.XtX <- crossprod(m.X)
    m.inv <- solve(m.XtX + lambda * diag(ncol(m.X)))
    v.beta <- m.inv %*% crossprod(m.X, v.y)
    v.mu <- as.numeric(m.X %*% v.beta)
    v.e <- v.y - v.mu
    v.h <- diag(m.X %*% m.inv %*% t(m.X))
    m.dfb <- sweep(m.X %*% m.inv, 1, v.e / (1 - v.h), "*")
    list(beta = as.numeric(v.beta), LS = 2 * max(sqrt(rowSums(m.dfb^2))))
}

# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

#' V3 full-DP ridge synthetic data generator (public covariates)
#'
#' Generates synthetic data using a fully differentially private ridge
#' regression mechanism under the assumption that the covariate matrix
#' \eqn{X} is public.  Both \eqn{\hat\beta_\mathrm{ridge}} and
#' \eqn{\hat\sigma^2} are privatised; the synthetic response is constructed
#' so that re-fitting ridge (with the same \code{lambda}) on the synthetic
#' data recovers \eqn{\hat\beta_\mathrm{DP}} to IEEE 754 double-precision
#' machine epsilon (\eqn{2^{-52} \approx 2.2 \times 10^{-16}}, measured
#' within 1 ulp).
#'
#' @details
#' Both channels are Gaussian mechanisms at a \strong{global} sensitivity, so
#' the release composes in one mu-GDP quadrature,
#' \eqn{\mu^2 = \mu_\beta^2 + \mu_v^2}, converted to
#' \eqn{(\varepsilon, \delta)} by the exact Gaussian-DP dual.
#' \itemize{
#'   \item \strong{Coefficient.}  With \eqn{A = X'X + \lambda I},
#'     \eqn{\hat\beta_\mathrm{ridge} = A^{-1}X'y}, and replacing one row's
#'     response gives
#'     \eqn{\Delta_\beta = 2 B_y \max_i \|A^{-1} x_i\|_2}: a function of the
#'     public design, the public \code{lambda} and the public bound.  This is
#'     the OLS global sensitivity with \eqn{A} in place of \eqn{X'X}, and it
#'     reduces to it exactly as \eqn{\lambda \to 0}.
#'   \item \strong{Dispersion.}  \eqn{\Delta_{\sigma^2} = 4 B_y^2 K_M / (n - p)},
#'     where \eqn{K_M = \max_i \max\{M_{ii}, \sum_{j \neq i} |M_{ij}|\}} is the
#'     row constant of \eqn{M = (I - H_\lambda)'(I - H_\lambda)} and
#'     \eqn{H_\lambda} is the ridge hat matrix of the public design.
#' }
#'
#' @section What changed, and why the old mechanism had no theorem:
#' Before version 0.2.0 the coefficient channel was calibrated to the ridge
#' DFBETA \strong{local} sensitivity, read off the realised residuals, on a
#' column-standardised scale and with a 3-sigma truncation.  A data-dependent
#' sensitivity gives neighbouring datasets different noise scales, so no bound
#' holds across the pair and the release was not private at any epsilon.
#'
#' V3 was the generator with no excuse.  The sensitive-design variants reached
#' for local sensitivity because a public \eqn{\Delta_\beta} is a function of
#' \eqn{X} and so would leak a private design.  V3's design is \strong{public},
#' so the global bound above was always available.  It was missed because an
#' earlier repair migrated only this function's variance channel and left the
#' coefficient channel behind.  The standardisation and truncation were props
#' for the DFBETA noise and went with it; the variance channel was Laplace,
#' which cannot enter a mu quadrature, and is Gaussian now.
#'
#' Corrected residual construction (deviation from source):
#' The research implementation \code{34-dp-ridge-projection.r} projects iid
#' residuals
#' via the ridge hat matrix, which does \emph{not} guarantee exact ridge
#' recovery.  This implementation uses the OLS null-space projection
#' \eqn{e_\mathrm{null} = (I - H_\mathrm{OLS}) e_\mathrm{iid}} (so that
#' \eqn{X' e_\mathrm{null} = 0} exactly) and a lambda-corrected fitted part
#' \eqn{X_\mathrm{fitted} = X (X'X)^{-1} A \hat\beta_\mathrm{DP}} (where
#' \eqn{A = X'X + \lambda I}), ensuring \eqn{X' y_\mathrm{syn} = A \hat\beta_\mathrm{DP}}.
#'
#' @param d Data frame containing response and predictors.
#' @param formula Model formula (e.g. \code{y ~ x1 + x2}).
#' @param lambda Ridge penalty parameter (positive scalar).
#' @param epsilon Total privacy budget \eqn{\varepsilon > 0}.
#' @param delta Privacy parameter \eqn{\delta \in (0, 1)}.
#'   Default \code{1e-6}.
#' @param eps.split Named numeric vector with entries \code{"beta"} and
#'   \code{"sigma"} summing to 1, giving the fraction of the \eqn{\mu^2} budget
#'   allocated to \eqn{\hat\beta} and \eqn{\hat\sigma^2} respectively.  These are
#'   fractions of \eqn{\mu^2}, not of \eqn{\varepsilon}: the two Gaussian
#'   channels compose in quadrature.  A scalar is interpreted as the beta
#'   fraction and \code{1 - scalar} as the sigma fraction.
#'   Default \code{c(beta = 0.8, sigma = 0.2)}.
#' @param B.y Public upper bound on \eqn{|y|}. \strong{Required.} There is no
#'   default and there cannot be one: \eqn{B_y} sets the noise scale of both
#'   released channels, the coefficient vector through \eqn{\Delta_\beta} and
#'   the residual variance through
#'   \eqn{\Delta_{\sigma^2} = 4 B_y^2 K_M / (n - p)}. A bound read off the
#'   sample makes those scales a function of the private data, and the
#'   guarantee then does not hold at any epsilon. Supply a bound from the
#'   domain (a physiological range, a reporting cap, a structurally bounded
#'   response). Where the response is genuinely unbounded, no \eqn{B_y} exists
#'   and this release is unavailable: use
#'   \code{\link{gen.syn.dp.projected.subagg}}, whose sensitivity comes from
#'   clipping the per-block coefficients into a public box and which therefore
#'   needs no bound on \eqn{|y|}.
#' @param seed Optional integer RNG seed. If \code{NULL} the global RNG
#'   state is used.
#' @return An object of class \code{synthexpln}: a list with
#'   \describe{
#'     \item{\code{$syn}}{Synthetic data frame with privatised response.}
#'     \item{\code{$beta.ridge}}{Original ridge coefficient vector.}
#'     \item{\code{$beta.dp}}{DP-perturbed ridge coefficient vector.}
#'     \item{\code{$sigma.beta}}{Gaussian sd on the raw coefficient scale,
#'       \eqn{\Delta_\beta / \mu_\beta}.  Isotropic across coefficients.}
#'     \item{\code{$sigma.v}}{Gaussian sd of the dispersion channel.}
#'     \item{\code{$sigma2.hat}}{Original ridge residual variance (NOT private).}
#'     \item{\code{$sigma2.dp}}{DP-perturbed residual variance.}
#'     \item{\code{$Delta.beta}, \code{$Delta.sigma}}{The two global
#'       sensitivities actually used.}
#'     \item{\code{$lambda}}{Ridge penalty used.}
#'     \item{\code{$epsilon}, \code{$delta}}{As supplied.}
#'     \item{\code{$mu.total}}{Total mu-GDP level, the exact dual of
#'       \eqn{(\varepsilon, \delta)}.}
#'     \item{\code{$mu.beta}, \code{$mu.v}}{Per-channel mu.}
#'     \item{\code{$mu.check}}{The two channels re-composed in quadrature.
#'       Equals \code{$mu.total}.}
#'     \item{\code{$B.y}}{Bound on \eqn{|y|} used.}
#'     \item{\code{$variant}}{\code{"V3"}.}
#'     \item{\code{$scope}}{\code{"joint-release"}.}
#'   }
#'   \code{$s.j} and \code{$LS} are gone at 0.2.0: there is no column
#'   standardisation and no local sensitivity.  \code{$epsilon.split} is
#'   replaced by the \code{$mu.*} fields, since the accounting is a mu
#'   quadrature and not an epsilon partition.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(200), x = rnorm(200))
#' syn <- gen.syn.dp.ridge.public(d, y ~ x, lambda = 1,
#'                                 epsilon = 1, delta = 1e-6, B.y = 5)
#' X.s <- model.matrix(~ x, syn$syn)
#' solve(crossprod(X.s) + diag(1, 2), crossprod(X.s, syn$syn$y))  # == syn$beta.dp
gen.syn.dp.ridge.public <- function(d, formula, lambda,
                                     epsilon, delta = 1e-6,
                                     eps.split = c(beta = 0.8, sigma = 0.2),
                                     B.y = NULL, seed = NULL) {
    if (lambda <= 0) stop("lambda must be > 0 for ridge")

    # ---- parse formula ------------------------------------------------------
    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)
    fml   <- stats::reformulate(x.nms, response = y.nm)

    # ---- B.y is required, and it must be public ------------------------------
    # It used to default to 1.5 * max|y|, the sample maximum, behind a message().
    # B.y sets the noise scale through Delta_beta and Delta_sigma2 (the latter via
    # 4 B_y^2 K_M / (n - p)), so a bound taken from the sample makes the noise a
    # function of the private data and voids the guarantee entirely.  See
    # gen.syn.dp.ols.public for the full statement.
    v.y <- d[[y.nm]]
    if (is.null(B.y))
        stop("B.y is required and must be PUBLIC.\n",
             "  It sets the noise scale, so a bound taken from the sample makes\n",
             "  the noise a function of the private data and the DP guarantee\n",
             "  does not hold at any epsilon.  Supply a bound from the DOMAIN.\n",
             "  If the response is unbounded, use gen.syn.dp.projected.subagg(),\n",
             "  which needs no bound on |y|.", call. = FALSE)

    # ---- allocate the mu-GDP budget -----------------------------------------
    # eps.split holds fractions of the mu^2 budget, not of epsilon.  Both
    # channels are Gaussian now and compose in one quadrature,
    # mu^2 = mu_beta^2 + mu_v^2, matching V2.  The sigma^2 channel was Laplace,
    # which is pure eps-DP and does not enter a mu quadrature at all, so the
    # release could not state the mu-GDP level it claimed.
    if (length(eps.split) == 1L) {
        eps.split <- c(beta = eps.split, sigma = 1 - eps.split)
    }
    eps.split <- eps.split / sum(eps.split)   # guard against non-unit sums
    mu.total  <- .gdp.mu.from.eps.delta(epsilon, delta)
    mu.beta   <- mu.total * sqrt(eps.split[["beta"]])
    mu.v      <- mu.total * sqrt(eps.split[["sigma"]])

    # ---- Step 1: design matrix and ridge fit --------------------------------
    m.X <- stats::model.matrix(fml, d)
    n   <- nrow(m.X)
    p   <- ncol(m.X)

    XtX     <- crossprod(m.X)
    A       <- XtX + lambda * diag(p)
    A.inv   <- solve(A)
    beta.ridge <- `names<-`(as.numeric(A.inv %*% crossprod(m.X, v.y)), colnames(m.X))

    v.mu.r     <- as.numeric(m.X %*% beta.ridge)
    sigma2.hat <- sum((v.y - v.mu.r)^2) / (n - p)

    # ---- Step 2: DP beta, Gaussian mechanism at the global sensitivity ------
    # Delta_beta = 2 B.y max_i ||(X'X + lambda I)^-1 x_i||, a function of the
    # public design, the public lambda and the public bound.  This was the ridge
    # DFBETA local sensitivity until 0.2.0, read off the realised residuals.
    #
    # The column standardisation went with it: it existed to balance the DFBETA
    # noise across columns of very different scale, and the global sensitivity is
    # already a bound on the raw-scale L2 movement of beta, so the noise is
    # isotropic on the raw scale exactly as V2 does it.  The 3-sigma truncation
    # went too, having been a post-processing adjustment to the local-sensitivity noise, and
    # a Gaussian mechanism at a global sensitivity needs no such assumption.
    Delta.beta <- .gs.beta.ridge.public(m.X, lambda, B.y)
    sigma      <- Delta.beta / mu.beta

    if (!is.null(seed)) set.seed(seed + 1L)
    beta.dp <- `names<-`(beta.ridge + stats::rnorm(p, 0, sigma), colnames(m.X))

    # ---- Step 3: DP sigma^2, Gaussian mechanism (was Laplace) ---------------
    # Delta_v = 4 B.y^2 K_M / (n - p), K_M the row constant of
    # M = (I - H_lambda)'(I - H_lambda), computed from the public X.
    Delta.sigma <- .gs.sigma2.ridge(m.X, lambda, B.y)
    sigma.v     <- Delta.sigma / mu.v

    if (!is.null(seed)) set.seed(seed + 2L)
    sigma2.dp <- max(sigma2.hat + stats::rnorm(1, 0, sigma.v), 1e-6)

    # ---- Step 6: corrected residual construction ----------------------------
    # Source's ridge-hat-matrix projection does NOT guarantee exact ridge
    # recovery; we instead use OLS null-space projection plus a
    # lambda-correction in the fitted part (derivation in task notes).
    #
    # Target: X' y_syn = A beta_dp so that analyst's ridge gives beta_dp.
    # X.fitted = X (X'X)^{-1} A beta_dp   =>  X' X.fitted = A beta_dp  [exact]
    # e.null   = (I - H_OLS) e_iid         =>  X' e.null   = 0          [exact]
    # => X' y_syn = A beta_dp + 0 = A beta_dp  QED
    XtX.inv <- solve(XtX)                              # OLS inverse
    X.fitted <- as.numeric(m.X %*% (beta.dp + lambda * XtX.inv %*% beta.dp))
    H.OLS    <- m.X %*% XtX.inv %*% t(m.X)

    if (!is.null(seed)) set.seed(seed + 3L)
    e.iid  <- stats::rnorm(n, 0, sqrt(sigma2.dp))
    e.null <- as.numeric((diag(n) - H.OLS) %*% e.iid)

    # Rescale residual so RSS = (n - p) * sigma2.dp
    rss <- sum(e.null^2)
    if (rss > 0)
        e.null <- e.null * sqrt((n - p) * sigma2.dp / rss)

    y.syn      <- X.fitted + e.null
    d.syn      <- d
    d.syn[[y.nm]] <- y.syn

    # ---- Return synthexpln object ----------------------------------------------
    # NOT private, evaluation only: beta.ridge and sigma2.hat are the full-data
    # fits.  The DP release is $syn.
    #
    # $LS is GONE with the local-sensitivity channel that produced it, and so are
    # $s.j and $epsilon.split: there is no standardisation and the accounting is
    # mu-GDP, not an epsilon partition.  mu.check re-composes the two channels.
    structure(list(
        syn           = d.syn,
        beta.hat      = beta.ridge,   # alias for synthexpln API consistency
        beta.ridge    = beta.ridge,
        beta.dp       = beta.dp,
        sigma.beta    = sigma,
        # raw-scale noise sd per coefficient (isotropic).  See gen.syn.dp.ols.public.
        sigma.beta.raw = `names<-`(rep(sigma, p), colnames(m.X)),
        sigma.v       = sigma.v,
        sigma2.hat    = sigma2.hat,
        sigma2.dp     = sigma2.dp,
        Delta.beta    = Delta.beta,
        Delta.sigma   = Delta.sigma,
        lambda        = lambda,
        epsilon       = epsilon,
        delta         = delta,
        mu.total      = mu.total,
        mu.beta       = mu.beta,
        mu.v          = mu.v,
        mu.check      = .gdp.compose(c(mu.beta, mu.v)),
        B.y           = B.y,
        variant       = "V3",
        scope         = "joint-release"
    ), class = "synthexpln")
}
