# dp-public-ols.R -- V2 full-DP OLS synthesiser for public covariates.
# Public function: gen.syn.dp.ols.public()
# Ported from lib/33-dp-projection-public-x.r lines 39-133.

#' V2 full-DP OLS synthetic data generator (public covariates)
#'
#' Generates synthetic data using a fully differentially private OLS
#' mechanism under the assumption that the covariate matrix \eqn{X} is
#' public (non-private).  Both \eqn{\hat\beta} and \eqn{\hat\sigma^2} are
#' privatised; the synthetic response is constructed as
#' \eqn{y_\mathrm{syn} = X\hat\beta_\mathrm{DP} + e_\mathrm{proj}}, where
#' \eqn{e_\mathrm{proj}} is projected onto the null space of \eqn{X} so
#' that re-fitting OLS on the synthetic data recovers \eqn{\hat\beta_\mathrm{DP}}
#' exactly.
#'
#' @details
#' Both channels are Gaussian mechanisms at a \strong{global} sensitivity, so
#' the release composes in one mu-GDP quadrature,
#' \eqn{\mu^2 = \mu_\beta^2 + \mu_v^2}, which is converted to
#' \eqn{(\varepsilon, \delta)} by the exact Gaussian-DP dual.  \code{eps.split}
#' gives the fractions of the \eqn{\mu^2} budget.
#' \itemize{
#'   \item \strong{Coefficient.}  \eqn{\hat\beta = (X'X)^{-1}X'y}.  Replacing
#'     one row's response moves \eqn{X'y} by \eqn{x_i(y_i' - y_i)} and nothing
#'     else, so
#'     \eqn{\Delta_\beta = 2 B_y \max_i \|(X'X)^{-1} x_i\|_2}, a function of the
#'     public design and the public bound alone.  Noise is isotropic on the raw
#'     coefficient scale at \eqn{\sigma_\beta = \Delta_\beta / \mu_\beta}.
#'   \item \strong{Dispersion.}  \eqn{\Delta_{\sigma^2} = 4 B_y^2 K_X / (n - p)}
#'     with \eqn{K_X = \max_i \max\{1 - h_{ii}, \sum_{j \neq i} |h_{ij}|\}} the
#'     hat-row constant of the public design.  The naive bound without
#'     \eqn{K_X} under-states the one-row worst case whenever \eqn{K_X > 1}.
#' }
#'
#' @section What changed, and why the old mechanism had no theorem:
#' Before version 0.2.0 the coefficient channel was calibrated to the DFBETA
#' \strong{local} sensitivity
#' \eqn{\widetilde{LS} = 2\max_i\|\mathrm{DFBETA}_i \cdot s\|_2}, on a
#' column-standardised scale and with a 3-sigma L2 truncation.  Local
#' sensitivity is a function of the realised data, so two neighbouring datasets
#' received two different noise scales, no sensitivity bound held across the
#' pair, and there was no theorem to instantiate.  The release was not private at
#' any epsilon, and supplying a public \code{B.y} did not repair it, since
#' \code{B.y} did not reach that channel: it fed only the variance channel.
#'
#' The column standardisation and the truncation went with it.  Both were props
#' for the DFBETA noise (the standardisation balanced it across columns of very
#' different scale, the truncation tamed its tail); a Gaussian mechanism at a
#' global sensitivity already bounds the raw-scale L2 movement of
#' \eqn{\hat\beta}, so it needs neither.  The variance channel was also Laplace,
#' which is pure eps-DP and cannot enter a mu-GDP quadrature at all; it is
#' Gaussian now, so the two channels compose.
#'
#' @param d Data frame containing response and predictors.
#' @param formula Model formula (e.g. \code{y ~ x1 + x2}).
#' @param epsilon Total privacy budget \eqn{\varepsilon > 0}.
#' @param delta Privacy parameter \eqn{\delta \in (0, 1)}.
#'   Default \code{1e-6}.
#' @param eps.split Named numeric vector with entries \code{"beta"} and
#'   \code{"sigma"} summing to 1, giving the fraction of the \eqn{\mu^2} budget
#'   allocated to \eqn{\hat\beta} and \eqn{\hat\sigma^2} respectively.  These are
#'   fractions of \eqn{\mu^2}, not of \eqn{\varepsilon}: the two Gaussian
#'   channels compose in quadrature, so it is \eqn{\mu^2} that is additive.
#'   Default \code{c(beta = 0.8, sigma = 0.2)}.
#' @param B.y Public upper bound on \eqn{|y|}. \strong{Required.} There is no
#'   default and there cannot be one: \eqn{B_y} defines
#'   \eqn{\Delta_\beta = 2 B_y \max_i \|(X'X)^{-1} x_i\|}, so it sets the noise
#'   scale. A bound read off the sample makes that scale a function of the
#'   private data, and the guarantee then does not hold at any epsilon. Supply a
#'   bound from the domain (a physiological range, a reporting cap, a
#'   structurally bounded response). Where the response is genuinely unbounded,
#'   no \eqn{B_y} exists and this release is unavailable: use
#'   \code{\link{gen.syn.dp.projected.subagg}}, whose sensitivity comes from
#'   clipping the per-block coefficients into a public box and which therefore
#'   needs no bound on \eqn{|y|}.
#' @param seed Optional integer RNG seed. If \code{NULL} the global RNG
#'   state is used.
#' @return An object of class \code{synthexpln}: a list with
#'   \describe{
#'     \item{\code{$syn}}{Synthetic data frame with privatised response.  This
#'       is the DP release.  Everything below is evaluation-only and must not
#'       be composed downstream.}
#'     \item{\code{$beta.hat}}{Original OLS coefficient vector (NOT private).}
#'     \item{\code{$beta.dp}}{DP-perturbed coefficient vector.}
#'     \item{\code{$sigma.beta}}{Gaussian sd on the raw coefficient scale,
#'       \eqn{\Delta_\beta / \mu_\beta}.  Isotropic across coefficients.}
#'     \item{\code{$sigma.v}}{Gaussian sd of the dispersion channel.}
#'     \item{\code{$sigma2.hat}}{Original OLS residual variance (NOT private).}
#'     \item{\code{$sigma2.dp}}{DP-perturbed residual variance.}
#'     \item{\code{$Delta.beta}, \code{$Delta.sigma}}{The two global
#'       sensitivities actually used.}
#'     \item{\code{$epsilon}, \code{$delta}}{As supplied.}
#'     \item{\code{$mu.total}}{Total mu-GDP level, the exact dual of
#'       \eqn{(\varepsilon, \delta)}.}
#'     \item{\code{$mu.beta}, \code{$mu.v}}{Per-channel mu.}
#'     \item{\code{$mu.check}}{The two channels re-composed in quadrature.
#'       Equals \code{$mu.total}.}
#'     \item{\code{$B.y}}{Bound on \eqn{|y|} used.}
#'     \item{\code{$variant}}{\code{"V2"}.}
#'     \item{\code{$scope}}{\code{"joint-release"}.}
#'   }
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(200), x = rnorm(200))
#' syn <- gen.syn.dp.ols.public(d, y ~ x, epsilon = 1, delta = 1e-6,
#'                               B.y = 5)
#' coef(lm(y ~ x, syn$syn))   # == syn$beta.dp (up to floating point)
gen.syn.dp.ols.public <- function(d, formula,
                                   epsilon, delta = 1e-6,
                                   eps.split = c(beta = 0.8, sigma = 0.2),
                                   B.y = NULL, seed = NULL) {
    # ---- parse formula ------------------------------------------------------
    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)
    fml   <- stats::reformulate(x.nms, response = y.nm)

    # ---- B.y is required, and it must be public ------------------------------
    # This used to read
    #     if (is.null(B.y)) { B.y <- max(abs(d[[y.nm]])) * 1.5; message(...) }
    # which silently substituted the sample maximum for the public bound the
    # theorem requires, behind a message() that is invisible in any batch run.
    # B.y defines Delta_beta = 2 B_y max_i ||(X'X)^-1 x_i||, so it sets the noise
    # scale: a bound read off the sample makes the noise scale a function of the
    # private data, no valid sensitivity bound exists, and there is no THEOREM.
    # Not a weaker epsilon, and not any epsilon.  A default that quietly voids the
    # guarantee is worse than no default, so the call is now refused.
    if (is.null(B.y))
        stop("B.y is required and must be PUBLIC.\n",
             "  It sets the noise scale (Delta_beta = 2 B_y max_i ||(X'X)^-1 x_i||),\n",
             "  so a bound taken from the sample makes the noise a function of the\n",
             "  private data and the DP guarantee does not hold at any epsilon.\n",
             "  Supply a bound from the DOMAIN, not from the data.\n",
             "  If the response is genuinely unbounded, no B.y exists and this\n",
             "  release is unavailable: use gen.syn.dp.projected.subagg(), which\n",
             "  needs no bound on |y|.", call. = FALSE)

    # ---- allocate the mu-GDP budget -----------------------------------------
    # eps.split holds fractions of the mu^2 budget, not of epsilon: the two
    # Gaussian channels compose in quadrature, so it is mu^2 that is additive.
    eps.split <- eps.split / sum(eps.split)   # guard against non-unit sums
    mu.total  <- .gdp.mu.from.eps.delta(epsilon, delta)
    mu.beta   <- mu.total * sqrt(eps.split[["beta"]])
    mu.v      <- mu.total * sqrt(eps.split[["sigma"]])

    # ---- Step 1: fit OLS ----------------------------------------------------
    m.X      <- stats::model.matrix(fml, d)
    n        <- nrow(m.X)
    p        <- ncol(m.X)
    mo       <- stats::lm(fml, d)
    beta.hat <- stats::coef(mo)
    sigma2.hat <- stats::summary.lm(mo)[["sigma"]]^2

    # ---- Step 2: DP beta, Gaussian mechanism at the global sensitivity ------
    # Delta_beta = 2 B.y max_i ||(X'X)^-1 x_i||, a function of the public design
    # and the public bound and of nothing else, so it holds across every
    # neighbouring pair.  This was the DFBETA local sensitivity until 0.2.0, on a
    # standardised scale and with a 3-sigma truncation; both were props for that
    # noise and neither is needed here.  The noise is isotropic on the raw scale.
    Delta.beta <- .gs.beta.ols.public(m.X, B.y)
    sigma      <- Delta.beta / mu.beta

    if (!is.null(seed)) set.seed(seed + 1L)
    beta.dp <- `names<-`(beta.hat + stats::rnorm(p, 0, sigma), names(beta.hat))

    # ---- Step 3: DP sigma^2, Gaussian mechanism (was Laplace) ---------------
    # Gaussianised so it enters the same mu quadrature as the beta channel: a
    # Laplace release is pure eps-DP and does not compose in mu at all.  The
    # 1e-6 floor is the numerical form of the truncation at zero.
    Delta.sigma <- .gs.sigma2.public(m.X, B.y)
    sigma.v     <- Delta.sigma / mu.v

    if (!is.null(seed)) set.seed(seed + 2L)
    sigma2.dp <- max(sigma2.hat + stats::rnorm(1, 0, sigma.v), 1e-6)

    # ---- Step 4: draw iid Gaussian residuals --------------------------------
    if (!is.null(seed)) set.seed(seed + 3L)
    e.iid <- stats::rnorm(n, 0, sqrt(sigma2.dp))

    # ---- Step 5: project onto null(X) ---------------------------------------
    # H = X (X'X)^{-1} X'; (I - H) e.iid is orthogonal to col(X)
    m.H    <- m.X %*% solve(crossprod(m.X)) %*% t(m.X)
    e.proj <- as.numeric((diag(n) - m.H) %*% e.iid)

    # Rescale so that RSS = (n - p) * sigma2.dp
    rss <- sum(e.proj^2)
    if (rss > 0)
        e.proj <- e.proj * sqrt((n - p) * sigma2.dp / rss)

    # ---- Step 6: construct synthetic response --------------------------------
    y.syn  <- as.numeric(m.X %*% beta.dp) + e.proj
    d.syn  <- d
    d.syn[[y.nm]] <- y.syn

    # ---- Return synthexpln object ----------------------------------------------
    # NOT private, evaluation only: beta.hat and sigma2.hat are the full-data
    # fits.  The DP release is $syn.  Do not compose the rest downstream.
    #
    # mu.check re-composes the two channels in quadrature and must reproduce
    # mu.total; it is asserted in the test suite rather than merely reported.
    structure(list(
        syn           = d.syn,
        beta.hat      = beta.hat,
        beta.dp       = beta.dp,
        sigma.beta    = sigma,
        # the noise sd on the raw scale of each coefficient.  Isotropic here, so
        # every coordinate gets the same one.  Used by calibrated.se(): there is
        # no column scale to divide by, and dividing by one would understate the
        # uncertainty of a private release.
        sigma.beta.raw = `names<-`(rep(sigma, p), names(beta.hat)),
        sigma.v       = sigma.v,
        sigma2.hat    = sigma2.hat,
        sigma2.dp     = sigma2.dp,
        Delta.beta    = Delta.beta,
        Delta.sigma   = Delta.sigma,
        epsilon       = epsilon,
        delta         = delta,
        mu.total      = mu.total,
        mu.beta       = mu.beta,
        mu.v          = mu.v,
        mu.check      = .gdp.compose(c(mu.beta, mu.v)),
        B.y           = B.y,
        variant       = "V2",
        scope         = "joint-release"
    ), class = "synthexpln")
}
