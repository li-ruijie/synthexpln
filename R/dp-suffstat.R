# dp-suffstat.R -- V6: the whitened sufficient-statistic release (Route 3)
#
# Mirrors lib/66-dp-suffstat.r of the replication repository,
# with the public constants supplied by the caller rather than read from the
# published simulation DGP: a package has no scenario constants to consult, so
# the whitening matrix, the whitened row bound, and the response
# standardisation and clip arrive as arguments, each required and documented
# as public. The ridge penalty defaults to the released-eigenmin rule
# (AdaSSP, Wang 2018), the deployment mode for data without a published
# error curve; a numeric lambda with public provenance may be supplied
# instead. The analyst-side inference pair (beta.infer, Cov.beta.infer)
# recentres the calibrated interval at the lambda-0 solve of the released
# statistics, the 2026-07-16 coverage finding: the first-order ridge-bias
# fold reads half the true shrinkage at lambda = n and under-covers.

#' Differentially private whitened sufficient-statistic release (V6)
#'
#' Releases the three OLS sufficient statistics \eqn{\mathrm{vech}(X_w'X_w)},
#' \eqn{X_w'y_s} and \eqn{y_s'y_s} on a whitened public scale under one
#' \eqn{\mu}-GDP quadrature, then solves through a ridge penalty. The
#' whitened scale is \eqn{x \mapsto W x} with \eqn{W} a fixed public
#' positive-definite matrix, and \eqn{y_s} the publicly standardised,
#' publicly clipped response. The covariates stay \strong{private}: every
#' sensitivity is a functional of \eqn{(W, B_x, B_y)}, not of the realised
#' design, which is what the public-covariate release
#' (\code{\link{gen.syn.dp.ols.public}}) cannot offer.
#'
#' The \eqn{L_2} sensitivities under one-row replacement are
#' \eqn{\sqrt{2} B_x^2} for the vech channel (a proven strict upper bound on
#' a box containing no pair of one-hot rows of norm \eqn{B_x} with disjoint
#' supports), \eqn{2 B_x B_y} for \eqn{X_w'y_s}, and \eqn{B_y^2} for
#' \eqn{y_s'y_s} (attained). The three Gaussian channels compose in
#' quadrature to \code{mu.total}, and everything after the release, the
#' ridge solve, the PD repair, the back-map, the inference pair, and the
#' optional synthetic draw, reads only released and public objects, so it is
#' post-processing.
#'
#' \strong{The inference pair.} The ridge solve \code{beta.dp} is biased
#' toward zero, and folding a first-order bias correction into the standard
#' error under-covers at strong shrinkage (as low as 66\% against the 95\%
#' claimed, measured 2026-07-16). The calibrated interval is therefore
#' centred at \code{beta.infer}, the lambda-0 solve of the released
#' statistics through eigenvalues floored at the public \eqn{10^{-6} n},
#' with the covariance plug-in \code{Cov.beta.infer} evaluated at the
#' de-shrunk coefficient,
#' \deqn{SE_j^2 = SE_{syn,j}^2 + \mathrm{Cov.beta.infer}_{jj}.}
#' \code{beta.dp} stays the point estimate. \code{bias.fold} is the exact
#' de-shrink \code{beta.infer - beta.dp}.
#'
#' @param d Data frame with the response and covariates.
#' @param formula Model formula, e.g. \code{y ~ x1 + x2}.
#' @param epsilon,delta Privacy budget; converted once to \code{mu.total} by
#'   the exact \eqn{\mu}-GDP dual.
#' @param mu.split Named fractions of \eqn{\mu^2} for the channels
#'   \code{c(vech, xty, yy)}, summing to 1. With a public dispersion the
#'   \code{yy} share is dropped and the split renormalised.
#' @param W Public positive-definite whitening matrix on the model-matrix
#'   basis (intercept included), with dimnames matching the model matrix
#'   columns. \strong{Required, and it must be public}: the natural choice
#'   is \eqn{\Sigma^{-1/2}} for a published or prior second-moment matrix.
#'   A whitening estimated from \code{d} voids the guarantee.
#' @param B.x Public bound on \eqn{\|W x\|_2} over the public covariate
#'   support box. \strong{Required.} For a box, the maximum is attained at a
#'   vertex, so it is computable in closed form from public constants. The
#'   premise is enforced rather than assumed: any whitened row outside the
#'   ball is norm-clipped onto it, a deterministic public map that is the
#'   identity wherever the declared box truly bounds the support.
#' @param y.centre,y.scale Public affine standardisation of the response.
#' @param y.clip Public length-2 clip range for the response, applied before
#'   standardisation. \strong{Required}, and it defines \eqn{B_y}. Where the
#'   response is genuinely unbounded no public clip exists and this release
#'   is unavailable: use \code{\link{gen.syn.dp.projected.subagg}}.
#' @param lambda Either \code{"released-eigenmin"} (the default), which
#'   spends \code{eig.split} of the vech channel's \eqn{\mu^2} share
#'   releasing \eqn{\lambda_{\min}(S)} at the Weyl sensitivity
#'   \eqn{B_x^2} and sets the penalty from the released value, or a
#'   non-negative numeric with \strong{public} provenance. The penalty
#'   moves utility only, not the guarantee.
#' @param eig.split Fraction of the vech share diverted to the eigenvalue
#'   release under \code{"released-eigenmin"}.
#' @param dispersion \code{NULL} (default) privatises the variance through
#'   the \eqn{y_s'y_s} channel. A \code{\link{dp.dispersion.public}} object
#'   declares \eqn{\sigma^2} public instead, dropping that channel.
#' @param x.public If \code{TRUE}, additionally emits synthetic data on the
#'   \emph{original} design via the projection construction (the premise of
#'   \code{\link{gen.syn.dp.ols.public}}), whose OLS refit recovers
#'   \code{beta.dp} exactly. The default \code{FALSE} emits the
#'   statistics-only release.
#' @param seed Optional integer RNG seed. If \code{NULL} the global RNG
#'   state is used.
#' @return An object of class \code{synthexpln}: a list with
#'   \describe{
#'     \item{\code{$syn}}{Synthetic data frame under \code{x.public = TRUE},
#'       otherwise \code{NULL}.}
#'     \item{\code{$beta.dp}}{The ridge solve of the released statistics on
#'       the analyst's raw scale, the point estimate.}
#'     \item{\code{$beta.infer}, \code{$Cov.beta.infer}}{The inference pair:
#'       the lambda-0 solve and its plug-in covariance, for the calibrated
#'       interval.}
#'     \item{\code{$bias.fold}}{The exact de-shrink
#'       \code{beta.infer - beta.dp}.}
#'     \item{\code{$Cov.beta.dp}}{Plug-in conditional covariance of the
#'       ridge solve.}
#'     \item{\code{$sigma2.dp}}{Dispersion from the released statistics
#'       (floored at \eqn{10^{-6}}), or the declared public value.}
#'     \item{\code{$S.dp}, \code{$b.dp}, \code{$q.dp}}{The released triple
#'       (whitened scale).}
#'     \item{\code{$lambda}, \code{$lambda.eff}}{The penalty and its
#'       PD-repaired effective value.}
#'     \item{\code{$mu.total}, \code{$mu.channels}}{The \eqn{\mu}-GDP
#'       accounting; the channels re-compose in quadrature to
#'       \code{mu.total}.}
#'     \item{\code{$sigma.channels}, \code{$Delta.channels}}{Noise scales
#'       and sensitivities per channel.}
#'     \item{\code{$epsilon}, \code{$delta}}{As supplied.}
#'     \item{\code{$variant}}{\code{"V6"}.}
#'     \item{\code{$scope}}{\code{"statistic-triple"}, or
#'       \code{"joint-release"} under \code{x.public = TRUE}.}
#'   }
#' @references Dwork, C., Talwar, K., Thakurta, A. and Zhang, L. (2014).
#'   Analyze Gauss: optimal bounds for privacy-preserving principal
#'   component analysis. \emph{STOC 2014}, 11--20.
#'   Wang, Y.-X. (2018). Revisiting differentially private linear
#'   regression: optimal and adaptive prediction & estimation in unbounded
#'   domain. \emph{UAI 2018}, 93--103.
#' @export
#' @examples
#' set.seed(1)
#' n <- 400
#' d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
#' W <- diag(3)
#' dimnames(W) <- list(c("(Intercept)", "x1", "x2"),
#'                     c("(Intercept)", "x1", "x2"))
#' syn <- gen.syn.dp.suffstats(d, y ~ x1 + x2, epsilon = 5,
#'     W = W, B.x = sqrt(1 + 2 * 6^2), y.centre = 0, y.scale = 1,
#'     y.clip = c(-8, 8), x.public = TRUE, seed = 1)
#' coef(lm(y ~ x1 + x2, syn$syn))   # == syn$beta.dp (up to floating point)
gen.syn.dp.suffstats <- function(d, formula, epsilon, delta = 1e-6,
                                 mu.split = c(vech = 0.45, xty = 0.45,
                                              yy = 0.10),
                                 W = NULL, B.x = NULL,
                                 y.centre = NULL, y.scale = NULL,
                                 y.clip = NULL,
                                 lambda = "released-eigenmin",
                                 eig.split = 0.10,
                                 dispersion = NULL,
                                 x.public = FALSE, seed = NULL) {
    y.nm  <- all.vars(formula)[1]
    x.nms <- setdiff(all.vars(formula), y.nm)
    fml   <- stats::reformulate(x.nms, response = y.nm)

    # ---- the public constants are required -------------
    # A generator does not invent its own bounds. Each of these sets a noise
    # scale, so a value derived from d would make the noise a function of the
    # private data and the guarantee would not hold at any epsilon.
    if (is.null(W) || is.null(B.x) || is.null(y.centre) || is.null(y.scale) ||
        is.null(y.clip))
        stop("W, B.x, y.centre, y.scale, and y.clip are required, and every ",
             "one must be PUBLIC.\n",
             "  They set the noise scales (sqrt(2) B_x^2, 2 B_x B_y, B_y^2), ",
             "so a value\n",
             "  derived from the sample voids the guarantee at any epsilon. ",
             "Supply them\n",
             "  from the domain or a published source, not from the data.",
             call. = FALSE)
    if (abs(sum(mu.split) - 1) > 1e-12 ||
        !identical(sort(names(mu.split)), c("vech", "xty", "yy")))
        stop("mu.split must be named fractions (vech, xty, yy) summing to 1",
             call. = FALSE)
    b.disp.public <- inherits(dispersion, "dp.dispersion")
    if (!is.null(dispersion) && !b.disp.public)
        stop("dispersion must be NULL or a dp.dispersion.public object",
             call. = FALSE)
    if (b.disp.public && dispersion[["mode"]] != "public")
        stop("dispersion must be NULL (private y'y channel) or a ",
             "dp.dispersion.public object", call. = FALSE)
    b.eigmin <- is.character(lambda)
    if (b.eigmin) lambda <- match.arg(lambda, "released-eigenmin")
    if (!b.eigmin && (!is.numeric(lambda) || length(lambda) != 1L ||
                      lambda < 0))
        stop("lambda must be \"released-eigenmin\" or a single non-negative ",
             "PUBLIC numeric", call. = FALSE)
    if (b.disp.public) {
        mu.split[["yy"]] <- 0
        mu.split <- mu.split / sum(mu.split)
    }

    # ---- the whitened, publicly clipped statistics ---------------------------
    m.X <- stats::model.matrix(fml, d)
    if (!identical(colnames(m.X), colnames(W)))
        stop("W's dimnames must match the model matrix columns: ",
             paste(colnames(m.X), collapse = ", "), call. = FALSE)
    n <- nrow(m.X); p <- ncol(m.X)
    m.Xw <- m.X %*% W
    # The bounded-universe premise is enforced, not assumed: any whitened row
    # outside the public ball is norm-clipped onto it, a deterministic public
    # per-row map costing no budget, exactly as the response clip below. On a
    # support the declared box truly bounds, it is the identity map, and the
    # estimand is the clipped one either way.
    v.nrm <- sqrt(rowSums(m.Xw^2))
    if (any(v.nrm > B.x))
        m.Xw <- m.Xw * pmin(1, B.x / v.nrm)
    y.c  <- pmin(pmax(d[[y.nm]], y.clip[[1L]]), y.clip[[2L]])
    y.s  <- (y.c - y.centre) / y.scale
    B.y  <- max(abs((y.clip - y.centre) / y.scale))
    m.S  <- crossprod(m.Xw)
    v.b  <- as.numeric(crossprod(m.Xw, y.s))
    q    <- sum(y.s^2)

    # ---- the mu accounting ---------------------------------------------------
    Delta.vech <- sqrt(2) * B.x^2
    Delta.xty  <- 2 * B.x * B.y
    Delta.yy   <- B.y^2
    mu.total <- .gdp.mu.from.eps.delta(epsilon, delta)
    v.f <- mu.split
    f.eig <- 0
    if (b.eigmin) {
        f.eig <- eig.split * v.f[["vech"]]
        v.f[["vech"]] <- v.f[["vech"]] - f.eig
    }
    sig.S <- Delta.vech / (sqrt(v.f[["vech"]]) * mu.total)
    sig.b <- Delta.xty  / (sqrt(v.f[["xty"]])  * mu.total)
    sig.q <- if (b.disp.public) NA_real_ else
        Delta.yy / (sqrt(v.f[["yy"]]) * mu.total)

    # ---- the release: three Gaussian channels, one quadrature ---------------
    unvech <- function(v) {
        m <- matrix(0, p, p)
        m[upper.tri(m, diag = TRUE)] <- v
        m + t(m) - diag(diag(m))
    }
    if (!is.null(seed)) set.seed(seed + 1L)
    m.S.dp <- m.S + unvech(stats::rnorm(p * (p + 1) / 2, 0, sig.S))
    if (!is.null(seed)) set.seed(seed + 2L)
    v.b.dp <- v.b + stats::rnorm(p, 0, sig.b)
    q.dp <- NA_real_
    if (!b.disp.public) {
        if (!is.null(seed)) set.seed(seed + 3L)
        q.dp <- q + stats::rnorm(1L, 0, sig.q)
    }

    # ---- the penalty (post-processing either way) ----------------------------
    s.lambda <- if (b.eigmin) {
        sig.eig <- B.x^2 / (sqrt(f.eig) * mu.total)
        if (!is.null(seed)) set.seed(seed + 4L)
        eigmin.dp <- min(eigen(m.S.dp, symmetric = TRUE,
                               only.values = TRUE)[["values"]]) +
            stats::rnorm(1L, 0, sig.eig)
        max(0, 2 * sig.S * sqrt(p) - max(eigmin.dp, 0))
    } else lambda

    # ---- PD repair and the ridge solve ---------------------------------------
    eig.min <- min(eigen(m.S.dp + s.lambda * diag(p), symmetric = TRUE,
                         only.values = TRUE)[["values"]])
    lambda.eff <- s.lambda + max(0, 1e-6 * n - eig.min)
    m.M <- solve(m.S.dp + lambda.eff * diag(p))
    v.beta.w <- as.numeric(m.M %*% v.b.dp)

    # Back to the analyst's raw scale: X W beta_w = X (W beta_w), so the raw
    # coefficient is W beta_w for any public W, and the response
    # standardisation undoes through the scale, its centre landing on the
    # intercept.
    back <- function(v.w) {
        v <- as.numeric(W %*% v.w) * y.scale
        v[[1L]] <- v[[1L]] + y.centre
        `names<-`(v, colnames(m.X))
    }
    beta.dp <- back(v.beta.w)

    # ---- dispersion -----------------------------------------------------------
    sigma2.dp <- if (b.disp.public) dispersion[["sigma2"]] else {
        rss.w <- q.dp - 2 * sum(v.beta.w * v.b.dp) +
            as.numeric(t(v.beta.w) %*% m.S.dp %*% v.beta.w)
        max(rss.w * y.scale^2 / (n - p), 1e-6)
    }

    # ---- the analyst-side error law and the inference pair -------------------
    cov.at <- function(m.inv, v.w) {
        m.cw <- m.inv %*% (sig.b^2 * diag(p) +
                           sig.S^2 * (sum(v.w^2) * diag(p) +
                                      tcrossprod(v.w))) %*% m.inv
        (y.scale^2) * W %*% m.cw %*% W
    }
    m.cov.raw <- cov.at(m.M, v.beta.w)
    l.e0 <- eigen(m.S.dp, symmetric = TRUE)
    m.S0inv <- l.e0[["vectors"]] %*%
        (t(l.e0[["vectors"]]) / pmax(l.e0[["values"]], 1e-6 * n))
    v.beta.w0 <- as.numeric(m.S0inv %*% v.b.dp)
    beta.infer <- back(v.beta.w0)
    m.cov.infer <- cov.at(m.S0inv, v.beta.w0)

    # ---- optional synthetic emission under the public-X premise --------------
    d.syn <- NULL
    if (isTRUE(x.public)) {
        if (!is.null(seed)) set.seed(seed + 5L)
        e.iid <- stats::rnorm(n, 0, sqrt(sigma2.dp))
        m.H <- m.X %*% solve(crossprod(m.X)) %*% t(m.X)
        e.proj <- as.numeric((diag(n) - m.H) %*% e.iid)
        rss <- sum(e.proj^2)
        if (rss > 0) e.proj <- e.proj * sqrt((n - p) * sigma2.dp / rss)
        d.syn <- d
        d.syn[[y.nm]] <- as.numeric(m.X %*% beta.dp) + e.proj
    }

    out <- list(
        syn            = d.syn,
        beta.dp        = beta.dp,
        beta.infer     = beta.infer,
        Cov.beta.infer = m.cov.infer,
        bias.fold      = beta.infer - beta.dp,
        Cov.beta.dp    = m.cov.raw,
        sigma2.dp      = sigma2.dp,
        S.dp           = m.S.dp,
        b.dp           = v.b.dp,
        q.dp           = q.dp,
        lambda         = s.lambda,
        lambda.eff     = lambda.eff,
        mu.total       = mu.total,
        mu.channels    = sqrt(c(v.f, eig = f.eig)) * mu.total,
        sigma.channels = c(vech = sig.S, xty = sig.b, yy = sig.q),
        Delta.channels = c(vech = Delta.vech, xty = Delta.xty, yy = Delta.yy),
        B.x            = B.x,
        B.y            = B.y,
        epsilon        = epsilon,
        delta          = delta,
        n              = n,
        p              = p,
        variant        = "V6",
        scope          = if (isTRUE(x.public)) "joint-release"
                         else "statistic-triple")
    class(out) <- "synthexpln"
    out
}
