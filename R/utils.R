# Internal utility helpers for synthexpln.
# None of these are exported; all are marked @noRd.

# +------------------------------------------------------------------------+
# | .fit.copula                                                            |
# +------------------------------------------------------------------------+

#' Fit a Gaussian copula to a data frame
#'
#' Computes pairwise Pearson correlations of all numeric columns, constructs
#' a \code{normalCopula} object, and returns both the copula and the
#' correlation matrix (with column names preserved).
#'
#' @param d Data frame with numeric columns.
#' @return List with \code{$copula} (normalCopula S4 object) and \code{$R}
#'   (correlation matrix with row and column names equal to \code{names(d)}).
#' @noRd
#' @importFrom copula normalCopula P2p
.fit.copula <- function(d) {
    nms <- names(d)
    p   <- length(nms)
    R   <- stats::cor(d, use = "pairwise.complete.obs")
    R   <- `rownames<-`(R, nms)
    R   <- `colnames<-`(R, nms)
    # P2p converts the lower-triangle of a correlation matrix to the
    # parameter vector expected by normalCopula.
    param <- copula::P2p(R)
    cop   <- copula::normalCopula(param = param, dim = p, dispstr = "un")
    list(copula = cop, R = R)
}

# +------------------------------------------------------------------------+
# | std.cumulants (private helper for .cf.quantile)                       |
# +------------------------------------------------------------------------+

#' Standardised cumulants of a sample
#'
#' Extracts standardised cumulants gamma_3, ..., gamma_order from a sample
#' vector \code{v}.  Returns a zero vector if the sample is too small or
#' degenerate.  Ported verbatim from source \code{gen-syn-family.r} (line 160).
#'
#' @param v Numeric vector.
#' @param order Integer; highest cumulant order requested (3--6).
#' @return Numeric vector of length \code{max(order - 2, 0)}.
#' @noRd
.std.cumulants <- function(v, order) {
    s.n     <- length(v)
    s.sigma <- stats::sd(v)
    if (s.sigma < .Machine[["double.eps"]] || s.n < order + 1)
        return(rep(0, max(order - 2, 0)))
    v.c  <- (v - mean(v)) / s.sigma
    v.mu <- mapply(function(r) mean(v.c^r), seq(3, order))
    v.gamma <- numeric(length(v.mu))
    if (order >= 3) v.gamma[1] <- v.mu[1]
    if (order >= 4) v.gamma[2] <- v.mu[2] - 3
    if (order >= 5) v.gamma[3] <- v.mu[3] - 10 * v.mu[1]
    if (order >= 6) v.gamma[4] <- v.mu[4] - 15 * v.mu[2] -
                                   10 * v.mu[1]^2 + 30
    v.gamma
}

# +------------------------------------------------------------------------+
# | .cf.quantile                                                           |
# +------------------------------------------------------------------------+

#' Cornish-Fisher quantile function factory
#'
#' Given a sample \code{x}, estimates its mean, standard deviation, and
#' standardised cumulants up to \code{order}, then returns a quantile
#' function (closure) that maps probabilities in (0, 1) to quantiles via
#' the Cornish-Fisher expansion (the partial moment matching theorem of
#' Li, 2026).
#'
#' The CF expansion reproduces the first \code{order} cumulants by
#' construction.  Ported from \code{cf.quantile} in source
#' \code{gen-syn-family.r} (line 194), wrapped as a closure over the
#' estimated parameters.
#'
#' @param x Numeric vector (the sample used to estimate parameters).
#' @param order Integer; expansion order (3--6).
#' @return A function \code{function(p)} that returns CF quantiles at
#'   probability vector \code{p}.
#' @noRd
.cf.quantile <- function(x, order) {
    mu    <- mean(x)
    sigma <- stats::sd(x)
    gamma <- .std.cumulants(x, order)

    function(p) {
        z <- stats::qnorm(p)
        w <- z
        if (order >= 3) {
            g1 <- gamma[1]
            w  <- w + (g1 / 6) * (z^2 - 1)
        }
        if (order >= 4) {
            g2 <- gamma[2]
            w  <- w + (g2 / 24) * (z^3 - 3 * z) -
                      (g1^2 / 36) * (2 * z^3 - 5 * z)
        }
        if (order >= 5) {
            g3 <- gamma[3]
            w  <- w + (g3 / 120) * (z^4 - 6 * z^2 + 3) -
                      (g1 * g2 / 24) * (z^4 - 5 * z^2 + 2) +
                      (g1^3 / 324) * (12 * z^4 - 53 * z^2 + 17)
        }
        if (order >= 6) {
            g4 <- gamma[4]
            w  <- w + (g4 / 720) * (z^5 - 10 * z^3 + 15 * z) -
                      (g1 * g3 / 180) * (2 * z^5 - 17 * z^3 + 21 * z) -
                      (g2^2 / 192) * (z^5 - 7 * z^3 + 9 * z) +
                      (g1^2 * g2 / 216) * (4 * z^5 - 32 * z^3 + 41 * z) -
                      (g1^4 / 3888) * (14 * z^5 - 118 * z^3 + 155 * z)
        }
        mu + sigma * w
    }
}
