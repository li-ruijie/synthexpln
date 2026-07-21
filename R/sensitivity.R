#' Local sensitivity of a fitted model (a diagnostic, NOT a DP calibration)
#'
#' Computes the DFBETA-based local sensitivity
#' \eqn{LS(\hat\beta, D) \le 2 \max_i \|\mathrm{DFBETA}_i\|_2}.
#'
#' @section Do not calibrate a DP release to this:
#' Local sensitivity is a function of the realised data.  Two neighbouring
#' datasets get two different values, hence two different noise scales, so no
#' sensitivity bound holds across the neighbouring pair and no DP theorem can
#' be instantiated.  A release noised against it is not weakly private, it is
#' \strong{not private at any epsilon}.  The analysis-level proposition that
#' rested on this quantity was retired on 2026-07-07, and every generator that
#' MAKES a formal DP claim was migrated off it: the public-design routes to a
#' global sensitivity (\code{\link{gen.syn.dp.ols.public}},
#' \code{\link{gen.syn.dp.ridge.public}}), the sensitive-design routes to
#' subsample-and-aggregate (\code{\link{gen.syn.dp.full}},
#' \code{\link{gen.syn.dp.logistic}}).
#'
#' \strong{Two callers remain, and neither claims a guarantee.}
#' \code{\link{gen.syn.dp.projected}} (V1) calibrates to this quantity still,
#' so it is the retired no-guarantee baseline and warns on every
#' call, and \code{sens.method = "local"} on V4 and V5 reaches the legacy path
#' so that pre-0.2.0 runs can be reproduced.  Read "migrated off it" as scoped
#' to the sound routes, rather than as a claim that this is now uncalled.
#'
#' The quantity itself is perfectly real and worth reporting.  It measures how
#' much one row of \emph{this} dataset moves the fit, which is a useful
#' influence diagnostic and the natural yardstick for the empirical-versus-
#' worst-case comparison.  This function is retained for that, and for that
#' only.
#'
#' @param fit A fitted \code{lm} or \code{glm} object.
#' @param bound Method: \code{"dfbeta"} (default, leverage-aware) or
#'   \code{"leverage"} (naive upper bound).
#' @return A positive scalar.
#' @seealso \code{\link{gen.syn.dp.projected.subagg}} for a coefficient
#'   release whose sensitivity is data-independent.
#' @export
#' @examples
#' fit <- lm(mpg ~ wt, mtcars)
#' sensitivity.local(fit)   # an influence diagnostic, not a noise scale
sensitivity.local <- function(fit, bound = c("dfbeta", "leverage")) {
    bound <- match.arg(bound)
    if (bound == "dfbeta") {
        m <- stats::dfbeta(fit)
        return(2 * max(sqrt(rowSums(m^2))))
    }
    # leverage bound
    h <- stats::hatvalues(fit)
    r <- stats::residuals(fit)
    X <- stats::model.matrix(fit)
    XtXi <- solve(crossprod(X))
    2 * max(sqrt(rowSums((X %*% XtXi)^2)) * abs(r) / (1 - h))
}

# ---------------------------------------------------------------------------
# Global sensitivities of the COEFFICIENT channel (public X)
# ---------------------------------------------------------------------------
# These are the two functions the package was missing, and their absence is why
# every generator here fell back on the DFBETA local sensitivity above.  Both
# are functions of the public design and a public bound on |y|, and of nothing
# else, so they hold across every neighbouring pair.  That is the whole
# difference between a bound and a statistic.

#' Global sensitivity of the OLS coefficient vector (public X)
#'
#' \eqn{\hat\beta = (X'X)^{-1} X'y}.  Replacing row i's response moves
#' \eqn{X'y} by \eqn{x_i (y_i' - y_i)} and nothing else, so
#' \deqn{\hat\beta' - \hat\beta = (X'X)^{-1} x_i (y_i' - y_i),}
#' whose norm is at most \eqn{2 B_y \|(X'X)^{-1} x_i\|_2} over the bounded
#' universe.  Hence
#' \deqn{\Delta_\beta = 2 B_y \max_i \|(X'X)^{-1} x_i\|_2,}
#' a function of the public design and the public bound alone.
#'
#' X must have full column rank, without which \eqn{(X'X)^{-1}} does not exist
#' and \eqn{\Delta_\beta} is undefined.  Since \eqn{\Delta_\beta} is the
#' sensitivity the noise is calibrated to, an undefined value is an undefined
#' noise scale rather than a loose one, so this function lets \code{solve()}
#' error rather than regularising.  The condition is a hypothesis of the
#' public-covariate exact release proposition (added 2026-07-31, review
#').
#'
#' Ported from \code{gs.beta.ols.public} in \file{lib/28-dp-projection.r}.
#'
#' @param m.X Design (model) matrix, n x p.  public by the route's premise, and
#'   of full column rank.
#' @param B.y Public bound on |y|.
#' @return Positive scalar L2 sensitivity of beta_hat.
#' @noRd
.gs.beta.ols.public <- function(m.X, B.y) {
    m.inv <- solve(crossprod(m.X))
    2 * B.y * max(sqrt(rowSums((m.X %*% m.inv)^2)))
}

#' Global sensitivity of the ridge coefficient vector (public X)
#'
#' With \eqn{A = X'X + \lambda I}, \eqn{\hat\beta_\mathrm{ridge} = A^{-1} X'y}.
#' The same one-row argument gives
#' \deqn{\Delta_\beta = 2 B_y \max_i \|A^{-1} x_i\|_2,}
#' which depends on X, lambda and B.y only.  It is
#' \code{.gs.beta.ols.public} with A in place of X'X, and reduces to it
#' exactly as lambda tends to 0.
#'
#' Ported from \code{gs.beta.ridge.public} in
#' \file{lib/34-dp-ridge-projection.r}, the fix.  V3 was left on the ridge
#' DFBETA local sensitivity when V2 was migrated on 2026-07-07, because the
#' earlier repair touched only its variance channel and left the coefficient
#' channel behind.  V3 had no excuse for the local route: its design is public,
#' so this global bound was always available, and it is two lines.
#'
#' @param m.X Design (model) matrix, n x p.  public by the route's premise.
#' @param lambda Ridge penalty (positive).
#' @param B.y Public bound on |y|.
#' @return Positive scalar L2 sensitivity of beta_ridge.
#' @noRd
.gs.beta.ridge.public <- function(m.X, lambda, B.y) {
    m.A.inv <- solve(crossprod(m.X) + lambda * diag(ncol(m.X)))
    2 * B.y * max(sqrt(rowSums((m.X %*% m.A.inv)^2)))
}

#' Worst-case sensitivity of the OLS variance estimate (public X)
#'
#' One-row replacement of y (with |y_i| <= B.y) changes
#' sigma_hat^2 = y'(I - H)y / (n - p) by at most 4 B.y^2 K_X / (n - p),
#' where K_X = max_i max(1 - h_ii, sum_{j != i} |h_ij|) is the hat-row
#' constant of the public design. The exact identity
#' RSS(y') - RSS(y) = 2 d r_i + d^2 (1 - h_ii) makes the bound attainable
#' over the bounded universe at an adversarial sign pattern, so K_X cannot
#' be dropped: the naive 4 B.y^2 / (n - p) under-states the worst case
#' whenever K_X > 1, which real designs exceed comfortably.
#'
#' @param m.X Design (model) matrix, n x p.
#' @param B.y Public bound on |y|.
#' @return Positive scalar sensitivity of sigma_hat^2.
#' @noRd
.gs.sigma2.public <- function(m.X, B.y) {
    n <- nrow(m.X)
    p <- ncol(m.X)
    m.H <- m.X %*% solve(crossprod(m.X)) %*% t(m.X)
    v.h <- diag(m.H)
    v.s <- rowSums(abs(m.H)) - abs(v.h)
    4 * B.y^2 * max(pmax(1 - v.h, v.s)) / (n - p)
}

#' Worst-case sensitivity of the ridge variance estimate (public X)
#'
#' The ridge RSS is y'My with M = (I - H_l)'(I - H_l) and
#' H_l = X (X'X + lambda I)^{-1} X'. With A = X'X + lambda I this expands
#' to M = I - X A^{-1} (X'X + 2 lambda I) A^{-1} X', avoiding the n x n
#' product. The OLS one-row identity holds with M in place of I - H,
#' giving |Delta sigma_hat^2| <= 4 B.y^2 K_M / (n - p) with
#' K_M = max_i max(M_ii, sum_{j != i} |M_ij|). Reduces to
#' .gs.sigma2.public() at lambda = 0.
#'
#' @param m.X Design (model) matrix, n x p.
#' @param lambda Ridge penalty (positive).
#' @param B.y Public bound on |y|.
#' @return Positive scalar sensitivity of sigma_hat^2.
#' @noRd
.gs.sigma2.ridge <- function(m.X, lambda, B.y) {
    n <- nrow(m.X)
    p <- ncol(m.X)
    m.XtX <- crossprod(m.X)
    m.A.inv <- solve(m.XtX + lambda * diag(p))
    m.M <- diag(n) - m.X %*% m.A.inv %*% (m.XtX + 2 * lambda * diag(p)) %*%
        m.A.inv %*% t(m.X)
    v.m <- diag(m.M)
    v.s <- rowSums(abs(m.M)) - abs(v.m)
    4 * B.y^2 * max(pmax(v.m, v.s)) / (n - p)
}
