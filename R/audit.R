# audit.R -- Empirical differential-privacy audit for synthexpln mechanisms.
# Public functions: dp.audit().

# +------------------------------------------------------------------------+
# | dp.audit -- Wasserman-Zhou membership-inference DP audit               |
# +------------------------------------------------------------------------+

#' Empirical differential-privacy audit via membership inference
#'
#' Constructs a pair of neighbouring datasets \eqn{D} and \eqn{D'} (differing
#' in one high-leverage row), generates \code{n.trials} mechanism outputs from
#' each, and measures whether an optimal linear classifier can distinguish the
#' two output distributions better than the Wasserman-Zhou theoretical bound
#' \deqn{(e^\varepsilon + \delta) / (1 + e^\varepsilon + \delta).}
#' If the empirical AUC exceeds this bound by more than \code{tolerance}, the
#' audit fails (\code{pass = FALSE}), suggesting the claimed
#' \eqn{(\varepsilon, \delta)} is under-stated for this mechanism.
#'
#' The neighbouring datasets are constructed by replacing the row with the
#' largest DFBETA norm (highest leverage under \code{lm(y ~ .)}) with the
#' row with the smallest DFBETA norm.  The first column of \code{data} is
#' treated as the response \eqn{y}.
#'
#' @param generator A function \code{function(data)} that accepts a data frame
#'   and returns a \code{synthexpln} object (or any list with a \code{$beta.dp}
#'   or \code{$beta} slot).
#' @param data A data frame.  Its first column is treated as the response.
#' @param epsilon Numeric scalar; the claimed privacy parameter
#'   \eqn{\varepsilon > 0}.
#' @param delta Numeric scalar; the claimed failure probability
#'   \eqn{\delta \ge 0}.
#' @param n.trials Integer; number of mechanism outputs to generate from each
#'   dataset.  Default 1000.
#' @param alpha Numeric scalar; significance level (currently unused;
#'   reserved for future hypothesis-test extensions).  Default 0.05.
#' @param tolerance Numeric scalar; slack added to the theoretical bound
#'   before declaring failure (accounts for Monte Carlo error).  Default 0.02.
#' @return A named list:
#'   \describe{
#'     \item{\code{pass}}{Logical; \code{TRUE} iff
#'       \code{auc <= theory.bound + tolerance}.}
#'     \item{\code{auc}}{Numeric; empirical AUC of the projection classifier.}
#'     \item{\code{theory.bound}}{Numeric; Wasserman-Zhou theoretical bound.}
#'     \item{\code{n.trials}}{Integer; requested number of trials.}
#'     \item{\code{epsilon}}{Numeric; \eqn{\varepsilon} as supplied.}
#'     \item{\code{delta}}{Numeric; \eqn{\delta} as supplied.}
#'   }
#'   If fewer than 10 trials succeed per dataset, returns
#'   \code{list(pass = NA, message = "insufficient successful trials")} instead.
#' @references
#'   Wasserman, L. and Zhou, S. (2010). A statistical framework for
#'   differential privacy. \emph{JASA}, 105(489), 375-389.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = rnorm(200), x = rnorm(200))
#' gen <- function(data) gen.syn.dp.projected(data, y ~ x,
#'                                             epsilon = 2, delta = 1e-6)
#' dp.audit(gen, d, epsilon = 2, delta = 1e-6, n.trials = 100)
#' }
dp.audit <- function(generator, data, epsilon, delta,
                     n.trials = 1000L, alpha = 0.05, tolerance = 0.02) {

    # ------------------------------------------------------------------
    # Build neighbouring datasets D and D'
    # ------------------------------------------------------------------
    y.nm  <- names(data)[1L]
    fml   <- stats::reformulate(setdiff(names(data), y.nm), response = y.nm)
    mo    <- stats::lm(fml, data = data)
    v.dfb <- sqrt(rowSums(stats::dfbeta(mo)^2))

    i.high  <- which.max(v.dfb)
    i.other <- which.min(v.dfb)

    d.prime           <- data
    d.prime[i.high, ] <- data[i.other, ]

    # ------------------------------------------------------------------
    # Helper: run generator n.trials times, extracting beta coefficients
    # ------------------------------------------------------------------
    .run.trials <- function(dataset) {
        l.beta <- Map(function(t) {
            res <- tryCatch(generator(dataset), error = function(e) NULL)
            if (is.null(res)) return(NULL)
            # Prefer beta.dp; fall back to beta.
            b <- if (!is.null(res[["beta.dp"]])) res[["beta.dp"]]
                 else res[["beta"]]
            if (!is.null(b)) as.numeric(b) else NULL
        }, seq_len(n.trials))
        Filter(Negate(is.null), l.beta)
    }

    l.D  <- .run.trials(data)
    l.Dp <- .run.trials(d.prime)

    # ------------------------------------------------------------------
    # Guard against too few successful trials
    # ------------------------------------------------------------------
    if (length(l.D) < 10L || length(l.Dp) < 10L)
        return(list(pass = NA,
                    message = "insufficient successful trials"))

    # ------------------------------------------------------------------
    # Empirical AUC via difference-of-means projection
    # (Wasserman & Zhou 2010; source 50-dp-audit.r lines 111-118)
    # ------------------------------------------------------------------
    m.D  <- do.call(rbind, l.D)
    m.Dp <- do.call(rbind, l.Dp)

    v.mean.diff <- colMeans(m.D) - colMeans(m.Dp)

    # If all columns of v.mean.diff are 0 (identical distributions), scores
    # are all 0 and AUC = 0.5 by the tie-correction term.
    scores <- rbind(m.D, m.Dp) %*% v.mean.diff
    labels <- c(rep(1L, nrow(m.D)), rep(0L, nrow(m.Dp)))

    v.pos <- scores[labels == 1L]
    v.neg <- scores[labels == 0L]
    s.auc <- mean(outer(v.pos, v.neg, ">") + 0.5 * outer(v.pos, v.neg, "=="))

    # ------------------------------------------------------------------
    # Theoretical bound and verdict
    # ------------------------------------------------------------------
    s.theory <- (exp(epsilon) + delta) / (1 + exp(epsilon) + delta)
    s.pass   <- s.auc <= s.theory + tolerance

    list(
        pass         = s.pass,
        auc          = s.auc,
        theory.bound = s.theory,
        n.trials     = n.trials,
        epsilon      = epsilon,
        delta        = delta
    )
}
