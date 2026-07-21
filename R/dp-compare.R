# dp-compare.R -- one-call comparison of a synthetic release against the
# original data: inferential fidelity, marginal fidelity, dependence fidelity
# (including the categorical-continuous associations the mixed-margin copula
# preserves), and the privacy accounting the release carries.
# Public functions: dp.compare(), print.dp.comparison()

#' Compare a synthetic release against the original data
#'
#' One call binding the package's fidelity metrics and the release's privacy
#' accounting into a single structured report:
#' \code{\link{fidelity.inferential}} (coefficient recovery, sign and
#' significance agreement, calibrated standard errors),
#' \code{\link{fidelity.distributional}} (per-column marginal distances), a
#' dependence section (continuous-continuous correlation error and the
#' categorical-continuous correlation ratios \eqn{\eta} that the mixed-margin
#' copula preserves), and the privacy accounting read off the release object.
#'
#' @section What this does and does not quantify:
#' Differential privacy is a property of the mechanism, quantified over all
#' neighbouring datasets, not a property of one realised output.  No
#' comparison of a synthetic dataset against the original can certify or
#' measure it.  What this report gives is (a) utility, how close the release
#' is to the original on the inferential, marginal, and dependence metrics,
#' and (b) the accounting the release carries, echoed from the mechanism: the
#' claimed \eqn{(\varepsilon, \delta)}, the composed \eqn{\mu}, and the
#' membership-inference ceiling \eqn{\Phi(\mu / \sqrt{2})}, the AUC that no
#' attack on a \eqn{\mu}-GDP mechanism can exceed.  An empirical audit
#' against that ceiling, which attacks the mechanism on neighbouring pairs
#' and catches implementation errors the composition cannot see, is
#' \code{\link{dp.audit}}; empirical attacks lower-bound leakage and not
#' certify its absence.  Releases produced by the retired
#' \code{sens.method = "local"} path are reported as carrying no guarantee.
#'
#' @param obj A \code{synthexpln} release object (from any generator).
#' @param orig The original data frame.
#' @param formula Model formula for the inferential comparison.  Required.
#' @param alpha Significance level for the agreement metrics.  Default
#'   \code{0.05}.
#' @param calibrate Logical; compute calibrated standard errors where the
#'   object carries a per-coefficient noise scale.  Default \code{TRUE}.
#' @return An object of class \code{dp.comparison}: a list with
#'   \describe{
#'     \item{\code{$inferential}}{The \code{\link{fidelity.inferential}}
#'       result.}
#'     \item{\code{$marginals}}{The \code{\link{fidelity.distributional}}
#'       vector (KS for numeric columns, total-variation for categorical).}
#'     \item{\code{$dependence}}{A list with \code{cor.frob} and
#'       \code{cor.max.abs.err} over the continuous block, and \code{eta}, a
#'       data frame of categorical-continuous correlation ratios (original
#'       against synthetic, one row per pair).}
#'     \item{\code{$privacy}}{The accounting echoed from the release:
#'       \code{guarantee} (a one-line statement), \code{epsilon},
#'       \code{delta}, \code{mu.total}, \code{mu.channels}, and
#'       \code{mi.ceiling}, the membership-inference AUC bound
#'       \eqn{\Phi(\mu / \sqrt 2)}.}
#'   }
#' @seealso \code{\link{dp.audit}} for the empirical membership-inference
#'   audit of a mechanism, \code{\link{privacy.metrics}} for per-record
#'   disclosure heuristics.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(200), age = rnorm(200, 50, 10),
#'                 sex = factor(sample(c("M", "F"), 200, replace = TRUE)))
#' x.spec <- list(age = list(type = "continuous", bounds = c(0, 100)),
#'                sex = list(type = "categorical", levels = c("M", "F")))
#' # The boxes are facts about this example's DGP (coefficients 0, sigma^2 1),
#' # fixed before the data are drawn.
#' syn <- gen.syn.dp.full(d, y ~ age + sex, x.spec = x.spec,
#'                         epsilon = 5, delta = 1e-6,
#'                         clip.lo = -1, clip.hi = 1,
#'                         disp.lo = 0, disp.hi = 4, seed = 1)
#' cmp <- dp.compare(syn, d, y ~ age + sex)
#' cmp$privacy$mi.ceiling
#' print(cmp)
dp.compare <- function(obj, orig, formula = NULL, alpha = 0.05,
                       calibrate = TRUE) {
    if (!inherits(obj, "synthexpln"))
        stop("obj must be a synthexpln release object", call. = FALSE)
    if (is.null(formula))
        stop("'formula' is required; cannot be NULL", call. = FALSE)

    l.inf  <- fidelity.inferential(obj, orig, formula, alpha = alpha,
                                   calibrate = calibrate)
    v.marg <- fidelity.distributional(obj, orig)
    l.dep  <- .dependence.fidelity(obj[["syn"]], orig)

    # ---- the accounting, echoed from the mechanism ---------------------------
    s.mu <- obj[["mu.total"]]
    s.guarantee <- if (identical(obj[["sens.method"]], "local")) {
        paste0("NONE. sens.method = \"local\" calibrates to a data-dependent ",
               "sensitivity, so no (epsilon, delta) statement holds (retired).")
    } else if (!is.null(s.mu) && is.finite(s.mu)) {
        sprintf(paste0("mu-GDP at mu = %.4f, the exact dual of the supplied ",
                       "(epsilon = %g, delta = %g)."),
                s.mu, obj[["epsilon"]], obj[["delta"]])
    } else if (!is.null(obj[["epsilon"]])) {
        sprintf("(epsilon = %g, delta = %g), as supplied by the mechanism.",
                obj[["epsilon"]], obj[["delta"]])
    } else {
        "Not a DP release. No privacy mechanism ran."
    }
    l.priv <- list(
        guarantee   = s.guarantee,
        epsilon     = obj[["epsilon"]],
        delta       = obj[["delta"]],
        mu.total    = s.mu,
        mu.channels = obj[["mu.channels"]],
        mi.ceiling  = if (!is.null(s.mu) && is.finite(s.mu))
                          stats::pnorm(s.mu / sqrt(2)) else NA_real_)

    structure(list(inferential = l.inf, marginals = v.marg,
                   dependence = l.dep, privacy = l.priv),
              class = "dp.comparison")
}

#' Dependence fidelity between two data frames
#'
#' Continuous-continuous correlation error over the shared numeric columns,
#' and the categorical-continuous correlation ratio \eqn{\eta} per
#' factor-numeric pair (original against synthetic).  \eqn{\eta} is the
#' square root of the between-level share of the numeric column's variance;
#' for a binary factor it is the absolute point-biserial correlation.  Under
#' independence \eqn{\eta} concentrates near \eqn{\sqrt{(K-1)/(n-1)}} rather
#' than at zero.
#'
#' @param syn,orig Data frames sharing column names.
#' @return list(cor.frob, cor.max.abs.err, eta = data.frame).
#' @noRd
.dependence.fidelity <- function(syn, orig) {
    v.cols <- intersect(names(orig), names(syn))
    v.num  <- Filter(function(nm) is.numeric(orig[[nm]]), v.cols)
    v.cat  <- setdiff(v.cols, v.num)

    l.cor <- if (length(v.num) >= 2L) {
        m.d <- stats::cor(syn[v.num]) - stats::cor(orig[v.num])
        list(cor.frob = sqrt(sum(m.d^2)), cor.max.abs.err = max(abs(m.d)))
    } else list(cor.frob = NA_real_, cor.max.abs.err = NA_real_)

    eta.ratio <- function(v.g, v.y) {
        s.gm <- mean(v.y)
        l.g  <- Filter(function(g) length(g) > 0L, split(v.y, v.g))
        ss.b <- sum(vapply(l.g, function(g) length(g) * (mean(g) - s.gm)^2,
                           numeric(1)))
        ss.t <- sum((v.y - s.gm)^2)
        if (ss.t <= 0) NA_real_ else sqrt(ss.b / ss.t)
    }
    d.eta <- if (length(v.cat) > 0L && length(v.num) > 0L) {
        d.grid <- expand.grid(cat = v.cat, cont = v.num,
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        d.grid[["eta.orig"]] <- mapply(function(g, y)
            eta.ratio(orig[[g]], orig[[y]]), d.grid[["cat"]], d.grid[["cont"]])
        d.grid[["eta.syn"]] <- mapply(function(g, y)
            eta.ratio(syn[[g]], syn[[y]]), d.grid[["cat"]], d.grid[["cont"]])
        d.grid
    } else data.frame(cat = character(0), cont = character(0),
                      eta.orig = numeric(0), eta.syn = numeric(0))

    c(l.cor, list(eta = d.eta))
}

#' Print a dp.comparison object
#'
#' @param x A \code{dp.comparison} object from \code{\link{dp.compare}}.
#' @param ... Unused, for S3 compatibility.
#' @return \code{x}, invisibly.  Called for its printed summary.
#' @export
print.dp.comparison <- function(x, ...) {
    cat("synthetic-vs-original comparison (synthexpln)\n")
    cat(strrep("-", 62), "\n")
    l.i <- x[["inferential"]]
    cat("Inferential fidelity\n")
    cat(sprintf("  max |%% coef deviation|   %.3f%%\n", l.i[["max.abs.pct"]]))
    cat(sprintf("  sign agreement           %.0f%%\n",
                100 * l.i[["sign.agreement"]]))
    cat(sprintf("  significance agreement   %.0f%%\n",
                100 * l.i[["sig.agreement"]]))
    cat("Marginal fidelity (KS numeric, TV categorical)\n")
    invisible(Map(function(nm, v) cat(sprintf("  %-22s %.4f\n", nm, v)),
                  names(x[["marginals"]]), x[["marginals"]]))
    l.d <- x[["dependence"]]
    cat("Dependence fidelity\n")
    if (is.finite(l.d[["cor.frob"]]))
        cat(sprintf("  continuous cor error     Frobenius %.4f, max %.4f\n",
                    l.d[["cor.frob"]], l.d[["cor.max.abs.err"]]))
    d.eta <- l.d[["eta"]]
    invisible(Map(function(i) cat(sprintf(
        "  eta(%s, %s)%s original %.4f, synthetic %.4f\n",
        d.eta[i, "cat"], d.eta[i, "cont"],
        strrep(" ", max(1, 10 - nchar(d.eta[i, "cat"]) -
                           nchar(d.eta[i, "cont"]))),
        d.eta[i, "eta.orig"], d.eta[i, "eta.syn"])),
        seq_len(nrow(d.eta))))
    l.p <- x[["privacy"]]
    cat("Privacy accounting (a property of the MECHANISM, echoed here;\n")
    cat("no dataset comparison can certify it. See ?dp.audit.)\n")
    cat(sprintf("  guarantee   %s\n", l.p[["guarantee"]]))
    if (is.finite(l.p[["mi.ceiling"]]))
        cat(sprintf("  MI ceiling  AUC <= %.3f (= Phi(mu / sqrt(2)))\n",
                    l.p[["mi.ceiling"]]))
    invisible(x)
}
