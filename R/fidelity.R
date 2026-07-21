# fidelity.R -- Fidelity metrics for synthexpln synthetic data.
# Public functions: fidelity.inferential(), fidelity.distributional(),
#                   calibrated.se().

# +------------------------------------------------------------------------+
# | calibrated.se -- DP-corrected standard error                          |
# +------------------------------------------------------------------------+

#' Calibrated standard error for DP projection
#'
#' Corrects synthetic-fit standard errors to account for the additional
#' uncertainty introduced by the Gaussian DP mechanism on the regression
#' coefficients.  The correction combines the analytic-fit variance
#' \eqn{\hat{\sigma}^2_j} with the noise variance
#' \eqn{(\sigma_\beta / s_j)^2} in quadrature.
#'
#' @param se.syn Numeric vector of synthetic-fit standard errors (e.g.,
#'   column 2 of \code{summary(lm(...))$coefficients}).
#' @param sigma.beta Scalar; Gaussian mechanism sd applied to the
#'   standardised beta (i.e., \code{obj$sigma.beta} from
#'   \code{\link{gen.syn.dp.projected}}).
#' @param s.j Numeric vector of column scales matching \code{se.syn}, needed
#'   \strong{only} when the DP noise was added on a column-standardised scale.
#'   Defaults to 1, which is the correct value whenever \code{sigma.beta} is
#'   already the noise sd on the raw coefficient scale.  Prefer passing
#'   \code{obj$sigma.beta.raw} and leaving this alone.
#' @return Numeric vector; element-wise
#'   \eqn{\sqrt{\hat{\sigma}_j^2 + (\sigma_\beta / s_j)^2}}.
#' @section Do not divide by a column scale that is not there:
#' The correction is \eqn{\sqrt{SE_\mathrm{syn}^2 + \sigma_{\beta,j}^2}}, where
#' \eqn{\sigma_{\beta,j}} is the sd of the DP noise \strong{on the raw scale of
#' coefficient j}.  There is no column scale in it.
#'
#' A \code{s.j} divisor is right only when the mechanism perturbed the
#' \emph{standardised} coefficients, since noise of sd \eqn{\sigma_\beta} on
#' \eqn{\beta_j s_j} lands as sd \eqn{\sigma_\beta / s_j} on \eqn{\beta_j}.  The
#' sound routes do not do that: they add isotropic noise to the raw
#' \eqn{\hat\beta} at a global sensitivity, and carry no column scale at all.
#' Dividing anyway shrinks the DP variance component by the column scale, which
#' on a wide-scale covariate (an income in the thousands, say) can be a factor of
#' several hundred, and the reported interval is then barely wider than the
#' uncorrected one.  That understates the uncertainty of a private release, which
#' is the only direction an error here must not take.
#'
#' Every DP generator reports \code{$sigma.beta.raw}, a \strong{named} vector of
#' raw-scale noise sds, one per coefficient.  Use it.
#' @export
#' @examples
#' calibrated.se(c(0.1, 0.2), sigma.beta = 0.05)
calibrated.se <- function(se.syn, sigma.beta, s.j = 1) {
    sqrt(se.syn^2 + (sigma.beta / s.j)^2)
}

# +------------------------------------------------------------------------+
# | .dp.noise.sd -- the raw-scale noise sd, per coefficient, or an error    |
# +------------------------------------------------------------------------+
# Take the noise sd DIRECTLY from the object and hard-error if it does not name
# every coefficient.  Do not infer it, and do not silently proceed without it.
#
# The old code asked for $sigma.beta and $s.j and, if either was missing, quietly
# returned se.corr = se.syn "because there is no DP noise to correct for".  Once
# the sound routes stopped carrying a column scale, that branch began firing on
# releases that do carry DP noise, and it reported an uncorrected interval for a
# private release.  Silence is the wrong failure mode here: an SE that is too
# small is an inferential claim the release cannot support.

#' Raw-scale DP noise sd per coefficient (internal)
#' @noRd
.dp.noise.sd <- function(obj, v.names) {
    v.raw <- obj[["sigma.beta.raw"]]
    if (is.null(v.raw)) return(NULL)          # not a DP object at all
    if (is.null(names(v.raw)) || !all(v.names %in% names(v.raw)))
        stop("the DP object's $sigma.beta.raw does not name every coefficient ",
             "being calibrated (missing: ",
             paste(setdiff(v.names, names(v.raw)), collapse = ", "),
             "). Refusing to guess a noise scale: an SE that is too small is an ",
             "inferential claim the release cannot support.", call. = FALSE)
    v.raw[v.names]
}

# +------------------------------------------------------------------------+
# | fidelity.inferential -- coefficient-level fidelity metrics            |
# +------------------------------------------------------------------------+

#' Inferential fidelity metrics for a synthexpln object
#'
#' Fits the same regression model on both the original and synthetic data,
#' then compares coefficient vectors on four criteria: maximum absolute
#' percentage deviation, sign agreement, significance agreement, and
#' (optionally) calibrated standard errors.
#'
#' For the non-DP \code{\link{projection}} output (no \code{$sigma.beta.raw}
#' slot), the \code{calibrate} argument is accepted but returns
#' \code{se.corr = se.syn} with a message because no DP noise was added.
#'
#' @param obj A \code{synthexpln} object, typically from
#'   \code{\link{projection}} or \code{\link{gen.syn.dp.projected}}.
#' @param orig Data frame; the original (non-synthetic) data.
#' @param formula A model formula.  Required.
#' @param alpha Significance level for p-value thresholding.
#'   Default \code{0.05}.
#' @param calibrate Logical; if \code{TRUE} compute calibrated SEs.  Requires
#'   \code{obj$sigma.beta.raw} alone, the named vector of raw-scale noise sds
#'   that every DP generator reports.  It needs no column scale, and asking for
#'   one was the defect: once the sound routes stopped carrying \code{$s.j},
#'   a requirement keyed on it fired on releases that do carry DP noise and
#'   returned an uncorrected interval for a private release.  See the
#'   \emph{Do not divide by a column scale that is not there} section of
#'   \code{\link{calibrated.se}}.  Default \code{FALSE}.
#' @return A list with components:
#'   \describe{
#'     \item{\code{max.abs.pct}}{Maximum absolute percentage deviation of
#'       synthetic vs original coefficients.}
#'     \item{\code{sign.agreement}}{Fraction of coefficients with matching
#'       sign (in \eqn{[0,1]}).}
#'     \item{\code{sig.agreement}}{Fraction of coefficients with matching
#'       significance at level \code{alpha} (in \eqn{[0,1]}).}
#'     \item{\code{coef.orig}}{Named vector of original coefficients.}
#'     \item{\code{coef.syn}}{Named vector of synthetic coefficients.}
#'     \item{\code{pct.dev}}{Percentage deviations per coefficient.}
#'     \item{\code{p.orig}}{p-values from the original model.}
#'     \item{\code{p.syn}}{p-values from the synthetic model.}
#'     \item{\code{se.syn}}{(Only when \code{calibrate = TRUE}) Synthetic
#'       model standard errors before correction.}
#'     \item{\code{se.corr}}{(Only when \code{calibrate = TRUE}) Calibrated
#'       standard errors after DP noise adjustment.}
#'   }
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(100), x = rnorm(100))
#' syn <- projection(d, y ~ x)
#' fidelity.inferential(syn, orig = d, formula = y ~ x)
fidelity.inferential <- function(obj, orig, formula = NULL,
                                  alpha = 0.05, calibrate = FALSE) {
    if (is.null(formula))
        stop("'formula' is required; cannot be NULL")

    # select model family: binomial glm for V5, OLS for all others
    b.logistic <- identical(obj[["variant"]], "V5")

    .fit.model <- function(data) {
        if (b.logistic) {
            stats::glm(formula, data = data,
                       family = stats::binomial())
        } else {
            stats::lm(formula, data = data)
        }
    }

    mo.orig <- .fit.model(orig)
    mo.syn  <- .fit.model(obj[["syn"]])

    # align common coefficient names
    v.common <- intersect(names(stats::coef(mo.orig)),
                          names(stats::coef(mo.syn)))

    coef.orig <- stats::coef(mo.orig)[v.common]
    coef.syn  <- stats::coef(mo.syn)[v.common]

    # percentage deviation with divide-by-zero guard
    pct.dev <- mapply(function(co, cs) {
        if (abs(co) < 1e-10) {
            abs(cs - co) * 100
        } else {
            100 * (cs - co) / abs(co)
        }
    }, coef.orig, coef.syn)
    pct.dev <- `names<-`(pct.dev, v.common)

    max.abs.pct    <- max(abs(pct.dev), na.rm = TRUE)
    sign.agreement <- mean(sign(coef.syn) == sign(coef.orig))

    # p-value column name differs between lm (Pr(>|t|)) and glm (Pr(>|z|))
    p.col <- if (b.logistic) "Pr(>|z|)" else "Pr(>|t|)"
    p.orig <- summary(mo.orig)[["coefficients"]][v.common, p.col]
    p.syn  <- summary(mo.syn)[["coefficients"]][v.common, p.col]

    sig.orig      <- p.orig < alpha
    sig.syn       <- p.syn  < alpha
    sig.agreement <- mean(sig.orig == sig.syn)

    out <- list(
        max.abs.pct    = max.abs.pct,
        sign.agreement = sign.agreement,
        sig.agreement  = sig.agreement,
        coef.orig      = coef.orig,
        coef.syn       = coef.syn,
        pct.dev        = pct.dev,
        p.orig         = p.orig,
        p.syn          = p.syn
    )

    if (calibrate) {
        se.syn <- summary(mo.syn)[["coefficients"]][v.common, 2]
        # Take the raw-scale per-coefficient noise sd straight off the object.
        # .dp.noise.sd() returns NULL only when there is genuinely no DP noise
        # (a plain projection), and ERRORS rather than guessing when a DP object
        # cannot name its own noise scale.
        v.sd <- .dp.noise.sd(obj, v.common)
        if (is.null(v.sd)) {
            message("calibrate = TRUE but obj carries no DP noise ",
                    "(no $sigma.beta.raw); returning se.corr = se.syn")
            se.corr <- se.syn
        } else {
            # No column scale: the correction is sqrt(SE_syn^2 + sigma_beta_j^2).
            se.corr <- calibrated.se(se.syn, v.sd)
        }
        out <- c(out, list(se.syn = se.syn, se.corr = se.corr))
    }

    out
}

# +------------------------------------------------------------------------+
# | fidelity.distributional -- marginal distributional fidelity           |
# +------------------------------------------------------------------------+

#' Distributional fidelity metrics per column
#'
#' Measures marginal distributional similarity between synthetic and original
#' data using the Kolmogorov-Smirnov statistic (for numeric columns) or
#' total-variation distance (for factor/character columns).
#'
#' @param syn A data frame of synthetic data, or a \code{synthexpln} object
#'   (in which case \code{syn$syn} is used automatically).
#' @param orig A data frame of original data.
#' @param type Metric to use.  Currently only \code{"ks"} is implemented.
#'   \code{"energy"} and \code{"mmd"} are reserved for future releases.
#' @return A named numeric vector with one entry per column that is present
#'   in both \code{syn} and \code{orig}:
#'   \itemize{
#'     \item For numeric columns: KS statistic in \eqn{[0, 1]}.
#'     \item For factor/character columns: total-variation distance
#'       \eqn{(0.5) \sum |p_i - q_i|} in \eqn{[0, 0.5]}.
#'   }
#' @note \code{type = "energy"} (energy distance) and \code{type = "mmd"}
#'   (maximum mean discrepancy) are not yet implemented.  Contributions are
#'   welcome; for now use \code{type = "ks"}.
#' @export
#' @examples
#' set.seed(1)
#' orig <- data.frame(a = rnorm(100), b = rnorm(100))
#' syn  <- data.frame(a = rnorm(100), b = rnorm(100))
#' fidelity.distributional(syn, orig, type = "ks")
fidelity.distributional <- function(syn, orig,
                                     type = c("ks", "energy", "mmd")) {
    type <- match.arg(type)

    if (type != "ks")
        stop("not implemented; use type = 'ks' for now")

    # accept synthexpln objects transparently
    if (inherits(syn, "synthexpln"))
        syn <- syn[["syn"]]

    v.cols <- intersect(names(orig), names(syn))

    .col.stat <- function(col) {
        v.o <- orig[[col]]
        v.s <- syn[[col]]
        if (is.numeric(v.o)) {
            suppressWarnings(
                stats::ks.test(v.o, v.s)[["statistic"]])
        } else {
            # total-variation distance for factor/character
            v.lvls <- union(as.character(v.o), as.character(v.s))
            t.o <- prop.table(table(factor(v.o, levels = v.lvls)))
            t.s <- prop.table(table(factor(v.s, levels = v.lvls)))
            0.5 * sum(abs(as.numeric(t.o) - as.numeric(t.s)))
        }
    }

    v.stats <- mapply(.col.stat, v.cols)
    `names<-`(as.numeric(v.stats), v.cols)
}
