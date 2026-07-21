# The public block-count rule ------------------------------------------------
# Public functions:  dp.rate.public(), dp.subagg.blocks()
#
# m = floor(n_eff / (5 p)): the largest block count leaving five rows per
# released coefficient per block.  n is public under the replace-one
# neighbouring convention, p is the number of released coefficients, and on a
# two-part positive channel n_eff scales n by the expected positive share
# n * (1 - pi0), with pi0 the published zero rate.  A realised rate
# (mean(y == 0)) would calibrate the mechanism to the very data it protects,
# so the rate travels as a provenance object on the dp.dispersion pattern.

#' Declare a public zero rate
#'
#' Wraps a \strong{published} zero rate \eqn{\pi_0} for the zero-inflated
#' generator's block-count rule, together with a required \code{source} string
#' naming where it comes from.  The positive channel of
#' \code{\link{gen.syn.dp.projected.zi}} fits the positive rows only, so its
#' block count under the rule is \eqn{\lfloor n (1 - \pi_0) / (5 p) \rfloor},
#' and \eqn{\pi_0} must be a design constant (prior literature, a published
#' registry rate, the data owner's published table).  A rate computed from the
#' sample, such as \code{mean(y == 0)}, makes the block count a function of the
#' private data and voids the calibration.
#'
#' @param pi0 The published zero rate, a single value in \eqn{[0, 1)}.
#' @param source Character string naming where the rate comes from.  Required,
#'   and printed alongside the release so the provenance is on the record.
#' @return An object of class \code{dp.rate}.
#' @seealso \code{\link{dp.subagg.blocks}},
#'   \code{\link{gen.syn.dp.projected.zi}}.
#' @export
#' @examples
#' dp.rate.public(0.6, source = "published registry share of zero-cost users")
dp.rate.public <- function(pi0, source) {
    if (missing(source) || !is.character(source) || length(source) != 1L ||
        !nzchar(source))
        stop("source is REQUIRED: name where this PUBLIC zero rate comes from ",
             "(a published design constant, not the sample).  If it came ",
             "from the sample, it is not public.", call. = FALSE)
    if (!is.numeric(pi0) || length(pi0) != 1L || is.na(pi0) || pi0 < 0 ||
        pi0 >= 1)
        stop("pi0 must be a single value in [0, 1)", call. = FALSE)
    structure(list(pi0 = pi0, source = source), class = "dp.rate")
}

#' @export
print.dp.rate <- function(x, ...) {
    cat(sprintf("<dp.rate: pi0 = %g>\n  source: %s\n",
                x[["pi0"]], x[["source"]]))
    invisible(x)
}

#' The public block-count rule for subsample-and-aggregate
#'
#' Returns \eqn{m = \lfloor n_{\mathrm{eff}} / (5 p) \rfloor}, the largest
#' block count leaving five rows per released coefficient per block.  This is
#' the default block count of every subsample-and-aggregate generator in the
#' package.  Its inputs are public: \code{n} is fixed under the replace-one
#' neighbouring convention, \code{p} is the number of released coefficients,
#' and \code{rate} (when given) scales \code{n} to the expected positive count
#' \eqn{n (1 - \pi_0)} using a \strong{published} zero rate declared through
#' \code{\link{dp.rate.public}}.  The rule is an a priori cap, not a
#' quantity tuned on the analysis sample.
#'
#' @param n Number of rows the partition covers.
#' @param p Number of released coefficients (the clip box length).
#' @param rate Optional \code{dp.rate} object from \code{\link{dp.rate.public}}
#'   scaling \code{n} on a two-part positive channel.
#' @return The block count as an integer.  Errors when the cap falls below 2,
#'   since so few rows cannot be subsampled-and-aggregated.
#' @seealso \code{\link{dp.rate.public}},
#'   \code{\link{gen.syn.dp.projected.subagg}},
#'   \code{\link{gen.syn.dp.projected.zi}}.
#' @export
#' @examples
#' dp.subagg.blocks(2000, 6)                                 # 66
#' dp.subagg.blocks(2000, 6, dp.rate.public(0.6, "design"))  # 26
dp.subagg.blocks <- function(n, p, rate = NULL) {
    if (!is.null(rate) && !inherits(rate, "dp.rate"))
        stop("rate must be built by dp.rate.public(), so that its provenance ",
             "is on the record", call. = FALSE)
    s.eff <- if (is.null(rate)) n else n * (1 - rate[["pi0"]])
    s.m <- floor(s.eff / (5 * p))
    if (s.m < 2L)
        stop(sprintf(paste0("dp.subagg.blocks: floor(n_eff / (5 p)) = %d < 2 ",
                            "(n_eff = %.1f, p = %d), too few rows to ",
                            "subsample-and-aggregate"), s.m, s.eff, p),
             call. = FALSE)
    as.integer(s.m)
}
