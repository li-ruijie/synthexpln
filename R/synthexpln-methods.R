#' S3 methods for synthexpln objects
#'
#' Print, summarise, and plot synthetic projection objects.
#'
#' @param x,object A \code{synthexpln} object.
#' @param ... Additional arguments (unused).
#' @return
#'   \code{summary.synthexpln} returns an object of class
#'   \code{summary.synthexpln}, a list holding the release variant, scope,
#'   privacy parameters \code{epsilon} and \code{delta}, the row count, a
#'   coefficient table, and the maximum absolute percentage coefficient
#'   deviation. The \code{print} methods return their argument invisibly. The
#'   \code{plot} method returns the \code{synthexpln} object invisibly, and all
#'   three are called for their side effect of printing a summary or drawing a
#'   diagnostic plot of synthetic against original coefficients.
#' @name synthexpln-methods
NULL

#' @rdname synthexpln-methods
#' @exportS3Method base::print
print.synthexpln <- function(x, ...) {
    cat("<synthexpln>\n")
    cat("  variant:", x[["variant"]], "\n")
    cat("  scope:  ", x[["scope"]], "\n")
    if (!is.null(x[["epsilon"]])) {
        cat(sprintf("  epsilon: %g, delta: %g\n", x[["epsilon"]], x[["delta"]]))
    }
    cat("  rows:   ", nrow(x[["syn"]]), "\n")
    invisible(x)
}

#' @rdname synthexpln-methods
#' @exportS3Method base::summary
summary.synthexpln <- function(object, ...) {
    # build a deviations table from beta.dp (DP projection) or beta (non-DP)
    beta.ref <- if (!is.null(object[["beta.dp"]])) object[["beta.dp"]] else object[["beta"]]
    coef.table <- data.frame(
        term = names(beta.ref),
        beta = unname(beta.ref),
        row.names = NULL
    )
    if (!is.null(object[["beta.hat"]]) && !is.null(object[["beta.dp"]])) {
        coef.table[["beta.hat"]] <- unname(object[["beta.hat"]])
        coef.table[["pct.dev"]]  <- 100 * (coef.table[["beta"]] - coef.table[["beta.hat"]]) /
                               pmax(abs(coef.table[["beta.hat"]]), .Machine[["double.eps"]])
        max.abs.pct <- max(abs(coef.table[["pct.dev"]]), na.rm = TRUE)
    } else {
        max.abs.pct <- 0
    }
    structure(list(
        variant     = object[["variant"]],
        scope       = object[["scope"]],
        epsilon     = object[["epsilon"]],
        delta       = object[["delta"]],
        rows        = nrow(object[["syn"]]),
        coef        = coef.table,
        max.abs.pct = max.abs.pct
    ), class = "summary.synthexpln")
}

#' @rdname synthexpln-methods
#' @exportS3Method base::print
print.summary.synthexpln <- function(x, ...) {
    cat("<summary.synthexpln>\n")
    cat("  variant:", x[["variant"]], "\n")
    cat("  scope:  ", x[["scope"]], "\n")
    if (!is.null(x[["epsilon"]]))
        cat(sprintf("  epsilon: %g, delta: %g\n", x[["epsilon"]], x[["delta"]]))
    cat("  rows:   ", x[["rows"]], "\n")
    cat("  max|%|: ", sprintf("%.4f", x[["max.abs.pct"]]), "\n")
    cat("\n  Coefficients:\n")
    print(x[["coef"]], row.names = FALSE)
    invisible(x)
}

#' @rdname synthexpln-methods
#' @importFrom graphics abline
#' @exportS3Method base::plot
plot.synthexpln <- function(x, ...) {
    # minimal base-R plot: syn coef (or beta.dp) vs beta.hat if available
    beta.ref <- if (!is.null(x[["beta.dp"]])) x[["beta.dp"]] else x[["beta"]]
    if (is.null(x[["beta.hat"]])) {
        message("No beta.hat slot; nothing to plot.")
        return(invisible(x))
    }
    # Let caller override main/xlab/ylab via `...`; supply defaults otherwise.
    dots <- list(...)
    defaults <- list(
        xlab = "beta.hat (original)",
        ylab = "beta.dp (synthetic)",
        main = paste0("synthexpln ", x[["variant"]])
    )
    args <- c(list(x = x[["beta.hat"]], y = beta.ref),
              dots,
              defaults[setdiff(names(defaults), names(dots))])
    do.call(plot, args)
    abline(0, 1, lty = 2, col = "red")
    invisible(x)
}
