# privacy.R -- Privacy metrics for synthexpln synthetic data.
# Public functions: privacy.metrics().

# +------------------------------------------------------------------------+
# | .is.numeric.col -- internal column-type classifier                    |
# +------------------------------------------------------------------------+

# Returns TRUE if the column should be treated as continuous (numeric with
# more than 10 unique values, not factor/character/logical).
.is.numeric.col <- function(v) {
    is.numeric(v) && !is.factor(v) && length(unique(v)) > 10L
}

# +------------------------------------------------------------------------+
# | privacy.metrics -- DCR, hit rate, membership AUC                      |
# +------------------------------------------------------------------------+

#' Privacy metrics for synthetic data
#'
#' Computes one or more privacy-risk metrics comparing a synthetic dataset
#' to the original data.  All metrics are interpretable as "lower is more
#' private", except membership AUC where values near 0.5 indicate strong
#' privacy (the synthetic data is statistically indistinguishable from the
#' original).
#'
#' @section Distance to Closest Record (DCR):
#' Mean Gower distance from each synthetic row to its nearest-neighbour
#' original row.  Gower distance normalises numeric columns by the observed
#' range of the original data and treats factor/character mismatches as 1.
#' For speed, at most 500 synthetic rows are sampled.  A high DCR indicates
#' that synthetic records are well separated from originals (low memorisation
#' risk).
#'
#' @section Hit Rate:
#' Fraction of synthetic rows that are exact duplicates of some original
#' row.  Computed by converting each row to a string key and checking
#' membership.  For \code{syn = orig} (perfect memorisation) the hit rate
#' equals 1.
#'
#' @section Membership AUC:
#' Cross-validated AUC of a logistic classifier trained to distinguish
#' original records (label 1) from synthetic records (label 0).  Values
#' near 0.5 indicate the classifier cannot distinguish the datasets
#' (desirable).  Computed with 5-fold CV; returns 0.5 on model failure.
#'
#' @param syn A data frame of synthetic data, or a \code{synthexpln} object
#'   (in which case \code{syn$syn} is used automatically).
#' @param orig A data frame of original (reference) data.
#' @param methods Character vector; one or more of \code{"dcr"},
#'   \code{"hit"}, \code{"membership"}.  Only the requested metrics are
#'   computed.
#' @return A named list containing only the requested metrics:
#'   \describe{
#'     \item{\code{dcr}}{Numeric scalar; mean Gower distance to nearest
#'       original record.  Only present when \code{"dcr"} is requested.}
#'     \item{\code{hit.rate}}{Numeric scalar in \eqn{[0,1]}; fraction of
#'       synthetic rows exactly matching some original row.  Only present
#'       when \code{"hit"} is requested.}
#'     \item{\code{membership.auc}}{Numeric scalar in \eqn{[0,1]}; AUC of
#'       membership-inference classifier.  Only present when
#'       \code{"membership"} is requested.}
#'   }
#' @export
#' @examples
#' set.seed(1)
#' orig <- data.frame(a = rnorm(100), b = rnorm(100))
#' syn  <- data.frame(a = rnorm(100), b = rnorm(100))
#' privacy.metrics(syn, orig, methods = c("dcr", "hit"))
privacy.metrics <- function(syn, orig,
                            methods = c("dcr", "hit", "membership")) {
    methods <- match.arg(methods, several.ok = TRUE)

    # Unwrap synthexpln objects transparently.
    if (inherits(syn, "synthexpln"))
        syn <- syn[["syn"]]

    # Align to common columns.
    v.cols <- intersect(names(syn), names(orig))
    if (length(v.cols) == 0L)
        stop("'syn' and 'orig' share no column names")
    syn  <- syn[, v.cols, drop = FALSE]
    orig <- orig[, v.cols, drop = FALSE]

    out <- list()

    # ------------------------------------------------------------------
    # DCR -- mean Gower distance to nearest original record
    # ------------------------------------------------------------------
    if ("dcr" %in% methods) {
        v.cnt <- v.cols[mapply(.is.numeric.col, orig)]
        v.cat <- setdiff(v.cols, v.cnt)

        # Range for normalisation; guard against zero-range.
        l.range <- `names<-`(
            Map(function(nm) {
                s.r <- diff(range(orig[[nm]], na.rm = TRUE))
                if (s.r > 0) s.r else 1
            }, v.cnt),
            v.cnt)

        s.p <- length(v.cnt) + length(v.cat)

        # Sample for speed (O(n_syn * n_orig) otherwise).
        s.ns  <- min(500L, nrow(syn))
        v.idx <- sample(nrow(syn), s.ns)

        v.dcr <- mapply(function(i) {
            # For each synthetic row i, compute min Gower distance to any original row.
            # Accumulate distance across columns via Reduce, then take min over rows.
            v.row.dists <- mapply(function(j) {
                s.d.cnt <- Reduce(function(acc, nm)
                    acc + abs(syn[[nm]][i] - orig[[nm]][j]) / l.range[[nm]],
                    v.cnt, init = 0)
                s.d.cat <- Reduce(function(acc, nm)
                    acc + as.integer(
                        as.character(syn[[nm]][i]) !=
                        as.character(orig[[nm]][j])),
                    v.cat, init = 0)
                (s.d.cnt + s.d.cat) / s.p
            }, seq_len(nrow(orig)))
            min(v.row.dists)
        }, v.idx)

        out[["dcr"]] <- mean(v.dcr)
    }

    # ------------------------------------------------------------------
    # Hit rate -- fraction of synthetic rows exactly in original
    # ------------------------------------------------------------------
    if ("hit" %in% methods) {
        # Convert each row to a string key for O(n log n) matching.
        key.orig <- apply(orig, 1L, paste, collapse = "|")
        key.syn  <- apply(syn,  1L, paste, collapse = "|")
        out[["hit.rate"]] <- mean(key.syn %in% key.orig)
    }

    # ------------------------------------------------------------------
    # Membership AUC -- 5-fold cross-validated logistic classifier
    # ------------------------------------------------------------------
    if ("membership" %in% methods) {
        s.n  <- min(nrow(syn), nrow(orig))
        d.both <- rbind(
            cbind(orig[sample(nrow(orig), s.n), , drop = FALSE], .real = 1L),
            cbind(syn[ sample(nrow(syn),  s.n), , drop = FALSE], .real = 0L))

        v.fold <- sample(rep(1:5, length.out = nrow(d.both)))

        # For each fold k, return a named list: idx -> predictions.
        # Then combine into v.pred by position.
        v.pred <- Reduce(function(acc, k) {
            b.test  <- v.fold == k
            b.train <- !b.test
            mo <- tryCatch(
                stats::glm(.real ~ ., data = d.both[b.train, ],
                           family = stats::binomial(),
                           control = stats::glm.control(maxit = 50L)),
                error = function(e) NULL)
            fold.pred <- if (!is.null(mo)) {
                tryCatch(
                    stats::predict(mo, d.both[b.test, ], type = "response"),
                    error = function(e) rep(0.5, sum(b.test)))
            } else {
                rep(0.5, sum(b.test))
            }
            acc[b.test] <- fold.pred
            acc
        }, 1:5, init = numeric(nrow(d.both)))

        v.real <- d.both[[".real"]]
        v.pos  <- v.pred[v.real == 1L]
        v.neg  <- v.pred[v.real == 0L]
        # Cap sample size for AUC computation to avoid O(n^2) blow-up.
        v.ps   <- sample(v.pos, min(2000L, length(v.pos)))
        v.ns   <- sample(v.neg, min(2000L, length(v.neg)))
        out[["membership.auc"]] <-
            mean(outer(v.ps, v.ns, ">") + 0.5 * outer(v.ps, v.ns, "=="))
    }

    out
}
