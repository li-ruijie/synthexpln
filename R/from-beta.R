# Beta-driven synthetic data generation for synthexpln.
# Public function: gen.syn.beta().
# Internal helpers: .gsb.norm.beta, .gsb.norm.family, .gsb.norm.dispersion,
#                   .gsb.build.X, .gsb.sample.y, .gsb.project.
# Companion to projection(): instead of estimating beta-hat from data and
# projecting onto its sufficient-statistics surface, the caller supplies
# beta directly and the function generates a dataset (with X drawn from
# user-specified marginals or supplied verbatim).  No real data needed.

# +------------------------------------------------------------------------+
# | gen.syn.beta() -- public API                                           |
# +------------------------------------------------------------------------+

#' Generate a synthetic dataset from a coefficient vector
#'
#' Produces synthetic data \eqn{(X, y)} whose generating model is the
#' user-supplied \code{beta}.  No real dataset is required: covariates are
#' drawn from caller-specified marginals (or default iid \eqn{N(0,1)}),
#' and the response is sampled from the requested family.
#'
#' The function is maximally polymorphic in its inputs.  \code{beta} is
#' the only required argument; every other slot has a sensible default
#' and accepts multiple input forms (see the parameter list).
#'
#' Two output modes are available:
#' \describe{
#'   \item{\code{method = "simulate"}}{Forward draw
#'     \eqn{y_i \sim \text{Family}(g^{-1}(x_i^\top \beta), \phi)}.
#'     A refit of \code{beta} on the synthetic data wobbles by the usual
#'     sampling variability (\eqn{O(1/\sqrt{n})}).}
#'   \item{\code{method = "project"}}{Forward simulate, then project
#'     \eqn{y} onto the sufficient-statistics surface of the supplied
#'     \code{beta} (the sufficient statistics projection theorem of Li, 2026).
#'     A refit
#'     reproduces \code{beta} at IEEE 754 double-precision machine epsilon
#'     (\eqn{2^{-52} \approx 2.2 \times 10^{-16}}, worst case 13.9 ulp
#'     measured) for Gaussian/identity.  Other GLM families are looser by
#'     orders of magnitude and improve with \eqn{n}, and Poisson is limited
#'     by integer rounding at a level set by the count magnitude.  See the
#'     Recovery accuracy section of \code{\link{projection}} for the
#'     measured grid, and paper Section 6.}
#' }
#'
#' @param beta Numeric coefficient vector.  Named (e.g.
#'   \code{c("(Intercept)" = 1, x1 = -0.5)}) is preserved; unnamed gets
#'   auto-naming as \code{(Intercept), x1, x2, ...}.  The first element
#'   is treated as the intercept unless \code{names(beta)[1]} is some
#'   other name.
#' @param n Sample size.  Default \code{1000L}.  Ignored if \code{X} is
#'   supplied (\code{nrow(X)} is used instead).
#' @param X Optional design data.frame or matrix.  When supplied, the
#'   columns must be a superset of the non-intercept names in \code{beta};
#'   the design matrix is built via \code{stats::model.matrix} from a
#'   formula reconstructed from \code{colnames(X)}.
#' @param marginals One of: \code{NULL} (default, iid \eqn{N(0,1)}); a
#'   named list of one-argument samplers \code{function(n)} keyed by
#'   predictor name; a data.frame with the predictor columns (rows are
#'   bootstrapped).  Ignored if \code{X} is supplied.
#' @param family A \code{family()} object, or a string like
#'   \code{"gaussian"}, \code{"poisson"}, \code{"binomial"},
#'   \code{"gamma"}, \code{"inv.gaussian"}.  A \code{"family/link"} form
#'   such as \code{"gamma/log"} or \code{"binomial/probit"} is also
#'   accepted.  Default \code{stats::gaussian()}.
#' @param dispersion Numeric scalar.  \code{NULL} (default) uses the
#'   family-appropriate convention: 1 for gaussian/gamma/inv.gaussian/
#'   poisson/binomial.  For projection mode, this is the target Pearson
#'   dispersion enforced after rescaling.
#' @param method One of \code{"simulate"} (default, forward draw) or
#'   \code{"project"} (post-hoc projection so refit \eqn{\hat\beta = \beta}).
#' @param seed Optional integer RNG seed.
#'
#' @return An object of class \code{synthexpln}: a list with
#'   \code{$syn} (data.frame, response in column \code{.y}),
#'   \code{$beta} (the supplied vector, named),
#'   \code{$formula} (\code{.y ~ predictors}, usable for refits),
#'   \code{$family} (the resolved family object),
#'   \code{$dispersion}, \code{$method},
#'   \code{$variant = "from-beta"}, \code{$scope = "none"}.
#'
#' @export
#' @examples
#' # Minimum call: just beta, defaults to gaussian/identity, n = 1000.
#' o <- gen.syn.beta(c(`(Intercept)` = 1, x1 = -0.5, x2 = 2), seed = 1)
#' coef(lm(o$formula, o$syn))                    # wobbles by ~SE/sqrt(n)
#'
#' # project mode reproduces beta to double-precision machine epsilon.
#' o <- gen.syn.beta(c(2, 0.8, -0.3), method = "project", seed = 7)
#' coef(lm(o$formula, o$syn))                    # equals c(2, 0.8, -0.3)
#'
#' # Custom marginals.
#' o <- gen.syn.beta(
#'     beta = c(`(Intercept)` = 30, age = 0.4, female = -2.1),
#'     marginals = list(age    = function(n) stats::rnorm(n, 50, 12),
#'                      female = function(n) stats::rbinom(n, 1L, 0.5)),
#'     seed = 5
#' )
#'
#' # Non-Gaussian families via string.
#' o <- gen.syn.beta(c(-1, 0.3, 0.7), family = "poisson/log", seed = 3)
gen.syn.beta <- function(
    beta,
    n          = 1000L,
    X          = NULL,
    marginals  = NULL,
    family     = stats::gaussian(),
    dispersion = NULL,
    method     = c("simulate", "project"),
    seed       = NULL
) {
    method <- match.arg(method)
    if (!is.null(seed)) set.seed(seed)

    beta       <- .gsb.norm.beta(beta)
    family     <- .gsb.norm.family(family)
    dispersion <- .gsb.norm.dispersion(dispersion, family)
    if (!is.null(X)) n <- nrow(X)
    df.X       <- .gsb.build.X(beta, n, X, marginals)

    # build the design matrix and align beta to its columns
    fml.rhs <- stats::reformulate(colnames(df.X), response = NULL)
    m.X     <- stats::model.matrix(fml.rhs, df.X)
    nm.X    <- colnames(m.X)
    miss.b  <- setdiff(nm.X, names(beta))
    if (length(miss.b) > 0L)
        stop("Design-matrix columns not in `beta`: ",
             paste(miss.b, collapse = ", "),
             ".  Check factor levels or supply matching names.")
    common  <- intersect(nm.X, names(beta))
    beta.a  <- `names<-`(rep(0, ncol(m.X)), nm.X)
    beta.a[common] <- beta[common]

    eta <- as.numeric(m.X %*% beta.a)
    mu  <- family[["linkinv"]](eta)
    y   <- .gsb.sample.y(family, mu, dispersion)

    if (method == "project")
        y <- .gsb.project(df.X, y, beta, dispersion, family)

    y.nm <- ".y"
    if (y.nm %in% colnames(df.X))
        stop("Response column `.y` conflicts with a predictor name.")
    d.syn <- cbind(`names<-`(list(y), y.nm), df.X) |> as.data.frame()
    fml   <- stats::reformulate(colnames(df.X), response = y.nm)

    structure(
        list(syn        = d.syn,
             beta       = beta,
             formula    = fml,
             family     = family,
             dispersion = dispersion,
             method     = method,
             variant    = "from-beta",
             scope      = "none"),
        class = "synthexpln"
    )
}

# +------------------------------------------------------------------------+
# | .gsb.norm.beta -- normalise `beta` argument                            |
# +------------------------------------------------------------------------+
# Auto-names unnamed numeric vectors as "(Intercept), x1, x2, ...".
# Named vectors are preserved verbatim.

#' Normalise beta argument (internal)
#' @noRd
.gsb.norm.beta <- function(beta) {
    if (!is.numeric(beta) || length(beta) < 1L)
        stop("`beta` must be a numeric vector of length >= 1.")
    if (is.null(names(beta))) {
        nm <- c("(Intercept)", paste0("x", seq_len(length(beta) - 1L)))
        beta <- `names<-`(beta, nm)
    }
    beta
}

# +------------------------------------------------------------------------+
# | .gsb.norm.family -- normalise `family` argument                        |
# +------------------------------------------------------------------------+
# Accepts (i) a family() object verbatim; (ii) a bare string giving the
# family name with the default link; (iii) a "family/link" string.
# Defers to stats::gaussian() etc. for canonical-link construction.

#' Normalise family argument (internal)
#' @noRd
.gsb.norm.family <- function(family) {
    if (inherits(family, "family"))
        return(family)
    if (!is.character(family) || length(family) != 1L)
        stop("`family` must be a family() object or a string.")
    parts <- strsplit(family, "/", fixed = TRUE)[[1L]]
    fam   <- tolower(parts[[1L]])
    link  <- `if`(length(parts) >= 2L, parts[[2L]], NULL)
    ctor  <- switch(fam,
        gaussian         = stats::gaussian,
        poisson          = stats::poisson,
        binomial         = stats::binomial,
        gamma            = stats::Gamma,
        inv.gaussian     = stats::inverse.gaussian,
        inverse.gaussian = stats::inverse.gaussian,
        stop("Unknown family: ", fam))
    `if`(is.null(link), ctor(), ctor(link = link))
}

# +------------------------------------------------------------------------+
# | .gsb.norm.dispersion -- normalise `dispersion` argument                |
# +------------------------------------------------------------------------+
# Defaults to 1 for every family (matches stats::glm convention, where
# poisson/binomial report dispersion = 1).  Validated as a positive scalar
# when explicitly supplied.

#' Normalise dispersion argument (internal)
#' @noRd
.gsb.norm.dispersion <- function(dispersion, family) {
    if (!is.null(dispersion)) {
        if (!is.numeric(dispersion) || length(dispersion) != 1L ||
            dispersion <= 0)
            stop("`dispersion` must be a positive scalar.")
        return(as.numeric(dispersion))
    }
    fam.nm <- family[["family"]]
    switch(fam.nm,
        gaussian            = 1,
        Gamma               = 1,
        `inverse.gaussian`  = 1,
        poisson             = 1,
        binomial            = 1,
        1)
}

# +------------------------------------------------------------------------+
# | .gsb.build.X -- assemble the predictor data.frame                      |
# +------------------------------------------------------------------------+
# Four precedence levels for sourcing covariates:
#   1. X supplied (data.frame or matrix)   -- used verbatim
#   2. marginals as data.frame             -- bootstrap rows
#   3. marginals as named list of samplers -- call each
#   4. neither                             -- iid N(0,1) per predictor

#' Build covariate data frame (internal)
#' @noRd
.gsb.build.X <- function(beta, n, X, marginals) {
    nm.pred <- setdiff(names(beta), "(Intercept)")

    if (!is.null(X)) {
        if (!is.null(marginals))
            message("`X` supplied; `marginals` ignored.")
        df.X <- as.data.frame(X)
        miss <- setdiff(nm.pred, colnames(df.X))
        if (length(miss) > 0L)
            stop("`X` is missing columns: ", paste(miss, collapse = ", "))
        return(df.X[, nm.pred, drop = FALSE])
    }

    if (is.data.frame(marginals)) {
        miss <- setdiff(nm.pred, colnames(marginals))
        if (length(miss) > 0L)
            stop("`marginals` data.frame is missing columns: ",
                 paste(miss, collapse = ", "))
        idx <- sample(nrow(marginals), n, replace = TRUE)
        return(`rownames<-`(marginals[idx, nm.pred, drop = FALSE], NULL))
    }

    if (is.list(marginals)) {
        miss <- setdiff(nm.pred, names(marginals))
        if (length(miss) > 0L)
            stop("`marginals` list is missing entries: ",
                 paste(miss, collapse = ", "))
        cols <- Map(function(nm) marginals[[nm]](n), nm.pred)
        return(as.data.frame(`names<-`(cols, nm.pred)))
    }

    cols <- Map(function(nm) stats::rnorm(n), nm.pred)
    as.data.frame(`names<-`(cols, nm.pred))
}

# +------------------------------------------------------------------------+
# | .gsb.sample.y -- forward draw from a GLM family                        |
# +------------------------------------------------------------------------+
# Conventions match stats::glm(): gaussian sd = sqrt(phi); Gamma shape
# = 1/phi, scale = mu * phi; poisson/binomial ignore phi.  Inverse
# Gaussian sampling requires the `statmod` package.

#' Sample response from family (internal)
#' @noRd
.gsb.sample.y <- function(family, mu, dispersion) {
    fam.nm <- family[["family"]]
    n      <- length(mu)
    switch(fam.nm,
        gaussian = stats::rnorm(n, mu, sqrt(dispersion)),
        poisson  = stats::rpois(n, pmax(mu, 0)),
        binomial = stats::rbinom(n, 1L, pmin(pmax(mu, 0), 1)),
        Gamma    = stats::rgamma(n, shape = 1 / dispersion,
                                 scale = pmax(mu, 1e-10) * dispersion),
        `inverse.gaussian` = {
            if (!requireNamespace("statmod", quietly = TRUE))
                stop("inverse.gaussian sampling needs the `statmod` package.")
            statmod::rinvgauss(n, mean = pmax(mu, 1e-10),
                               dispersion = dispersion)
        },
        stop("Sampling not implemented for family: ", fam.nm))
}

# +------------------------------------------------------------------------+
# | .gsb.project -- project y onto the score surface at supplied beta      |
# +------------------------------------------------------------------------+
# Mirrors .project.ols / .project.glm in R/projection.R but accepts the
# supplied (X, beta, phi, family) directly: no refit, since beta and phi
# are the targets, not estimates from data.
#
# Constraints enforced:
#   (C1) X' W (y - mu(X beta)) = 0      score equation at the supplied beta
#   (C2) Pearson dispersion = phi       (saturated LL / Pearson scaling)
#
# Three paths, selected by support of the response:
#   gaussian/identity      closed-form OLS-style projection (linear in y)
#   real-support GLM       single weighted score correction + Pearson scale
#   positive-support GLM   Newton on z = log(y) -> y = exp(z) > 0
# Poisson gets probabilistic rounding to integers after the Newton path.

#' Beta-driven projection (internal)
#' @noRd
.gsb.project <- function(df.X, y.prelim, beta, phi, family) {
    fam.nm  <- family[["family"]]
    link.nm <- family[["link"]]

    fml <- stats::reformulate(colnames(df.X), response = NULL)
    m.X <- stats::model.matrix(fml, df.X)
    nm.X    <- colnames(m.X)
    common  <- intersect(nm.X, names(beta))
    miss.b  <- setdiff(nm.X, names(beta))
    if (length(miss.b) > 0L)
        stop("Design columns not in `beta`: ",
             paste(miss.b, collapse = ", "))
    beta.a <- `names<-`(rep(0, ncol(m.X)), nm.X)
    beta.a[common] <- beta[common]

    n <- nrow(m.X)
    p <- ncol(m.X)

    eta   <- as.numeric(m.X %*% beta.a)
    mu    <- family[["linkinv"]](eta)
    dmude <- family[["mu.eta"]](eta)
    V     <- family[["variance"]](pmax(mu, 1e-10))
    w     <- dmude / V

    pos.fam <- fam.nm %in% c("Gamma", "inverse.gaussian", "poisson") ||
               (fam.nm == "quasi" && link.nm %in% c("log", "inverse"))
    int.fam <- fam.nm == "poisson"

    # ---- PATH 1: gaussian / identity (closed form) -----------------------
    if (fam.nm == "gaussian" && link.nm == "identity") {
        l.qr  <- qr(m.X)
        rk    <- l.qr[["rank"]]
        e.cop <- y.prelim - qr.fitted(l.qr, y.prelim)
        rss.c <- sum(e.cop^2)
        df.r  <- n - rk
        scal  <- `if`(rss.c > 0 && df.r > 0,
                      sqrt(phi * df.r / rss.c), 1)
        return(mu + e.cop * scal)
    }

    # ---- PATH 2: real-support GLM (e.g. gaussian/log) --------------------
    if (!pos.fam) {
        wX     <- sweep(m.X, 1, w, `*`)
        XtwX   <- crossprod(wX)
        diag(XtwX) <- diag(XtwX) + 1e-10
        score  <- as.numeric(crossprod(wX, y.prelim - mu))
        corr   <- as.numeric(wX %*% solve(XtwX, score))
        y.pr   <- y.prelim - corr
        e.pr   <- y.pr - mu
        sd.s   <- sqrt(pmax(V, 1e-10))
        pe     <- e.pr / sd.s
        phi.pr <- sum(pe^2) / max(n - p, 1L)
        if (phi.pr > 1e-10) {
            sc   <- sqrt(phi / phi.pr)
            y.pr <- mu + sd.s * pe * sc
        }
        return(y.pr)
    }

    # ---- PATH 3: positive-support GLM, Newton on log-scale ---------------
    z       <- log(pmax(y.prelim, 1e-10))
    maxit   <- 30L
    tol     <- 1e-8
    it      <- 0L
    repeat {
        it <- it + 1L
        y  <- exp(z)
        sc <- as.numeric(crossprod(m.X, w * (y - mu)))
        if (max(abs(sc)) < tol || it > maxit) break
        wy   <- w * y
        WYX  <- sweep(m.X, 1, wy, `*`)
        JJt  <- crossprod(WYX)
        diag(JJt) <- diag(JJt) * (1 + 1e-10)
        lam  <- tryCatch(solve(JJt, sc), error = function(e) rep(0, p))
        dz   <- as.numeric(WYX %*% lam)
        step <- min(1, 3 / (max(abs(dz)) + 1e-16))
        z    <- z - step * dz
    }
    y.nw <- exp(z)

    sd.s   <- sqrt(pmax(V, 1e-10))
    pe     <- (y.nw - mu) / sd.s
    phi.nw <- sum(pe^2) / max(n - p, 1L)
    y.pr   <- `if`(phi.nw > 1e-10,
                   pmax(mu + sd.s * pe * sqrt(phi / phi.nw), 1e-10),
                   y.nw)

    # Poisson -> probabilistic rounding to non-negative integers.
    # The integer constraint is the source of the ~2.3% residual error
    # documented in Section 6 of the paper (also affects projection()).
    if (int.fam) {
        fr   <- y.pr - floor(y.pr)
        y.pr <- pmax(ifelse(stats::runif(n) < fr,
                            ceiling(y.pr), floor(y.pr)), 0)
    }
    y.pr
}
