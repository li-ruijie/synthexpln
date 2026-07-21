# Stratified Gaussian-copula base generator for synthexpln.
# Internal functions: .id.cat, .detect.support, .make.transform,
#                     .cb.std.cumulants, .cb.cf.quantile, .cf.is.monotone,
#                     .cf.safe.order, .shrink.moments, .knn.strata,
#                     .copula.base.
# Ported from lib/gen-syn-family.r (read-only source), pseudo-random
# method only (the Sobol branch and its 'qrng' dependency are not
# ported). .copula.base replaces .bootstrap.rows as the preliminary
# synthetic data when projection(base = "copula") is requested, and
# strata.exact = TRUE draws the stratum labels as a permutation of the
# original labels so every cross-classification cell keeps its exact
# original count (the design-Gram projection's condition (i)).

# Classify a variable as categorical ("cat"), binary ("bin"), or
# continuous ("cnt").
.id.cat <- function(v, threshold = 10) {
    b.cat <- "cnt"
    v.ulen <- length(unique(v))
    if (is.factor(v) | is.character(v) | is.logical(v))
        b.cat <- "cat"
    if (v.ulen <= threshold)
        b.cat <- "cat"
    if ((b.cat == "cat") && (v.ulen == 2))
        b.cat <- "bin"
    b.cat
}

# Infer the support of a continuous variable.
.detect.support <- function(v) {
    if (all(v > 0))                                         return("positive")
    if (any(v == 0) && all(v >= 0) && mean(v == 0) > 0.01)  return("zeroinfl")
    "real"
}

# Support transforms h and h^{-1}.
.make.transform <- function(support, v = NULL) {
    switch(support,
        real = list(
            h    = identity,
            hinv = identity),
        positive = list(
            h    = log,
            hinv = exp),
        zeroinfl = list(
            h    = log,
            hinv = exp),
        count = {
            s.shift <- 0.5
            list(h    = function(x) log(x + s.shift),
                 hinv = function(x) pmax(round(exp(x) - s.shift), 0))
        },
        bounded = {
            s.lo  <- min(v, na.rm = TRUE)
            s.hi  <- max(v, na.rm = TRUE)
            s.eps <- 1 / (2 * length(v))
            list(h = function(x) {
                     x.sc <- (x - s.lo) / (s.hi - s.lo)
                     x.sc <- pmin(pmax(x.sc, s.eps), 1 - s.eps)
                     stats::qlogis(x.sc)
                 },
                 hinv = function(x) {
                     s.lo + (s.hi - s.lo) * stats::plogis(x)
                 })
        },
        stop(paste("Unknown support:", support)))
}

# Standardised cumulants gamma_3, ..., gamma_k of a sample.
.cb.std.cumulants <- function(v, order) {
    s.n <- length(v)
    s.sigma <- stats::sd(v)
    if (s.sigma < .Machine[["double.eps"]] || s.n < order + 1)
        return(rep(0, max(order - 2, 0)))
    v.c <- (v - mean(v)) / s.sigma
    v.mu <- mapply(function(r) mean(v.c^r), seq(3, order))
    v.gamma <- numeric(length(v.mu))
    if (order >= 3) v.gamma[1] <- v.mu[1]
    if (order >= 4) v.gamma[2] <- v.mu[2] - 3
    if (order >= 5) v.gamma[3] <- v.mu[3] - 10 * v.mu[1]
    if (order >= 6) v.gamma[4] <- v.mu[4] - 15 * v.mu[2] -
                                   10 * v.mu[1]^2 + 30
    v.gamma
}

# Cornish-Fisher quantile function.
.cb.cf.quantile <- function(p, mu, sigma, gamma, order) {
    z <- stats::qnorm(p)
    w <- z
    if (order >= 3) {
        g1 <- gamma[1]
        w <- w + (g1 / 6) * (z^2 - 1)
    }
    if (order >= 4) {
        g2 <- gamma[2]
        w <- w + (g2 / 24) * (z^3 - 3 * z) -
                 (g1^2 / 36) * (2 * z^3 - 5 * z)
    }
    if (order >= 5) {
        g3 <- gamma[3]
        w <- w + (g3 / 120) * (z^4 - 6 * z^2 + 3) -
                 (g1 * g2 / 24) * (z^4 - 5 * z^2 + 2) +
                 (g1^3 / 324) * (12 * z^4 - 53 * z^2 + 17)
    }
    if (order >= 6) {
        g4 <- gamma[4]
        w <- w + (g4 / 720) * (z^5 - 10 * z^3 + 15 * z) -
                 (g1 * g3 / 180) * (2 * z^5 - 17 * z^3 + 21 * z) -
                 (g2^2 / 192) * (z^5 - 7 * z^3 + 9 * z) +
                 (g1^2 * g2 / 216) * (4 * z^5 - 32 * z^3 + 41 * z) -
                 (g1^4 / 3888) * (14 * z^5 - 118 * z^3 + 155 * z)
    }
    mu + sigma * w
}

# Monotonicity check for the CF expansion over (0, 1).
.cf.is.monotone <- function(mu, sigma, gamma, order) {
    v.p <- seq(0.001, 0.999, length.out = 500)
    v.q <- .cb.cf.quantile(v.p, mu, sigma, gamma, order)
    all(diff(v.q) > 0)
}

# Highest order <= `order` at which the CF expansion stays monotone.
.cf.safe.order <- function(mu, sigma, gamma, order) {
    if (sigma < .Machine[["double.eps"]]) return(2)
    while (order > 2) {
        b.mono <- tryCatch(.cf.is.monotone(mu, sigma, gamma, order),
                           error = function(e) FALSE)
        if (isTRUE(b.mono)) break
        order <- order - 1
    }
    order
}

# Empirical Bayes shrinkage of per-stratum means toward the grand mean.
.shrink.moments <- function(l.strata.cf, l.strata.rows, v.cnt.nm) {
    Map(function(nm) {
        v.mu <- mapply(function(l) l[[nm]][["mu"]], l.strata.cf)
        v.sigma <- mapply(function(l) l[[nm]][["sigma"]], l.strata.cf)
        v.ns <- mapply(length, l.strata.rows)
        b.valid <- !is.na(v.mu) & v.sigma > 0

        if (sum(b.valid) < 3) return(l.strata.cf)

        s.grand <- mean(v.mu[b.valid])
        s.sigma2.pool <- mean(v.sigma[b.valid]^2)
        s.tau2 <- max(0, stats::var(v.mu[b.valid]) -
                       s.sigma2.pool * mean(1 / v.ns[b.valid]))
        v.w <- mapply(function(ns) {
            `if`(s.tau2 + s.sigma2.pool / ns > 0,
                 s.tau2 / (s.tau2 + s.sigma2.pool / ns),
                 1)
        }, v.ns)

        Map(function(s.key, w) {
            s.mu.old <- l.strata.cf[[s.key]][[nm]][["mu"]]
            if (!is.na(s.mu.old))
                l.strata.cf[[s.key]][[nm]][["mu"]] <<-
                    w * s.mu.old + (1 - w) * s.grand
        }, names(l.strata.cf), v.w)
    }, v.cnt.nm)
    l.strata.cf
}

# K-means strata when no categorical variables exist.
.knn.strata <- function(d.cnt, k = 30) {
    s.n <- nrow(d.cnt)
    s.k <- min(k, s.n - 1)
    m.sc <- scale(d.cnt)
    m.sc[is.na(m.sc)] <- 0
    s.n.clusters <- max(1, round(s.n / s.k))
    l.km <- stats::kmeans(m.sc, centers = s.n.clusters,
                          nstart = 5, iter.max = 50)
    factor(l.km[["cluster"]])
}

# The stratified Gaussian-copula base generator. Faithful port of
# gen.syn.family() minus the Sobol branch; the caller manages the RNG
# seed. strata.exact = TRUE draws the stratum labels without
# replacement (a permutation, exact cell counts), and requires the
# synthetic size to equal nrow(d).
.copula.base <- function(d, order = 4, supports = NULL, shrink = TRUE,
                         strata.exact = FALSE) {
    s.n <- nrow(d)
    v.types <- mapply(.id.cat, d)
    v.cat.nm <- names(v.types[v.types != "cnt"])
    v.cnt.nm <- names(v.types[v.types == "cnt"])

    # Stage 1: supports and transforms
    if (is.null(supports)) {
        supports <- `names<-`(mapply(.detect.support, d[v.cnt.nm]), v.cnt.nm)
    }
    v.missing <- v.cnt.nm[!v.cnt.nm %in% names(supports)]
    if (length(v.missing) > 0)
        supports <- c(supports,
                      `names<-`(mapply(.detect.support, d[v.missing]),
                                v.missing))
    l.transforms <- Map(function(nm) .make.transform(supports[[nm]], d[[nm]]),
                        v.cnt.nm)

    # Stage 1b: zero-inflated indicator joins the stratification
    v.zi.nm <- names(supports[supports == "zeroinfl"])
    d.aug <- d
    if (length(v.zi.nm) > 0) {
        Map(function(nm) {
            s.ind.nm <- paste0(nm, ".pos")
            d.aug[[s.ind.nm]] <<- factor(as.integer(d[[nm]] > 0))
        }, v.zi.nm)
        v.cat.nm <- c(v.cat.nm, paste0(v.zi.nm, ".pos"))
    }

    # Stage 2: stratification and the label draw
    if (length(v.cat.nm) > 0) {
        v.strata <- interaction(d.aug[v.cat.nm], drop = TRUE)
    } else {
        v.strata <- .knn.strata(d[v.cnt.nm])
    }
    l.strata.rows <- split(seq_len(s.n), v.strata)

    if (strata.exact) {
        if (s.n != nrow(d))
            stop("strata.exact = TRUE requires synthetic size == nrow(d)")
        v.strata.sampled <- sample(v.strata, s.n, replace = FALSE)
    } else {
        v.strata.sampled <- sample(v.strata, s.n, replace = TRUE)
    }

    if (length(v.cat.nm) > 0) {
        l.strata.cat <- Map(function(rows)
                            d.aug[rows[1], v.cat.nm, drop = FALSE],
                            l.strata.rows)
        d.syn.cat <- `rownames<-`(Reduce(rbind,
                                          Map(function(k)
                                              l.strata.cat[[as.character(k)]],
                                              v.strata.sampled)), NULL)
    } else {
        d.syn.cat <- NULL
    }

    # Stage 3: per-stratum CF parameters on the h-scale
    l.strata.cf <- Map(function(s.key, v.rows) {
        Map(function(nm) {
            v.raw <- d[[nm]][v.rows]
            l.tr <- l.transforms[[nm]]
            if (supports[[nm]] == "zeroinfl") {
                v.raw <- v.raw[v.raw > 0]
                if (length(v.raw) == 0)
                    return(list(mu = NA, sigma = 0,
                                gamma = rep(0, max(order - 2, 0)),
                                order = 2, skip = TRUE))
            }
            v.h <- l.tr[["h"]](v.raw)
            v.h <- v.h[is.finite(v.h)]
            if (length(v.h) < 2)
                return(list(mu = `if`(length(v.h) == 1, v.h, 0),
                            sigma = 0,
                            gamma = rep(0, max(order - 2, 0)),
                            order = 2, skip = FALSE))
            s.mu <- mean(v.h)
            s.sigma <- stats::sd(v.h)
            v.gamma <- .cb.std.cumulants(v.h, order)
            s.order <- .cf.safe.order(s.mu, s.sigma, v.gamma, order)
            list(mu = s.mu, sigma = s.sigma,
                 gamma = v.gamma, order = s.order, skip = FALSE)
        }, v.cnt.nm)
    }, names(l.strata.rows), l.strata.rows)

    # Stage 4: empirical Bayes shrinkage
    if (shrink)
        l.strata.cf <- .shrink.moments(l.strata.cf, l.strata.rows, v.cnt.nm)

    # Stage 5: pooled within-stratum correlation on the h-scale
    m.h <- mapply(function(nm) {
        l.tr <- l.transforms[[nm]]
        v.raw <- d[[nm]]
        if (supports[[nm]] == "zeroinfl")
            v.raw <- ifelse(v.raw > 0, v.raw, NA)
        v.h <- l.tr[["h"]](v.raw)
        v.h[!is.finite(v.h)] <- NA
        v.h
    }, v.cnt.nm)
    m.resid <- mapply(function(nm, j) {
        v <- m.h[, j]
        v.fitted <- mapply(function(rows) mean(v[rows], na.rm = TRUE),
                           l.strata.rows)
        v - v.fitted[as.integer(v.strata)]
    }, v.cnt.nm, seq_along(v.cnt.nm))
    s.p.cnt <- length(v.cnt.nm)
    m.cor.within <- `if`(s.p.cnt > 1,
                         stats::cor(m.resid, use = "pairwise.complete.obs"),
                         matrix(1, 1, 1))
    m.cor.within[is.na(m.cor.within)] <- 0
    diag(m.cor.within) <- 1
    v.eigen <- eigen(m.cor.within, only.values = TRUE)[["values"]]
    if (any(v.eigen < 0))
        m.cor.within <- as.matrix(Matrix::nearPD(m.cor.within,
                                                  corr = TRUE)[["mat"]])

    # Stage 6: correlated uniforms via the Gaussian copula (pseudo)
    m.z <- matrix(stats::rnorm(s.n * s.p.cnt), nrow = s.n, ncol = s.p.cnt)
    if (s.p.cnt > 1) {
        m.z <- m.z %*% solve(chol(stats::cov(m.z)))
        m.z <- m.z %*% chol(m.cor.within)
    }
    m.u <- stats::pnorm(m.z)

    # Stage 7: CF quantile then h^{-1}
    l.syn.idx <- split(seq_len(s.n), v.strata.sampled)
    m.syn.cnt <- matrix(NA_real_, nrow = s.n, ncol = s.p.cnt)
    m.syn.cnt <- `colnames<-`(m.syn.cnt, v.cnt.nm)
    Map(function(s.key, v.syn.idx) {
        l.cf <- l.strata.cf[[s.key]]
        Map(function(nm, j) {
            l.par <- l.cf[[nm]]
            l.tr  <- l.transforms[[nm]]
            if (supports[[nm]] == "zeroinfl") {
                s.ind.nm <- paste0(nm, ".pos")
                b.pos <- d.syn.cat[[s.ind.nm]][v.syn.idx] == "1"
                m.syn.cnt[v.syn.idx[!b.pos], j] <<- 0
                if (!any(b.pos)) return(NULL)
                v.syn.pos <- v.syn.idx[b.pos]
            } else {
                v.syn.pos <- v.syn.idx
            }
            if (isTRUE(l.par[["skip"]]) || l.par[["sigma"]] == 0) {
                m.syn.cnt[v.syn.pos, j] <<-
                    l.tr[["hinv"]](`if`(is.na(l.par[["mu"]]), 0,
                                        l.par[["mu"]]))
                return(NULL)
            }
            v.h <- .cb.cf.quantile(m.u[v.syn.pos, j],
                                l.par[["mu"]], l.par[["sigma"]],
                                l.par[["gamma"]], l.par[["order"]])
            m.syn.cnt[v.syn.pos, j] <<- l.tr[["hinv"]](v.h)
        }, v.cnt.nm, seq_along(v.cnt.nm))
    }, names(l.syn.idx), l.syn.idx)

    # Stage 8: assemble, dropping the internal *.pos indicator columns
    v.orig.cat <- names(v.types[v.types != "cnt"])
    d.syn <- `if`(!is.null(d.syn.cat),
                  cbind(d.syn.cat[v.orig.cat],
                        as.data.frame(m.syn.cnt)),
                  as.data.frame(m.syn.cnt))
    `rownames<-`(d.syn[colnames(d)], NULL)
}
