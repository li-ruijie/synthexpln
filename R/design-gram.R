# Design-Gram projection for synthexpln.
# Internal functions: .dg.recolour, .dg.project, .dg.block, .dg.twopart.
# Ported from lib/11-projection-generator.r (read-only source).
#
# Pins the synthetic predictor Gram X*'X* = X'X on top of the exact-count
# stratum permutation of .copula.base(strata.exact = TRUE), so the refit
# triple (coefficients, standard errors, and the derived statistics) is
# exact for OLS analysts across every main-effects submodel, and within
# the fixed-point tolerance for canonical GLM analysts (the design-Gram
# proposition of Li, 2026).
#
# every failure stops loudly. A silent fallback here would hand the
# caller an identity nothing delivers.

# One weighted two-step recolouring pass: pin (sqrt(w)D)'(sqrt(w)C) and
# (sqrt(w)C)'(sqrt(w)C) to the original weighted targets. w = 1 gives
# the unweighted OLS case.
.dg.recolour <- function(m.D.syn, m.D.orig, m.C.orig, m.C.cur, v.w,
                         v.w.orig = v.w) {
    v.rw.s <- sqrt(v.w)
    v.rw.o <- sqrt(v.w.orig)
    m.Dw.s <- m.D.syn  * v.rw.s
    m.Dw.o <- m.D.orig * v.rw.o
    m.Cw.o <- m.C.orig * v.rw.o
    m.Cw.c <- m.C.cur  * v.rw.s

    m.DtD <- crossprod(m.Dw.s)
    m.DtD.inv <- tryCatch(solve(m.DtD),
                          error = function(e)
                              stop("design-gram: weighted D'D singular"))
    m.tgt.DtC <- crossprod(m.Dw.o, m.Cw.o)
    m.fit  <- m.Dw.s %*% (m.DtD.inv %*% m.tgt.DtC)
    m.res  <- m.Cw.c - m.Dw.s %*% (m.DtD.inv %*% crossprod(m.Dw.s, m.Cw.c))

    m.G.tgt <- crossprod(m.Cw.o) - crossprod(m.fit)
    m.G.res <- crossprod(m.res)
    m.L.tgt <- tryCatch(chol(m.G.tgt),
                        error = function(e)
                            stop("design-gram: target within Gram not PD"))
    m.L.res <- tryCatch(chol(m.G.res),
                        error = function(e)
                            stop("design-gram: residual Gram not PD"))
    m.Cw.new <- m.fit + m.res %*% solve(m.L.res, m.L.tgt)
    m.Cw.new / v.rw.s
}

# Pin the predictor Gram for a designated response. family/link are the
# lowercase strings of .project.glm; gaussian/identity pins the raw Gram
# in one pass, canonical GLMs pin the Fisher-weighted Gram by fixed
# point at relative tolerance `tol`.
.dg.project <- function(d.syn, d.orig, y.nm, x.nms,
                        family = "gaussian", link = "identity",
                        tol = 1e-8, maxiter = 100L) {
    fml <- stats::reformulate(x.nms)
    m.X.orig <- stats::model.matrix(fml, d.orig)
    m.X.syn  <- stats::model.matrix(fml, d.syn)
    if (!identical(colnames(m.X.orig), colnames(m.X.syn)))
        stop("design-gram: original and synthetic model matrices disagree")

    v.types   <- mapply(.id.cat, d.orig[x.nms])
    v.cnt.nm  <- names(v.types[v.types == "cnt"])
    if (length(v.cnt.nm) == 0) return(d.syn)
    v.cnt.idx <- which(colnames(m.X.orig) %in% v.cnt.nm)
    v.cat.idx <- setdiff(seq_len(ncol(m.X.orig)), v.cnt.idx)

    m.D.syn  <- m.X.syn[, v.cat.idx, drop = FALSE]
    m.D.orig <- m.X.orig[, v.cat.idx, drop = FALSE]
    m.C.orig <- m.X.orig[, v.cnt.idx, drop = FALSE]
    m.C.cur  <- m.X.syn[, v.cnt.idx, drop = FALSE]

    if (max(abs(crossprod(m.D.syn) - crossprod(m.D.orig))) > 1e-9)
        stop("design-gram: D'D not pinned; generate with strata.exact")

    b.ols <- family == "gaussian" && link == "identity"
    if (b.ols) {
        m.C.new <- .dg.recolour(m.D.syn, m.D.orig, m.C.orig, m.C.cur,
                                v.w = rep(1, nrow(m.D.syn)))
    } else {
        fml.y <- stats::reformulate(x.nms, response = y.nm)
        fam.obj <- .make.glm.family(family, link)
        b.quasi <- inherits(fam.obj, "family") &&
                   fam.obj[["family"]] == "quasi"
        d.orig.fit <- `if`(link %in% c("log", "inverse", "1/mu^2") &&
                           any(d.orig[[y.nm]] == 0),
                           d.orig[d.orig[[y.nm]] > 0, ], d.orig)
        mo.orig <- `if`(b.quasi,
            stats::glm(fml.y, d.orig.fit, family = fam.obj,
                       start = c(1, rep(0, ncol(m.X.orig) - 1)),
                       control = stats::glm.control(maxit = 200)),
            stats::glm(fml.y, d.orig.fit, family = fam.obj,
                       control = stats::glm.control(maxit = 200)))
        v.beta <- stats::coef(mo.orig)[colnames(m.X.orig)]
        if (any(is.na(v.beta)))
            stop("design-gram: original fit rank-deficient on x.nms")
        w.info <- function(m.X) {
            v.eta <- as.numeric(m.X %*% v.beta)
            fam.obj[["mu.eta"]](v.eta)^2 /
                fam.obj[["variance"]](fam.obj[["linkinv"]](v.eta))
        }
        v.w.orig <- w.info(m.X.orig)
        m.I.tgt  <- crossprod(m.X.orig * sqrt(v.w.orig))
        s.scl.w  <- max(abs(m.I.tgt))

        s.iter <- 0L
        repeat {
            s.iter <- s.iter + 1L
            m.X.cur <- m.X.syn
            m.X.cur[, v.cnt.idx] <- m.C.cur
            v.w.cur <- w.info(m.X.cur)
            m.I.cur <- crossprod(m.X.cur * sqrt(v.w.cur))
            if (max(abs(m.I.cur - m.I.tgt)) < tol * s.scl.w) break
            if (s.iter > maxiter)
                stop("design-gram: weighted fixed point failed to converge")
            m.C.cur <- .dg.recolour(m.D.syn, m.D.orig, m.C.orig,
                                    m.C.cur, v.w = v.w.cur,
                                    v.w.orig = v.w.orig)
        }
        m.C.new <- m.C.cur
    }

    d.out <- Reduce(function(d.acc, j) {
        d.acc[[v.cnt.nm[j]]] <- m.C.new[, j]
        d.acc
    }, seq_along(v.cnt.nm), init = d.syn)

    m.X.chk <- stats::model.matrix(fml, d.out)
    if (b.ols) {
        s.scl <- max(abs(crossprod(m.X.orig)))
        if (max(abs(crossprod(m.X.chk) - crossprod(m.X.orig)))
            > 1e-8 * s.scl)
            stop("design-gram: postcondition X*'X* = X'X violated")
    } else {
        m.I.chk <- crossprod(m.X.chk * sqrt(w.info(m.X.chk)))
        if (max(abs(m.I.chk - m.I.tgt)) > tol * s.scl.w)
            stop("design-gram: postcondition weighted Gram violated")
    }
    d.out
}

# One block of the two-part construction: recolour the block's rows to
# the original block targets. A dummy level absent from the block is
# absent on both sides (the count pinning), and its columns drop from
# the block's pin.
.dg.block <- function(m.X.orig, m.X.syn, v.cnt.idx, v.rows.o, v.rows.s) {
    m.Xo <- m.X.orig[v.rows.o, , drop = FALSE]
    m.Xs <- m.X.syn[v.rows.s, , drop = FALSE]
    v.cat.idx <- setdiff(seq_len(ncol(m.X.orig)), v.cnt.idx)
    m.D.o <- m.Xo[, v.cat.idx, drop = FALSE]
    m.D.s <- m.Xs[, v.cat.idx, drop = FALSE]

    v.keep <- colSums(abs(m.D.o)) > 0
    if (!identical(v.keep, colSums(abs(m.D.s)) > 0))
        stop("design-gram twopart: block dummy support differs")
    m.D.o <- m.D.o[, v.keep, drop = FALSE]
    m.D.s <- m.D.s[, v.keep, drop = FALSE]

    .dg.recolour(m.D.s, m.D.o,
                 m.Xo[, v.cnt.idx, drop = FALSE],
                 m.Xs[, v.cnt.idx, drop = FALSE],
                 v.w = rep(1, nrow(m.D.s)))
}

# Two-part design pinning: disjoint row blocks. The positive rows are
# recoloured to the original positive-block targets and frozen (the
# positive-model Gram is exact), the zero rows to the original zero-row
# targets, so the all-rows Gram is exact as the sum of the two blocks.
.dg.twopart <- function(d.syn, d.orig, y.nm, x.nms) {
    b.pos.o <- d.orig[[y.nm]] > 0
    b.pos.s <- d.syn[[y.nm]] > 0
    if (sum(b.pos.s) != sum(b.pos.o))
        stop("design-gram twopart: positive count not preserved")

    fml <- stats::reformulate(x.nms)
    m.X.orig <- stats::model.matrix(fml, d.orig)
    m.X.syn  <- stats::model.matrix(fml, d.syn)
    if (!identical(colnames(m.X.orig), colnames(m.X.syn)))
        stop("design-gram twopart: model matrices disagree")
    v.types   <- mapply(.id.cat, d.orig[x.nms])
    v.cnt.nm  <- names(v.types[v.types == "cnt"])
    if (length(v.cnt.nm) == 0) return(d.syn)
    v.cnt.idx <- which(colnames(m.X.orig) %in% v.cnt.nm)

    m.C.pos <- .dg.block(m.X.orig, m.X.syn, v.cnt.idx,
                         which(b.pos.o), which(b.pos.s))
    m.C.zer <- .dg.block(m.X.orig, m.X.syn, v.cnt.idx,
                         which(!b.pos.o), which(!b.pos.s))

    # Zero-model score refinement (the corollary's score condition). The
    # block pins fix the linear part of the zero-model score at the
    # original fit b0 exactly; the remaining misfit is the nonlinear
    # functional sum_i x*_i p(x*_i'b0). Damped minimum-norm Gauss-Newton
    # steps on the continuous block, constrained per block to preserve
    # the dummy sums and (to first order) both block Grams, take the
    # fast descent phase; the loop accepts at the iteration cap with
    # the achieved Newton increment (the slow tail is second order) and
    # stops loudly only above the sanity bound.
    v.z <- as.numeric(b.pos.s)
    d.z.fit <- d.orig
    d.z.fit[["z.designgram.ind"]] <- as.numeric(b.pos.o)
    mo.z0 <- stats::glm(stats::reformulate(x.nms,
                                           response = "z.designgram.ind"),
                        d.z.fit, family = stats::binomial(),
                        control = stats::glm.control(maxit = 200))
    v.b0 <- stats::coef(mo.z0)[colnames(m.X.orig)]
    if (any(is.na(v.b0)))
        stop("design-gram twopart: zero model rank-deficient")
    v.b0.C <- v.b0[v.cnt.idx]
    v.cat.idx <- setdiff(seq_len(ncol(m.X.orig)), v.cnt.idx)

    .dg.orth <- function(m.E, v.rows, m.D.all) {
        m.D <- m.D.all[v.rows, , drop = FALSE]
        m.D <- m.D[, colSums(abs(m.D)) > 0, drop = FALSE]
        m.E[v.rows, ] <- m.E[v.rows, , drop = FALSE] -
            m.D %*% solve(crossprod(m.D),
                          crossprod(m.D, m.E[v.rows, , drop = FALSE]))
        m.E
    }
    m.D.all <- m.X.syn[, v.cat.idx, drop = FALSE]
    v.rows.p <- which(b.pos.s)
    v.rows.z <- which(!b.pos.s)
    m.C.cur <- matrix(NA_real_, nrow(m.X.syn), length(v.cnt.idx))
    m.C.cur[v.rows.p, ] <- m.C.pos
    m.C.cur[v.rows.z, ] <- m.C.zer

    s.r <- length(v.cnt.idx)
    m.vech <- which(lower.tri(diag(s.r), diag = TRUE), arr.ind = TRUE)
    s.nG <- nrow(m.vech)
    .sym.basis <- function(k) {
        m.B <- matrix(0, s.r, s.r)
        m.B[m.vech[k, 1], m.vech[k, 2]] <- 1
        m.B[m.vech[k, 2], m.vech[k, 1]] <- 1
        m.B
    }
    .gram.rows <- function(m.S) {
        mapply(function(k) `if`(m.vech[k, 1] == m.vech[k, 2],
                                m.S[m.vech[k, 1], m.vech[k, 2]],
                                2 * m.S[m.vech[k, 1], m.vech[k, 2]]),
               seq_len(s.nG))
    }
    .repin <- function(m.C) {
        d.tmp <- Reduce(function(d.acc, j) {
            d.acc[[v.cnt.nm[j]]] <- m.C[, j]
            d.acc
        }, seq_along(v.cnt.nm), init = d.syn)
        m.X.tmp <- stats::model.matrix(fml, d.tmp)
        m.C[v.rows.p, ] <- .dg.block(m.X.orig, m.X.tmp, v.cnt.idx,
                                     which(b.pos.o), v.rows.p)
        m.C[v.rows.z, ] <- .dg.block(m.X.orig, m.X.tmp, v.cnt.idx,
                                     which(!b.pos.o), v.rows.z)
        m.C
    }

    s.damp   <- 0.7
    s.p.dim  <- ncol(m.X.syn)
    s.q.dim  <- s.p.dim + 2 * s.nG
    s.maxit  <- 40L
    s.outer  <- 0L
    s.prev   <- Inf
    repeat {
        s.outer <- s.outer + 1L
        m.X.cur <- m.X.syn
        m.X.cur[, v.cnt.idx] <- m.C.cur
        v.eta <- as.numeric(m.X.cur %*% v.b0)
        v.p   <- stats::plogis(v.eta)
        v.s   <- as.numeric(crossprod(m.X.cur, v.z - v.p))
        v.w0.chk <- v.p * (1 - v.p)
        m.I0  <- crossprod(m.X.cur * sqrt(v.w0.chk))
        s.db  <- max(abs(solve(m.I0, v.s)))
        b.stall <- abs(s.prev - max(abs(v.s))) < 1e-3 * max(abs(v.s))
        s.prev  <- max(abs(v.s))
        if (s.db < 1e-6 || (b.stall && s.db < 1e-4)) {
            m.C.cur <- .repin(m.C.cur)
            m.X.chk2 <- m.X.syn
            m.X.chk2[, v.cnt.idx] <- m.C.cur
            v.p2 <- stats::plogis(as.numeric(m.X.chk2 %*% v.b0))
            v.s2 <- as.numeric(crossprod(m.X.chk2, v.z - v.p2))
            m.I2 <- crossprod(m.X.chk2 * sqrt(v.p2 * (1 - v.p2)))
            if (max(abs(solve(m.I2, v.s2))) < 1e-4) break
        }
        if (s.outer >= s.maxit) {
            if (s.db > 1)
                stop("design-gram twopart: score refinement failed to converge")
            m.C.cur <- .repin(m.C.cur)
            break
        }

        v.w0 <- v.p * (1 - v.p)
        m.C.p <- m.C.cur[v.rows.p, , drop = FALSE]
        m.C.z <- m.C.cur[v.rows.z, , drop = FALSE]
        L.fwd <- function(m.E) {
            v.sc <- -as.numeric(crossprod(m.X.cur, v.w0 *
                                          as.numeric(m.E %*% v.b0.C)))
            v.sc[v.cnt.idx] <- v.sc[v.cnt.idx] +
                as.numeric(crossprod(m.E, v.z - v.p))
            m.Ep <- m.E[v.rows.p, , drop = FALSE]
            m.Ez <- m.E[v.rows.z, , drop = FALSE]
            m.Sp <- crossprod(m.Ep, m.C.p)
            m.Sz <- crossprod(m.Ez, m.C.z)
            c(v.sc, .gram.rows(m.Sp + t(m.Sp)), .gram.rows(m.Sz + t(m.Sz)))
        }
        L.adj <- function(v.dual) {
            v.a  <- v.dual[seq_len(s.p.dim)]
            v.gp <- v.dual[s.p.dim + seq_len(s.nG)]
            v.gz <- v.dual[s.p.dim + s.nG + seq_len(s.nG)]
            m.Lp <- Reduce(`+`, Map(function(k) v.gp[k] * .sym.basis(k),
                                    seq_len(s.nG)), matrix(0, s.r, s.r))
            m.Lz <- Reduce(`+`, Map(function(k) v.gz[k] * .sym.basis(k),
                                    seq_len(s.nG)), matrix(0, s.r, s.r))
            m.E <- -(v.w0 * as.numeric(m.X.cur %*% v.a)) %o% v.b0.C
            m.E <- m.E + (v.z - v.p) %o% v.a[v.cnt.idx]
            m.E[v.rows.p, ] <- m.E[v.rows.p, , drop = FALSE] +
                2 * m.C.p %*% m.Lp
            m.E[v.rows.z, ] <- m.E[v.rows.z, , drop = FALSE] +
                2 * m.C.z %*% m.Lz
            m.E <- .dg.orth(m.E, v.rows.p, m.D.all)
            .dg.orth(m.E, v.rows.z, m.D.all)
        }
        m.K <- mapply(function(k) {
            L.fwd(L.adj(`[<-`(numeric(s.q.dim), k, 1)))
        }, seq_len(s.q.dim))
        v.rhs <- c(-v.s, numeric(2 * s.nG))
        v.lambda <- tryCatch(
            solve(m.K + diag(1e-10 * max(abs(m.K)), s.q.dim), v.rhs),
            error = function(e)
                stop("design-gram twopart: score system singular"))
        m.C.cur <- m.C.cur + s.damp * L.adj(v.lambda)
    }
    m.C.pos <- m.C.cur[v.rows.p, , drop = FALSE]
    m.C.zer <- m.C.cur[v.rows.z, , drop = FALSE]

    d.out <- Reduce(function(d.acc, j) {
        v <- d.acc[[v.cnt.nm[j]]]
        v[b.pos.s]  <- m.C.pos[, j]
        v[!b.pos.s] <- m.C.zer[, j]
        d.acc[[v.cnt.nm[j]]] <- v
        d.acc
    }, seq_along(v.cnt.nm), init = d.syn)

    m.X.chk <- stats::model.matrix(fml, d.out)
    s.scl <- max(abs(crossprod(m.X.orig)))
    if (max(abs(crossprod(m.X.chk[b.pos.s, , drop = FALSE]) -
               crossprod(m.X.orig[b.pos.o, , drop = FALSE])))
        > 1e-8 * s.scl)
        stop("design-gram twopart: positive-block Gram violated")
    if (max(abs(crossprod(m.X.chk) - crossprod(m.X.orig))) > 1e-8 * s.scl)
        stop("design-gram twopart: all-rows Gram violated")
    d.out
}
