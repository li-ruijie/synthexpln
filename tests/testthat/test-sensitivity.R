test_that("sensitivity.local returns a scalar bound >= 0", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  fit <- lm(y ~ x, d)
  ls <- sensitivity.local(fit, bound = "dfbeta")
  expect_length(ls, 1)
  expect_gte(ls, 0)
})

test_that("sensitivity.local bound is attained for leave-one-out", {
  set.seed(1)
  d <- data.frame(y = rnorm(50), x = rnorm(50))
  fit <- lm(y ~ x, d)
  ls <- sensitivity.local(fit, bound = "dfbeta")
  # For OLS the two bounds are mathematically identical: dfbeta_i is exactly
  # (X'X)^-1 x_i r_i / (1 - h_i), which is the leverage formula. The two code paths differ
  # only by floating-point rounding whose sign is platform-dependent (dfbeta < leverage on
  # x86, dfbeta > leverage on ARM), so neither strict inequality is portable; assert
  # equality up to tolerance. The dfbeta path is smaller only for GLMs.
  expect_equal(ls, sensitivity.local(fit, bound = "leverage"), tolerance = 1e-6)
})

test_that("sensitivity.local works for GLM log-link", {
  set.seed(1)
  d <- data.frame(y = rpois(100, 2), x = rnorm(100))
  fit <- glm(y ~ x, family = poisson("log"), d)
  ls <- sensitivity.local(fit, bound = "dfbeta")
  expect_gt(ls, 0)
})

test_that(".gs.sigma2.public bounds and attains the worst-case sigma2 change", {
  set.seed(11)
  n <- 30
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  X[1, 2] <- 8                       # leverage row pushes K_X above 1
  B.y <- 4
  p <- ncol(X)
  H <- X %*% solve(crossprod(X)) %*% t(X)
  v.h <- diag(H)
  v.s <- rowSums(abs(H)) - abs(v.h)
  expect_gt(max(pmax(1 - v.h, v.s)), 1)   # design where the naive bound fails

  # brute-force worst case over the adversarial sign family, actual refits
  s2.of <- function(v) sum((v - H %*% v)^2) / (n - p)
  worst <- 0
  for (i in seq_len(n)) {
    v.adv <- -B.y * sign(H[i, ])
    v.adv[H[i, ] == 0] <- B.y
    y1 <- v.adv; y1[i] <- B.y
    y2 <- v.adv; y2[i] <- -B.y
    worst <- max(worst, abs(s2.of(y1) - s2.of(y2)))
  }
  gs <- synthexpln:::.gs.sigma2.public(X, B.y)
  expect_lte(worst, gs + 1e-10)            # valid upper bound
  expect_gte(worst, gs - 1e-8)             # attained: max s_i > 1 >= max(1 - h)
  expect_gt(worst, 4 * B.y^2 / (n - p))    # naive pre-K_X constant under-states
})

test_that(".gs.sigma2.ridge bounds the ridge worst case and has the OLS limit", {
  set.seed(12)
  n <- 30
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  X[1, 2] <- 8
  B.y <- 4
  p <- ncol(X)
  lambda <- 1
  A.inv <- solve(crossprod(X) + lambda * diag(p))
  M <- crossprod(diag(n) - X %*% A.inv %*% t(X))   # (I - H_l)'(I - H_l)
  v.m <- diag(M)
  v.s <- rowSums(abs(M)) - abs(v.m)
  expect_gt(max(pmax(v.m, v.s)), 1)

  # brute-force worst case via actual ridge refits
  s2r.of <- function(v) {
    b <- as.numeric(A.inv %*% crossprod(X, v))
    sum((v - X %*% b)^2) / (n - p)
  }
  worst <- 0
  for (i in seq_len(n)) {
    v.adv <- -B.y * sign(M[i, ])
    v.adv[M[i, ] == 0] <- B.y
    y1 <- v.adv; y1[i] <- B.y
    y2 <- v.adv; y2[i] <- -B.y
    worst <- max(worst, abs(s2r.of(y1) - s2r.of(y2)))
  }
  gs <- synthexpln:::.gs.sigma2.ridge(X, lambda, B.y)
  expect_lte(worst, gs + 1e-10)
  expect_gte(worst, gs - 1e-8)             # attained: max s_i > 1 >= max M_ii
  expect_gt(worst, 4 * B.y^2 / (n - p))

  # lambda -> 0 recovers the OLS K_X constant
  expect_equal(synthexpln:::.gs.sigma2.ridge(X, 1e-10, B.y),
               synthexpln:::.gs.sigma2.public(X, B.y),
               tolerance = 1e-6)
})
