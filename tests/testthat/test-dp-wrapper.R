test_that("dp.project REFUSES when it has no public bounds to dispatch on", {
  # This is the headline change.  The wrapper used to fall through to V1, the
  # retired analysis-scope route, whenever nothing else matched, so
  # dp.project(d, y ~ x, epsilon = 10) returned a release with no formal DP
  # guarantee, silently, from the package's most prominent entry point.
  #
  # Every sound route rests on something fixed before the data are seen.  A
  # wrapper cannot invent one, and guessing is the defect itself, so it refuses.
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))

  expect_error(dp.project(d, y ~ x, epsilon = 10, delta = 1e-6),
               "no public bounds supplied")

  # and the refusal names what each route would need: it is a signpost, not a
  # dead end
  e <- tryCatch(dp.project(d, y ~ x, epsilon = 10), error = conditionMessage)
  expect_match(e, "Route 1")
  expect_match(e, "Route 2")
  expect_match(e, "B\\.y")
  expect_match(e, "clip\\.lo")
  # it must not quietly send anyone to the route with no guarantee
  expect_match(e, "NO formal DP guarantee")
})

test_that("dp.project auto-selects V2 for public-X OLS", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  # y ~ N(0, 1) by construction, so six sd bound it: B.y = 6. A design constant
  # of the DGP on the line above, not a statistic of d$y. The wrapper forwards
  # B.y to the variant it dispatches to, and V2 requires it.
  syn <- dp.project(d, y ~ x, epsilon = 10, delta = 1e-6,
                    x.public = TRUE, B.y = 6)
  expect_equal(syn$variant, "V2")
})

test_that("dp.project auto-selects V3 when lambda specified", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))
  # y ~ N(0, 1) by construction, so B.y = 6 (six sd) from the DGP above.
  syn <- dp.project(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6,
                    x.public = TRUE, lambda = 1, B.y = 6)
  expect_equal(syn$variant, "V3")
})

test_that("dp.project auto-selects V4 when x.spec provided", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  # y is independent of x by construction, so the population coefficients are 0
  # and sigma^2 is 1.  Both boxes come from that DGP.
  syn <- dp.project(d, y ~ x, epsilon = 10, delta = 1e-6,
                    x.spec = list(x = list(type = "continuous",
                                            bounds = c(-5, 5))),
                    clip.lo = -5, clip.hi = 5, disp.lo = 0, disp.hi = 4)
  expect_equal(syn$variant, "V4")
})

test_that("dp.project auto-selects V5 for binomial family", {
  set.seed(1)
  d <- data.frame(y = rbinom(400, 1, 0.5), x = rnorm(400))
  syn <- dp.project(d, y ~ x, epsilon = 10, delta = 1e-6,
                    family = binomial("logit"), x.public = TRUE,
                    clip.lo = -3, clip.hi = 3)
  expect_equal(syn$variant, "V5")
})

test_that("dp.project routes a public design with only a coefficient box to Route 2", {
  set.seed(2)
  n <- 800
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  # A public design but no bound on |y|: y is Gaussian, hence unbounded, so no
  # finite B.y exists and Route 1 does not exist at all.  Route 2 needs no such
  # bound.  This is the route that exists precisely because the other one dies
  # here, and it is why the wrapper checks Route 1 first but does not stop there.
  syn <- dp.project(d, y ~ x1 + x2, epsilon = 10, delta = 1e-6,
                    x.public = TRUE, clip.lo = -5, clip.hi = 5,
                    dispersion = dp.dispersion.public(
                      1, source = "unit residual variance, by the DGP"),
                    seed = 1)
  expect_equal(syn$variant, "D2")
  expect_equal(unname(coef(lm(y ~ x1 + x2, syn$syn))), unname(syn$beta.dp),
               tolerance = 1e-8)
})

test_that("dp.project will not publish a SENSITIVE design verbatim", {
  # Both public-covariate routes release X as it stands.  The wrapper's default
  # is x.public = FALSE (assume the covariates are sensitive), so it must NOT
  # quietly hand a sensitive design to a route that publishes it.  Saying
  # nothing about x.public is not consent.
  set.seed(3)
  n <- 400
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)

  # a coefficient box alone does not buy a route: X would still be published
  expect_error(
    dp.project(d, y ~ x1 + x2, epsilon = 10, clip.lo = -5, clip.hi = 5,
               dispersion = dp.dispersion.public(1, source = "declared")),
    "no public bounds supplied")

  # the binomial route likewise refuses to publish a design it was told is
  # sensitive, and names the missing piece
  db <- data.frame(y = rbinom(n, 1, 0.5), x = rnorm(n))
  expect_error(
    dp.project(db, y ~ x, epsilon = 10, family = binomial("logit"),
               clip.lo = -3, clip.hi = 3),
    "x.spec")
})

# ---- Route 3, the opt-in selector ---------
#
# Before this the wrapper could not reach V6 at all, so a caller with sensitive
# covariates and a publicly clippable response, the case the whitened
# sufficient-statistic release proposition exists for, was sent to V4 instead.
# The selector is opt-in and not inferred, because choosing Route 3 over Route
# 2 is budget-dependent.

W3 <- local({
  m <- diag(3)
  dimnames(m) <- list(c("(Intercept)", "x1", "x2"),
                      c("(Intercept)", "x1", "x2"))
  m
})

mk.suff.d <- function(n = 400L) {
  set.seed(42)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
  d
}

test_that("route = 'suffstat' reaches V6 and agrees with the direct call", {
  d <- mk.suff.d()
  args <- list(B.x = sqrt(1 + 2 * 6^2), y.centre = 0, y.scale = 1,
               y.clip = c(-8, 8))

  via.wrapper <- do.call(dp.project,
    c(list(d = d, formula = y ~ x1 + x2, epsilon = 10, route = "suffstat",
           W = W3), args, list(seed = 7)))
  direct <- do.call(gen.syn.dp.suffstats,
    c(list(d = d, formula = y ~ x1 + x2, epsilon = 10, W = W3), args,
      list(seed = 7)))

  expect_equal(via.wrapper$variant, "V6")
  # the wrapper must FORWARD, not re-derive: same seed, byte-identical release
  expect_identical(via.wrapper$beta.dp, direct$beta.dp)
  expect_identical(via.wrapper$S.dp, direct$S.dp)
  expect_identical(via.wrapper$mu.channels, direct$mu.channels)
})

test_that("route = 'suffstat' forwards missing constants to V6's own refusal", {
  # The wrapper defaults nothing.  A missing W/B.x/y.clip must reach V6's
  # refusal, naming the public constant, and not on the wrapper's generic
  # "no public bounds supplied" message, which would misdescribe the problem.
  d <- mk.suff.d()
  e <- tryCatch(dp.project(d, y ~ x1 + x2, epsilon = 10, route = "suffstat"),
                error = conditionMessage)
  expect_match(e, "PUBLIC")
  expect_false(grepl("no public bounds supplied", e, fixed = TRUE))
})

test_that("Route 3 is not INFERRED from the constants alone", {
  # Supplying every Route 3 constant without asking for the route must not
  # select it.  Route 3's lead over Route 2 reverses with the budget, so
  # inferring it would pick a privacy-utility trade-off for the caller.
  d <- mk.suff.d()
  expect_error(
    dp.project(d, y ~ x1 + x2, epsilon = 10,
               W = W3, B.x = sqrt(1 + 2 * 6^2),
               y.centre = 0, y.scale = 1, y.clip = c(-8, 8)),
    "no public bounds supplied")

  # and the refusal points at the selector by name rather than at a direct call
  e <- tryCatch(dp.project(d, y ~ x1 + x2, epsilon = 10),
                error = conditionMessage)
  expect_match(e, 'route = "suffstat"', fixed = TRUE)
})

test_that("an unknown route is rejected rather than silently ignored", {
  d <- mk.suff.d()
  expect_error(dp.project(d, y ~ x1 + x2, epsilon = 10, route = "route3"))
})

# ---- the family guard ------------------------------------
#
# V2, V3, V4 and V6 release the Gaussian linear model and none of them takes a
# family argument, so before this guard a non-gaussian family handed to any of
# them was dropped in silence and the caller got a Gaussian refit of a model
# they did not ask for.  Measured at the time: family = poisson() with
# x.public = TRUE, B.y = 6, lambda = 1 returned variant "V3" and a continuous
# response.  The analysis tree's own wrapper has refused this since it was
# written ("Gamma/Poisson only supported under analysis scope",
# lib/28-dp-projection.r:1542), and this rewrite lost the guard.

mk.pois.d <- function(n = 400L) {
  set.seed(11)
  data.frame(y = rpois(n, 4), x1 = rnorm(n), x2 = rnorm(n))
}

test_that("a non-gaussian family is REFUSED by every gaussian-only route", {
  d <- mk.pois.d()
  spec <- list(x1 = list(type = "continuous", bounds = c(-5, 5)))

  # V3, the ridge branch: keyed on lambda alone, so it did not see the family.
  expect_error(
    dp.project(d, y ~ x1 + x2, epsilon = 10, family = poisson(),
               x.public = TRUE, B.y = 20, lambda = 1),
    "is not supported by Route 1 ridge")

  # V4, the x.spec branch.
  expect_error(
    dp.project(d, y ~ x1, epsilon = 10, family = poisson(), x.spec = spec,
               clip.lo = -5, clip.hi = 5, disp.lo = 0.1, disp.hi = 50),
    "is not supported by the x.spec route")

  # V6, reached only by the opt-in selector, which returned before the family
  # was even normalised until the guard was hoisted above it.
  expect_error(
    dp.project(d, y ~ x1 + x2, epsilon = 10, family = poisson(),
               route = "suffstat", B.x = 3, y.centre = 0, y.scale = 1,
               y.clip = c(-8, 8)),
    "is not supported by Route 3")

  # a non-identity link on gaussian() is the same defect wearing a family the
  # guard would otherwise wave through
  expect_error(
    dp.project(d, y ~ x1 + x2, epsilon = 10, family = gaussian(link = "log"),
               x.public = TRUE, B.y = 20, lambda = 1),
    "link = 'log'")
})

test_that("the refusal names the family-aware route rather than misdescribing", {
  # A public design and a public B.y under a poisson response used to fall all
  # the way through to "no public bounds supplied", which is false about a
  # caller who supplied B.y.  The package's own suffstat test already holds
  # refusals to this standard.
  d <- mk.pois.d()
  e <- tryCatch(
    dp.project(d, y ~ x1 + x2, epsilon = 10, family = poisson(),
               x.public = TRUE, B.y = 20),
    error = conditionMessage)

  expect_false(grepl("no public bounds supplied", e, fixed = TRUE))
  expect_match(e, "poisson")
  expect_match(e, "Route 2")
  expect_match(e, "clip\\.lo")
})

test_that("the two family-aware routes still accept a non-gaussian family", {
  # The guard must not close the routes that genuinely fit the family.
  d <- mk.pois.d()

  syn <- dp.project(d, y ~ x1 + x2, epsilon = 10, family = poisson(),
                    x.public = TRUE, clip.lo = -5, clip.hi = 5, seed = 1)
  expect_equal(syn$variant, "D2")
  expect_true(all(syn$syn$y >= 0))
  expect_equal(syn$syn$y, round(syn$syn$y))

  db <- data.frame(y = rbinom(400, 1, 0.5), x = rnorm(400))
  syn.b <- dp.project(db, y ~ x, epsilon = 10, family = binomial("logit"),
                      x.public = TRUE, clip.lo = -3, clip.hi = 3)
  expect_equal(syn.b$variant, "V5")
})

test_that("gaussian identity still reaches every route it did before", {
  # The guard is a refusal added to branches that previously returned, so the
  # regression that matters is whether the supported shapes are untouched.
  set.seed(1)
  d <- data.frame(y = rnorm(200), x1 = rnorm(200), x2 = rnorm(200))

  expect_equal(dp.project(d, y ~ x1, epsilon = 10, x.public = TRUE,
                          B.y = 6)$variant, "V2")
  expect_equal(dp.project(d, y ~ x1 + x2, epsilon = 10, x.public = TRUE,
                          lambda = 1, B.y = 6)$variant, "V3")
  expect_equal(dp.project(d, y ~ x1, epsilon = 10,
                          x.spec = list(x1 = list(type = "continuous",
                                                  bounds = c(-5, 5))),
                          clip.lo = -5, clip.hi = 5,
                          disp.lo = 0, disp.hi = 4)$variant, "V4")
})
