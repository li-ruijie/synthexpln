#' Unified DP projection wrapper
#'
#' Dispatches to the DP route implied by what you can actually assume, and
#' \strong{refuses} when you have assumed nothing.  The routes are selected by
#' the shape of your problem, not by a mechanism name:
#' \itemize{
#'   \item binomial family -> V5 (\code{\link{gen.syn.dp.logistic}}), needs a
#'     public coefficient box.
#'   \item \code{x.spec} supplied (the covariates are sensitive, so they must be
#'     released under DP too) -> V4 (\code{\link{gen.syn.dp.full}}), needs a
#'     public coefficient box and a public \eqn{\sigma^2} box.
#'   \item \code{x.public = TRUE} with \code{B.y} -> Route 1, exact and least noisy:
#'     V2 (\code{\link{gen.syn.dp.ols.public}}), or V3
#'     (\code{\link{gen.syn.dp.ridge.public}}) when \code{lambda} is given.
#'   \item \code{x.public = TRUE} with a coefficient box and \emph{no}
#'     \code{B.y} -> Route 2 (\code{\link{gen.syn.dp.projected.subagg}}), which
#'     needs no bound on \eqn{|y|} and so survives an unbounded response, where
#'     no finite \eqn{B_y} exists and Route 1 does not exist at all.
#' }
#'
#' Route 1 is preferred when a \code{B.y} is genuinely available, because it is
#' exact and less noisy.  Both public-covariate routes publish \eqn{X} as it stands,
#' so they sit behind \code{x.public = TRUE}, that flag being the assumption,
#' stated.
#'
#' @section Which routes are family-aware:
#' Only two of them.  Route 2 fits the \code{family} you pass, and a binomial
#' response dispatches to V5.  V2, V3, V4 and V6 release the \strong{Gaussian}
#' linear model and take no \code{family} argument at all, so this wrapper
#' \strong{refuses} rather than dispatching a non-gaussian family to one of
#' them.  Dispatching would drop the family in silence and return a Gaussian
#' refit of a model you did not ask for, which the refusal exists to prevent.
#' A non-identity link on \code{gaussian()} is refused on the same ground.
#'
#' A Poisson or Gamma response reaches DP through Route 2, which needs a
#' public design and a public coefficient box and no bound on \eqn{|y|}.
#'
#' @section Route 3 is opt-in, not inferred:
#' The third sound route, V6 (\code{\link{gen.syn.dp.suffstats}}), covers the
#' intermediate case none of the branches above match: the covariates are
#' \strong{sensitive}, as with \code{x.spec}, but the response admits a public
#' clip.  It keeps \eqn{X} private while releasing the whitened sufficient
#' statistics.  Its standing against Route 2 depends on the budget and on the
#' width of the public coefficient box, and it is ahead of the central box at
#' the smallest budget and at every budget on a structurally bounded response.
#'
#' It is reached by asking for it, \code{route = "suffstat"}, and not by
#' inference from the supplied constants:
#' \preformatted{dp.project(d, formula, epsilon, route = "suffstat",
#'            W = , B.x = , y.centre = , y.scale = , y.clip = )}
#' Every one of \eqn{(W, B_x, y.clip)} is required and must be public, on the
#' same footing as \code{B.y} everywhere else here.  They are forwarded, not
#' defaulted, so a missing one reaches V6's own refusal.
#'
#' The selector is opt-in rather than automatic because choosing Route 3 over
#' Route 2 or the \code{x.spec} route is \strong{budget-dependent}, not implied
#' by the shape of the problem: Route 3 leads at small budgets and Route 2
#' overtakes it from \eqn{\varepsilon \ge 5} where the public bound is a clip
#' on an unbounded-support response, while
#' Route 3 stays ahead at every budget only where the response is structurally
#' bounded.  A wrapper that dispatched on the constants alone would be picking a
#' privacy-utility trade-off for you, which is the choice the 0.2.0 refusal
#' design exists to leave with the caller.  For the tuning knobs V6 exposes
#' beyond these constants (\code{mu.split}, \code{lambda}, \code{eig.split}),
#' call \code{\link{gen.syn.dp.suffstats}} directly.
#'
#' @section Why there is no default route:
#' Before version 0.2.0 this function fell through to V1
#' (\code{\link{gen.syn.dp.projected}}) whenever nothing else matched, so
#' \code{dp.project(d, y ~ x, epsilon = 10)} returned a release, silently, from
#' the retired analysis-scope route that carries \strong{no formal DP guarantee}.
#' That was the package's headline entry point and its documented example.
#'
#' Every sound route rests on something fixed before the data are seen: a public
#' bound on the response, or a public box on the coefficients.  A wrapper cannot
#' invent one.  Guessing (\code{B.y = 1.5 * max(abs(y))}, say) is the
#' defect: a bound read off the sample makes the noise scale a function of the
#' private data, no sensitivity bound holds across neighbouring datasets, and the
#' release is not private at any epsilon.  Not weakly private.  Not private.
#'
#' The wrapper therefore refuses and names what each route would need.  The refusal
#' \emph{is} the result.
#'
#' @param d Data frame.
#' @param formula Model formula.
#' @param epsilon Privacy budget.
#' @param delta Privacy parameter.
#' @param family GLM family (default \code{gaussian()}).  Honoured by Route 2,
#'   and by the binomial dispatch to V5.  Every other route is Gaussian and
#'   refuses a non-gaussian family rather than dropping it: see \emph{Which
#'   routes are family-aware}.
#' @param route Which route to take.  \code{"auto"} (default) dispatches on the
#'   shape of the problem exactly as described above.  \code{"suffstat"} selects
#'   Route 3 (V6, \code{\link{gen.syn.dp.suffstats}}) explicitly and forwards
#'   \code{W}, \code{B.x}, \code{y.centre}, \code{y.scale} and \code{y.clip}.
#'   There is no value that infers Route 3, by design.
#' @param W,B.x Route 3 only: the public whitening matrix and the public bound on
#'   \eqn{\|W x\|} over the support box.  Forwarded to
#'   \code{\link{gen.syn.dp.suffstats}}, which refuses a NULL rather than
#'   inventing one.
#' @param y.centre,y.scale,y.clip Route 3 only: the public centring, scaling and
#'   clip for the response, forwarded to \code{\link{gen.syn.dp.suffstats}}.
#'   Like \code{B.y}, they must come from the domain and not from the sample.
#' @param x.public Logical; treat covariates as public.
#' @param lambda Ridge penalty (triggers V3 when combined with x.public = TRUE).
#' @param x.spec Named list of DP marginal specs (triggers V4).
#' @param score.method Copula score transform for V4/V5 when the copula is
#'   active: \code{"adaptive"} (default), \code{"normal"}, or \code{"box"}
#'   (legacy).
#' @param B.y Public upper bound on \eqn{|y|}, forwarded to V2 and V3.  It must
#'   come from the domain and not from the sample: see
#'   \code{\link{gen.syn.dp.ols.public}} for why a data-dependent bound voids the
#'   guarantee outright.
#' @param clip.lo,clip.hi public per-coefficient bounds, forwarded to the
#'   subsample-and-aggregate routes (V4, V5, Route 2).
#' @param disp.lo,disp.hi public bounds on \eqn{\sigma^2}, forwarded to V4.
#' @param dispersion A \code{dp.dispersion} object, forwarded to Route 2.
#' @param seed Optional RNG seed.
#' @return A \code{synthexpln} object as returned by the selected variant.
#' @seealso \code{\link{gen.syn.dp.projected.subagg}} for the route that needs
#'   the least.
#' @export
#' @examples
#' d <- data.frame(y = rnorm(100), x = rnorm(100))
#' # A public bound on |y| buys the exact public-design route.  Here y is
#' # standard normal by construction, so 6 sd bounds it: a fact about the
#' # example's DGP, not a quantity read off d$y.
#' syn <- dp.project(d, y ~ x, epsilon = 10, x.public = TRUE, B.y = 6)
#' coef(lm(y ~ x, syn$syn))   # == syn$beta.dp
#'
#' # With no public bound of any kind, the wrapper refuses rather than
#' # dispatching to a route that carries no guarantee:
#' try(dp.project(d, y ~ x, epsilon = 10))
dp.project <- function(d, formula, epsilon, delta = 1e-6,
                       family = stats::gaussian(),
                       route = c("auto", "suffstat"),
                       x.public = FALSE, lambda = NULL, x.spec = NULL,
                       score.method = "adaptive", B.y = NULL,
                       clip.lo = NULL, clip.hi = NULL,
                       disp.lo = NULL, disp.hi = NULL,
                       W = NULL, B.x = NULL,
                       y.centre = NULL, y.scale = NULL, y.clip = NULL,
                       dispersion = NULL, seed = NULL) {
    route <- match.arg(route)

    # Normalise the family name the same way as projection() and the V1 variant
    # so that future extensions dispatch on a single canonical lowercase key.
    # Hoisted above every branch so the opt-in Route 3 below is covered by the
    # same family guard as the shape-based routes.
    fam.nm.raw <- family[["family"]]
    fam.nm <- switch(fam.nm.raw,
        Gamma                = "gamma",
        `Negative Binomial`  = "neg.binomial",
        `inverse.gaussian`   = "inv.gaussian",
        `Inverse Gaussian`   = "inv.gaussian",
        tolower(fam.nm.raw))

    # V2, V3, V4 and V6 are the Gaussian linear-model releases.  Not one of them
    # takes a `family` argument, so a non-gaussian family handed to them is
    # dropped in silence and the caller receives a Gaussian refit of a model
    # they did not ask for.  Route 2 is the family-general route and V5 owns
    # binomial, so those two are the only ones that may see a general family.
    #
    # The analysis tree's own wrapper has refused here since it was written
    # ("Gamma/Poisson only supported under analysis scope",
    # lib/28-dp-projection.r), and this rewrite lost the guard.
    gaussian.only <- function(route.nm) {
        if (fam.nm == "gaussian" && family[["link"]] == "identity")
            return(invisible(NULL))
        stop(sprintf("family = %s(link = '%s') is not supported by %s.\n",
                     fam.nm.raw, family[["link"]], route.nm),
             "  That route releases the GAUSSIAN linear model and takes no\n",
             "  family argument at all, so dispatching to it would drop yours\n",
             "  and hand you a refit of a model you did not ask for.\n",
             "\n",
             "  Two routes are family-aware:\n",
             "    Route 2 fits your family directly, given a PUBLIC design and\n",
             "      a PUBLIC coefficient box:\n",
             "        dp.project(..., x.public = TRUE, clip.lo = , clip.hi = )\n",
             "    binomial responses dispatch to V5 automatically.\n",
             call. = FALSE)
    }

    # ROUTE 3, OPT-in.  Checked before every shape-based branch, so an explicit
    # request is not overridden by one of them.  The constants are forwarded,
    # not defaulted, so a missing W/B.x/y.clip reaches V6's own refusal, the
    # same discipline the shape-based branches use for B.y and the clip box.
    #
    # There is deliberately no branch that INFERS this route from the presence
    # of (W, B.x, y.clip): Route 3's lead over Route 2 reverses with the budget,
    # so inferring it would pick a privacy-utility trade-off on the caller's
    # behalf.  author decision 2026-07-30.
    if (route == "suffstat") {
        gaussian.only("Route 3 (V6, gen.syn.dp.suffstats)")
        return(gen.syn.dp.suffstats(d, formula, epsilon = epsilon,
                                    delta = delta,
                                    W = W, B.x = B.x,
                                    y.centre = y.centre, y.scale = y.scale,
                                    y.clip = y.clip,
                                    dispersion = dispersion,
                                    x.public = x.public, seed = seed))
    }

    if (fam.nm == "binomial") {
        return(gen.syn.dp.logistic(d, formula, epsilon, delta,
                                   x.public = x.public, x.spec = x.spec,
                                   clip.lo = clip.lo, clip.hi = clip.hi,
                                   score.method = score.method, seed = seed))
    }
    # Bounds are forwarded and not defaulted.  A NULL therefore produces an
    # explicit refusal for the user rather than being suppressed here.
    #
    # sensitive covariates.  x.spec is how you say "X must itself be released
    # under DP", so it selects the route that does that.
    if (!is.null(x.spec)) {
        gaussian.only("the x.spec route (V4, gen.syn.dp.full)")
        return(gen.syn.dp.full(d, formula, x.spec = x.spec,
                               epsilon = epsilon, delta = delta,
                               clip.lo = clip.lo, clip.hi = clip.hi,
                               disp.lo = disp.lo, disp.hi = disp.hi,
                               B.y = B.y,
                               score.method = score.method, seed = seed))
    }

    # public covariates.  Both routes below publish X as it stands, which is
    # sound only under that premise, so they sit behind x.public = TRUE.
    if (x.public) {
        # Route 1, exact and least noisy.  Defined THROUGH B.y, and the callees refuse
        # a NULL bound rather than inventing one.
        if (!is.null(lambda)) {
            gaussian.only("Route 1 ridge (V3, gen.syn.dp.ridge.public)")
            return(gen.syn.dp.ridge.public(d, formula, lambda = lambda,
                                           epsilon = epsilon, delta = delta,
                                           B.y = B.y, seed = seed))
        }
        if (!is.null(B.y) && fam.nm == "gaussian")
            return(gen.syn.dp.ols.public(d, formula,
                                         epsilon = epsilon, delta = delta,
                                         B.y = B.y, seed = seed))

        # Route 2.  Needs no bound on |y|: its sensitivity comes from clipping
        # the per-block coefficients into the public box, so it survives an
        # unbounded response, where no finite B.y exists and Route 1 does not
        # exist at all.  It is checked after Route 1 because Route 1 is exact
        # and less noisy when a B.y is genuinely available.
        if (!is.null(clip.lo) && !is.null(clip.hi))
            return(gen.syn.dp.projected.subagg(d, formula, epsilon = epsilon,
                                               delta = delta, family = family,
                                               clip.lo = clip.lo,
                                               clip.hi = clip.hi,
                                               dispersion = dispersion,
                                               seed = seed))

        # Public design, gaussian, and no bound of either kind: let Route 1's own
        # refusal explain what B.y is and why it cannot be defaulted.
        if (fam.nm == "gaussian")
            return(gen.syn.dp.ols.public(d, formula,
                                         epsilon = epsilon, delta = delta,
                                         B.y = B.y, seed = seed))

        # A public design under a non-gaussian family reaches here only when no
        # coefficient box was supplied, so Route 2, the one route that fits this
        # family, was unavailable.  Say that.  The generic refusal below would
        # otherwise report "no public bounds supplied" to a caller who supplied
        # B.y, which is true of no bound they gave and false of the situation.
        gaussian.only("Route 1 (V2, gen.syn.dp.ols.public)")
    }

    stop("no public bounds supplied, so no DP route is available.\n",
         "  Every sound route rests on something fixed BEFORE the data are\n",
         "  seen.  A bound read off the sample makes the noise scale a\n",
         "  function of the private data, and the release is then not private\n",
         "  at any epsilon, so pick the assumption you can actually defend:\n",
         "\n",
         "  Route 1, exact and least noisy, needs a PUBLIC design and a\n",
         "    PUBLIC bound on |y| from the domain (a physiological range, a\n",
         "    reporting cap):\n",
         "      dp.project(..., x.public = TRUE, B.y = <bound>)\n",
         "    Add lambda = <penalty> for the ridge variant.\n",
         "\n",
         "  Route 2 needs NO bound on |y| at all, only a PUBLIC box on the\n",
         "    coefficients.  Where the response is unbounded (any Gaussian\n",
         "    outcome), no finite B.y exists and this is the only route there\n",
         "    is.  It publishes X as it stands, so it too assumes x.public:\n",
         "      dp.project(..., x.public = TRUE, clip.lo = <lo>, clip.hi = <hi>,\n",
         "                 dispersion = dp.dispersion.public(s2, source = ))\n",
         "\n",
         "  SENSITIVE covariates must be released under DP too, which needs\n",
         "    x.spec, a coefficient box and a sigma^2 box:\n",
         "      dp.project(..., x.spec = , clip.lo = , clip.hi = ,\n",
         "                 disp.lo = , disp.hi = )\n",
         "\n",
         "  Route 3 keeps SENSITIVE covariates private where the response does\n",
         "    admit a public clip.  Against Route 2 it is ahead of the central\n",
         "    box at the smallest budget and at every budget on a structurally\n",
         "    bounded response.  It is OPT-IN and not inferred, since\n",
         "    choosing it over Route 2 or the x.spec route is budget-dependent.\n",
         "    Ask for it by name:\n",
         "      dp.project(..., route = \"suffstat\", W = , B.x = ,\n",
         "                 y.centre = , y.scale = , y.clip = )\n",
         "\n",
         "  gen.syn.dp.projected() (V1) carries NO formal DP guarantee and is\n",
         "  deliberately not dispatched to.  Call it directly if you want the\n",
         "  empirical-DP baseline, and read its warning.",
         call. = FALSE)
}
