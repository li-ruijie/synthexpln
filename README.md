# synthexpln

<!-- badges: start -->
[![R-CMD-check](https://github.com/li-ruijie/synthexpln/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/li-ruijie/synthexpln/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

An R package for synthetic data that reproduces your regression coefficients.

The name is short for **synth**etic data for **expl**a**n**atory models.

## What it does

You have sensitive data you cannot share, and a regression that others should be able to check and
build on. This package builds a synthetic dataset with a specific property. Fit the same model on
it and you get back the same coefficients. How close depends on the family. It is machine
precision for a Gaussian response with the identity link, about 1e-7 for other GLM families, and
about 2.3% for Poisson, where integer rounding sets the limit. Under differential privacy, noise
is added to the coefficients first, and the synthetic data matches the published noisy ones.

The guarantee covers the coefficients. The package does not try to make the synthetic data
resemble the original in its distribution, which is a separate goal.

## Install

Requires R 4.1.0 or later.

```` r
# install.packages("synthexpln")  # once on CRAN
remotes::install_github("li-ruijie/synthexpln")
````

## Quick example

```` r
library(synthexpln)
d <- data.frame(y = rnorm(500), x = rnorm(500))
# B.y is a bound on |y| that you must know before looking at the data. Here y is
# standard normal by construction, so 6 covers it. Do not compute a bound from
# d$y, since that reads the very data you are protecting.
syn <- dp.project(d, y ~ x, epsilon = 10, x.public = TRUE, B.y = 6)
coef(lm(y ~ x, syn$syn))  # equals syn$beta.dp to double-precision machine epsilon
````

## Routes

A route is one way of producing the private release. Which routes are open to you depends on what
you can state publicly before looking at the data, and `dp.project()` picks one from the arguments
you supply. There are three.

- Your predictors are public and you can bound the response. The most efficient of the three.
- Your predictors are public and you cannot bound the response. A public range for the
  coefficients stands in for the bound.
- Your predictors are sensitive and the response admits a public cutoff.

## Why it refuses

Call `dp.project()` with no public bound of any kind and it stops with an error, naming what each
route would need. Every sound route needs one quantity that was fixed before anyone looked at the
data. If you cannot say where a bound came from, it was read off the data, and a bound read off the
data carries no privacy guarantee.

## Where to go next

Three vignettes:

- `vignette("getting-started", package = "synthexpln")` to start.
- `vignette("dp-guide", package = "synthexpln")` to choose among the three sound routes.
- `vignette("calibrated-inference", package = "synthexpln")` for standard errors that account for
  the DP noise.

## License

GNU Affero General Public License v3.0 (see `LICENSE`).
Copyright (c) 2026 Ruijie Li.
