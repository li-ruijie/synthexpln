test_that("print.synthexpln shows variant and epsilon", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  syn <- suppressWarnings(
    gen.syn.dp.projected(d, y ~ x, epsilon = 10, delta = 1e-6))
  out <- capture.output(print(syn))
  expect_true(any(grepl("V1", out)))
  expect_true(any(grepl("epsilon", out) | grepl("ε", out)))
})

test_that("summary.synthexpln returns coefficient deviations", {
  set.seed(1)
  d <- data.frame(y = rnorm(100), x = rnorm(100))
  syn <- projection(d, y ~ x)
  s <- summary(syn)
  expect_s3_class(s, "summary.synthexpln")
  # projection() output has $beta but no $beta.hat / $beta.dp, so
  # max.abs.pct is exactly 0 and the coef table has term + beta columns.
  expect_equal(s$max.abs.pct, 0)
  expect_s3_class(s$coef, "data.frame")
  expect_true(all(c("term", "beta") %in% names(s$coef)))
  expect_equal(nrow(s$coef), length(coef(lm(y ~ x, d))))
})

test_that("summary.synthexpln populates deviation table for V3 (ridge)", {
  # regression test for the previous slot-name drift: V3 exposes
  # $beta.hat as an alias for $beta.ridge so summary.synthexpln's
  # deviation branch activates for V3 too.
  set.seed(1)
  d <- data.frame(y = rnorm(200), x1 = rnorm(200), x2 = rnorm(200))
  # y ~ N(0, 1) by construction, so six sd bound it: B.y = 6. A design constant
  # of the DGP on the line above, not a statistic of d$y.
  syn <- gen.syn.dp.ridge.public(d, y ~ x1 + x2, lambda = 1,
                                  epsilon = 10, delta = 1e-6, B.y = 6,
                                  seed = 1)
  s <- summary(syn)
  expect_true("pct.dev" %in% names(s$coef))
  expect_true("beta.hat" %in% names(s$coef))
})
