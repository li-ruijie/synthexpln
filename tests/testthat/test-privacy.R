test_that("privacy.metrics computes DCR without errors", {
  set.seed(1)
  orig <- data.frame(a = rnorm(100), b = rnorm(100))
  syn  <- data.frame(a = rnorm(100), b = rnorm(100))
  m <- privacy.metrics(syn, orig, methods = "dcr")
  expect_true(is.numeric(m$dcr))
  expect_gt(m$dcr, 0)
})

test_that("privacy.metrics reports hit rate", {
  set.seed(1)
  orig <- data.frame(a = rnorm(100), b = rnorm(100))
  syn  <- orig  # perfect memorisation
  m <- privacy.metrics(syn, orig, methods = "hit")
  expect_equal(m$hit.rate, 1, tolerance = 1e-6)
})
