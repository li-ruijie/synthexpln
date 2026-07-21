test_that(".fit.copula returns correlation matrix and preserves column order", {
  set.seed(1)
  d <- data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50))
  cop <- synthexpln:::.fit.copula(d)
  expect_s4_class(cop$copula, "normalCopula")
  expect_equal(ncol(cop$R), 3)
  expect_equal(rownames(cop$R), c("a", "b", "c"))
})

test_that(".cf.quantile reproduces marginals up to 4th order", {
  set.seed(1)
  x <- rgamma(500, shape = 2, rate = 1)
  q <- synthexpln:::.cf.quantile(x, order = 4)
  xs <- q(ppoints(500))
  expect_lt(abs(mean(xs) - mean(x)), 0.1)
  expect_lt(abs(sd(xs) - sd(x)), 0.1)
})
