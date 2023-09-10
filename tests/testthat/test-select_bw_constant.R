# second_derivative
test_that("second_derivatiev works", {
  f_1 <- second_derivative(exp)
  f_2 <- second_derivative(log)
  f_3 <- second_derivative(function(x) x^3 / 6)
  x <- seq(-3, 3, 0.1)
  expect_equal(f_1(x), exp(x))
  expect_equal(f_2(1+abs(x)), -1/abs(1+abs(x))^2)
  expect_equal(f_3(x), x)
})

# l2norm
test_that("l2norm works", {
  expect_equal(l2norm(dexp), 1/2)
  expect_equal(l2norm(function(x) 1/sqrt(x), lower = 1, upper = exp(1)), 1)
})

# moment
test_that("moment works", {
  expect_equal(moment(dnorm, c(1,2,3,4)), c(0,1,0,3))
  expect_equal(moment(dexp, c(1,2)), c(1,2))
  expect_equal(moment(dunif, c(1,2)), c(1/2,1/3))
})

# select_bw_constant
test_that("select_bw_constant works", {
  expect_equal(select_bw_constant(dnorm, dnorm), (4/3)^(1/5))
})
