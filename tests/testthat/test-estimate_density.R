
# iter_bw_est
# We dont test validity of output, because it is really hard to test

test_that("Output has the correct length", {
  expect_length(iter_bw_est(rnorm(1000)), 2)
})

# compile_density

test_that("Density compilaion works", {
  set.seed(0)
  x <- rnorm(1000)
  dens_hat <- compile_density(0, x, 0.1)
  set.seed(NULL)
  expect_equal(dens_hat, 0.412944886551152)
})

# get_kernel

test_that("get_kernel returns the right functions", {
  f <- get_kernel("e")$kernel
  x <- c(seq(-2,1, 0.01))
  expect_equal(f(x), 3/4*(1-x^2) * (abs(x) < 1))
})








