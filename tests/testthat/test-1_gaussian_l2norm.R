test_that("The integral is correct for a singleton", {
  expect_equal(gaussian_l2norm_matrix(c(0,0,0), 1, debug=F), 3 / 8 / sqrt(pi))
})

test_that("The integral is correct for multiple copies of the same value", {
  expect_equal(gaussian_l2norm_matrix(rep(0,10), 1), 3 / 8 / sqrt(pi))
})
