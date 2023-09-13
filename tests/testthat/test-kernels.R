test_that("Approximation of L2 norm runs", {
  expect_equal(estimate_l2norm(1, "e", 1), 9/4)
})
