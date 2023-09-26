test_that("Rejection sampler produces the correct distribution", {
  set.seed(0)
  enve <- LogLinearEnvelope(get_rv("n"))
  sampler <- rejection_sampler_factory(enve)
  sim <- sampler(5000)
  set.seed(NULL)
  expect_equal(shapiro.test(sim)$p > 0.05, TRUE)
})
