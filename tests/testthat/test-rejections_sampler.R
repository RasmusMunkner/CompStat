test_that("Rejection sampler produces the correct distribution", {
  set.seed(0)
  enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
  sampler <- rejection_sample_factory(enve)
  sim <- sampler(5000)
  set.seed(NULL)
  expect_equal(shapiro.test(sim)$p > 0.05, TRUE)
})
