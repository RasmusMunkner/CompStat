#Lazy test - Using one of the tangent points for evaluation
test_that("Envelope evaluation works", {
  enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
  expect_equal(enve$f(0), dnorm(0))
})

test_that("An error is thrown for invalid tangent points",{
  expect_error(LogLinearEnvelope(
    function(z){dgamma(z,2)},
    function(z){dgamma(z,2) / z - dgamma(z, 2)},
    c(-1, 0.5, 2, 4.5)
  ))
})

test_that("Quantile function is increasing",{
  enve <- LogLinearEnvelope(
    function(z){dgamma(z,2)},
    function(z){dgamma(z,2) / z - dgamma(z, 2)},
    c(0.5, 2, 4.5)
  )
  expect_equal(sum(diff(enve$quant(1:99 / 100)) == 0), 0)
})
