
test_that("Envelope evaluation works", {
  enve <- LogLinearEnvelope(named_rv("n"), tangent_points = c(-1,1))
  expect_equal(enve$base_rv$f(0), dnorm(0))
})

test_that("Tangent point selection works", {
  enve <- LogLinearEnvelope(named_rv("n"), autoselection_msg = F)
  expect_equal(enve$base_rv$f(0), dnorm(0))
})

test_that("An error is thrown for invalid tangent points",{
  expect_error(LogLinearEnvelope(
    named_rv("g"),
    c(-1, 0.5, 2, 4.5)
  ))
})

test_that("Quantile function is increasing",{
  enve <- LogLinearEnvelope(
    named_rv("g"),
    c(0.5, 2, 4.5)
  )
  expect_equal(sum(diff(enve$quant(1:99 / 100)) == 0), 0)
})
