
# I will assume that the loglikelihood is calculated correctly. Its really hard
# to mess that one up.

test_that("gradient calculation works", {
  n <- 1000
  p <- 10
  t <- 20
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  beta <- rnorm(p)
  y <- rbinom(n, size = 1, p = exp(X %*% beta) / (1 + exp(X %*% beta)))
  ll <- make_logistic_loglikelihood(X, y, diagmat(beta^2), lambda = 0.05)

  random_coef <- 1:t %>% purrr::map(.f = function(x) rnorm(p))

  for (i in 1:t){
    expect_equal(
      all(abs(numDeriv::grad(ll$objective, random_coef[[i]])-ll$grad(random_coef[[i]])) < 1e-4), #Precision is about 1e-6
      TRUE
      )
  }

})
