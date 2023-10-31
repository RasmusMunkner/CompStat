


test_that("cpp batch gradient results are consistent with r results", {

  t <- 20
  p <- 4

  sll <- simple_logistic_loglikelihood(n = 100, p = p)

  optfun_loglik <- make_logistic_loglikelihood(
    design = sll$X, response = sll$y,
    penalty_matrix = matrix(0, nrow = ncol(sll$X), ncol = ncol(sll$X)),
    lambda = 0
    )

  random_coef <- 1:t %>% purrr::map(.f = function(i) rnorm(p))

  for (i in 1:t){
    expect_equal(
      batch_gradient(sll$X, random_coef[[i]], sll$y) %>% as.vector(),
      optfun_loglik$grad(random_coef[[i]])
      )
  }

})
