



test_that("sgd converges to the right values ", {
  set.seed(0)
  targets <- rnorm(10)
  init_param <- rnorm(10)
  sc <- stopping_criterion(maxiter = 250)
  lr <- polynomial_schedule(0.4, 0.05, K=100, p=1)
  opt_target <- optimizable_parabola(targets)
  trace_vanilla <- SGD(
    opt_target, Vanilla_Optimizer(lr),
    init_param = init_param, stop_crit = sc
    )
  trace_momentum <- SGD(
    opt_target, Momentum_Optimizer(lr, beta_1 = 0.8),
    init_param = init_param, stop_crit = sc
  )
  trace_adam <- SGD(
    opt_target, Adam_Optimizer(lr, beta_1 = 0.8),
    init_param = init_param, stop_crit = sc
  )
  expect_equal(max(abs(tail(trace_vanilla, 1) - targets)) < 1e-6, TRUE)
  expect_equal(max(abs(tail(trace_momentum, 1) - targets)) < 1e-6, TRUE)
  expect_equal(max(abs(tail(trace_adam, 1) - targets)) < 1e-6, TRUE)


  #trace_vanilla %>% plot()
  #trace_momentum %>% plot()
  #trace_adam %>% plot()
  set.seed(NULL)
})




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
