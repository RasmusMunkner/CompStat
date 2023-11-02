
test_that("sgd converges to the right values ", {
  set.seed(0)
  targets <- rnorm(10)
  init_par <- rnorm(10)
  sc <- stopping_criterion(maxiter = 250)
  lr <- polynomial_schedule(0.4, 0.05, later=100, p=1)
  opt_target <- optimizable_parabola(targets)
  trace_vanilla <- SGD(
    opt_target, Vanilla_Optimizer(lr),
    init_par = init_par, stop_crit = sc
    )
  trace_momentum <- SGD(
    opt_target, Momentum_Optimizer(lr, beta_1 = 0.8),
    init_par = init_par, stop_crit = sc
  )
  trace_adam <- SGD(
    opt_target, Adam_Optimizer(lr, beta_1 = 0.8),
    init_par = init_par, stop_crit = sc
  )
  expect_equal(max(abs(tail(trace_vanilla) - targets)) < 1e-6, TRUE)
  expect_equal(max(abs(tail(trace_momentum) - targets)) < 1e-6, TRUE)
  expect_equal(max(abs(tail(trace_adam) - targets)) < 1e-6, TRUE)

  set.seed(NULL)
})

test_that("cpp batch gradient results are consistent with r results", {

  t <- 20
  p <- 4

  sll <- simple_logistic_loglikelihood(n = 100, p = p)

  optfun_loglik <- logistic_loglikelihood(
    design = sll$X, response = sll$y,
    penalty_matrix = matrix(0, nrow = ncol(sll$X), ncol = ncol(sll$X)),
    lambda = 0
    )

  random_coef <- 1:t %>% purrr::map(.f = function(i) rnorm(p))

  for (i in 1:t){
    expect_equal(
      lll_gradC(
        sll$X, random_coef[[i]], sll$y, optfun_loglik$penalty_matrix, 0) %>%
        as.vector(),
      optfun_loglik$grad(random_coef[[i]])
      )
  }

})

test_that("C++ sgd converges to the right value", {

  set.seed(0)
  t <- 10
  p <- 8

  targets <- rnorm(p)
  sll <- simple_logistic_loglikelihood(n = 647, p = p, beta = targets)

  opt_target <- logistic_loglikelihood(
    design = sll$X, response = sll$y,
    penalty_matrix = matrix(0, nrow = ncol(sll$X), ncol = ncol(sll$X)),
    lambda = 0
  )

  random_coef <- 1:t %>% purrr::map(.f = function(i) rnorm(p))

  epochs <- 50
  batch_size <- 24
  lr <- polynomial_schedule(0.1, 0.001, 50)
  opt <- Adam_Optimizer(
    lr, batch_size, beta_1 = 0.9, beta_2 = 0.95, eps = 1e-8, amsgrad = T)

  for (i in 1:t){
    par_r <- SGD(
      opt_target, opt,
      init_par = random_coef[[i]],
      stop_crit = epochs,
      seed = 0
      ) %>% tail()

    par_c <- SGD_CPP(
      opt_target, opt,
      init_par = random_coef[[i]],
      stop_crit = epochs,
      seed = 0
    ) %>% tail()

    expect_equal(par_c, par_r)
  }

  set.seed(NULL)
})










