
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

test_that("C++ sgd converges to the right value for a simple case", {

  set.seed(0)
  t <- 1
  p <- 1

  targets <- rnorm(p)
  sll <- simple_logistic_loglikelihood(n = 2, p = p, beta = targets)

  opt_target <- logistic_loglikelihood(
    design = sll$X, response = sll$y,
    penalty_matrix = matrix(0, nrow = ncol(sll$X), ncol = ncol(sll$X)),
    lambda = 0
  )

  random_coef <- 1:t %>% purrr::map(.f = function(i) rnorm(p))

  epochs <- 1
  batch_size <- 1
  lr <- polynomial_schedule(1, 0.5, 100)
  opt <- Adam_Optimizer(
    lr, batch_size = batch_size, beta_1 = 0.9, beta_2 = 0.95, eps = 1e-8, amsgrad = T)
  opt <- Vanilla_Optimizer(lr, batch_size)

  for (i in 1:t){
    par_r <- SGD(
      opt_target, opt,
      init_par = random_coef[[i]],
      stop_crit = epochs,
      seed = 0,
      debug = F
    )

    par_c <- SGD_CPP(
      opt_target, opt,
      init_par = random_coef[[i]],
      stop_crit = epochs,
      seed = 0,
      debug = F
    )

    expect_equal(object=par_c, expected=par_r)
  }

  set.seed(NULL)
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

test_that("EM algorithm works", {


  bg_theta <- list(
    p = c(0.25),
    mu = c(-10,10),
    sigma2 = c(1,9),
    nu = c(2,10)
  )
  theta <- bg_theta %>% unlist()
  theta_prime <- list(
    p = c(0.5),
    mu = c(-6,12),
    sigma2 = c(4,4),
    nu = c(2,10)
  ) %>% unlist()

  y <- rtmix(10000, bg_theta$p, bg_theta$mu, bg_theta$sigma2, bg_theta$nu)
  #y <- rtmix(50, bg_theta$p, bg_theta$mu, bg_theta$sigma2, bg_theta$nu)

  Estep <- MinimalEstep(y, init_par = theta_prime)

  plot(Estep)

  lr <- polynomial_schedule(1e-1, 1e-2, 500)
  opt <- Adam_Optimizer(lr = lr, batch_size = ceiling(sqrt(length(y))))
  sc <- stopping_criterion(100)
  trace <- SGD(
    Estep,
    opt,
    stop_crit = sc,
    init_par = c(
      c(0.5), quantile(y, probs = c(0.25, 0.75)),
      c(2, 3), c(2,10)) %>% setNames(unlist(theta) %>% names())
  )

  plot(trace)

  plot(Estep)

  optpar <- trace %>% tail() %>% magrittr::set_names(Estep$par_alloc)

  second_trace <- SGD(
    Estep,
    opt,
    stop_crit = 100,
    init_par = trace %>% tail()
  )

  plot(second_trace)

  optpar_2 <- second_trace %>% tail() %>% magrittr::set_names(Estep$par_alloc)

})










