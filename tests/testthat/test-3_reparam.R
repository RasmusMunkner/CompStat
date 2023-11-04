

test_that("reparametrized t-mixtures work", {

  set.seed(0)
  theta <- list(
    p = c(0.2, 0.4),
    mu = c(-3, 1, 8),
    sigma2 = c(1, 8, 2),
    nu = c(2, 3, 6)
  )
  y <- rtmix(1e3, theta$p, theta$mu, theta$sigma2, theta$nu)
  theta_newpar <- list(
    be = theta$p %>% p_() %>% inverse_softmax() %>% `[`(-length(.)),
    mu = theta$mu,
    kp = log(theta$sigma2),
    nu = theta$nu
  ) %>% unlist()
  theta <- theta %>% unlist()

  Qfunc_old <- MinimalEstep(y, theta)
  Qfunc_new <- EstepReparam(y, theta_newpar)

  expect_equal(
    Qfunc_new$f_y_by_z(theta_newpar),
    Qfunc_old$f_y_by_z(theta)
    )

  expect_equal(
    Qfunc_new$f_y(theta_newpar),
    Qfunc_old$f_y(theta)
  )

  expect_equal(
    Qfunc_new$cond_p(theta_newpar, mode = "r", cache = F),
    Qfunc_old$cond_p(theta, mode = "r")
  )

  expect_equal(
    Qfunc_new$objective(theta_newpar, cache = F),
    Qfunc_old$objective(theta)
  )

  # Note that the gradients are explicitly NOT equal,
  # since we are using different parametrizations
  # and hence the gradients are with respect to different parameters

  expect_equal(
    Qfunc_old$loglikelihood(theta),
    Qfunc_new$loglikelihood(theta_newpar)
  )

})

test_that("Gradients are correct", {

  set.seed(0)
  theta <- list(
    p = c(0.2, 0.4),
    mu = c(-3, 1, 8),
    sigma2 = c(1, 8, 2),
    nu = c(2, 3, 6)
  )
  y <- rtmix(1e3, theta$p, theta$mu, theta$sigma2, theta$nu)
  theta_newpar <- list(
    be = theta$p %>% p_() %>% inverse_softmax() %>% `[`(-length(.)),
    mu = theta$mu,
    kp = log(theta$sigma2),
    nu = theta$nu
  )
  theta_newpar_0 <- list(
    be = theta_newpar$be %>% magrittr::add(rnorm(length(.), 0, 1)),
    mu = theta_newpar$mu %>% magrittr::add(rnorm(length(.), 0, 1)),
    kp = theta_newpar$kp %>% magrittr::add(rnorm(length(.), 0, 1)),
    nu = theta_newpar$nu
  )
  t_old <- theta %>% unlist()
  t <- theta_newpar %>% unlist()
  t0 <- theta_newpar_0 %>% unlist()

  Qfunc <- EstepReparam(y, t)

  A <- numDeriv::grad(Qfunc$objective, t0)

  B <- Qfunc$grad(t0)

  expect_equal((A-B) %>% magrittr::set_names(NULL), rep(0, length(A)))

  set.seed(NULL)

})

test_that("fisher information is calculated correctly", {

  set.seed(0)
  t <- list(
    be = c(0.2, 0.4) %>% p_() %>% inverse_softmax() %>% `[`(-length(.)),
    mu = c(-3, 1, 8),
    kp = log(c(1, 8, 2)),
    nu = c(2, 3, 6)
  ) %>% unlist()
  Qfunc <- EstepReparam(0, t)
  y <- rtmix(1e4, Qfunc$get$p(t), Qfunc$get$mu(t),
             Qfunc$get$s2(t), Qfunc$get$nu(t))
  Qfunc <- EstepReparam(y, t) # Wraps the previously shown implementations

  A <- Qfunc$fisher(t, method = 1)

  B <- Qfunc$fisher(t, method = 2)

  C <- Qfunc$fisher(t, method = 3)

  expect_equal(A, B)
  expect_equal(max(abs(A - C)) < 1, TRUE) # Not exact for finite data

  set.seed(NULL)

})

test_that("EM for reparametrized mixtures works", {

  set.seed(0)
  theta <- list(
    p = c(0.2, 0.4),
    mu = c(-3, 1, 8),
    sigma2 = c(1, 8, 2),
    nu = c(2, 3, 6)
  )
  y <- rtmix(1e3, theta$p, theta$mu, theta$sigma2, theta$nu)
  theta_newpar <- list(
    be = theta$p %>% p_() %>% inverse_softmax() %>% `[`(-length(.)),
    mu = theta$mu,
    kp = log(theta$sigma2),
    nu = theta$nu
  )
  theta_newpar_0 <- list(
    be = theta_newpar$be %>% magrittr::add(rnorm(length(.), 0, 1)),
    mu = theta_newpar$mu %>% magrittr::add(rnorm(length(.), 0, 1)),
    kp = theta_newpar$kp %>% magrittr::add(rnorm(length(.), 0, 1)),
    nu = theta_newpar$nu
  )
  t_old <- theta %>% unlist()
  t <- theta_newpar %>% unlist()
  t0 <- theta_newpar_0 %>% unlist()
  Qfunc_old <- MinimalEstep(y, t_old)
  Qfunc_new <- EstepReparam(y, t, cache_by_default = T)

  Qfunc_new$set_w(t)
  plot(Qfunc_new)
  Qfunc_new$set_w(t0)
  plot(Qfunc_new)

  Qfunc_new$set_w(t0)
  trace <- GD(Qfunc_new, stop_crit = 500, gamma0 = 1, debug = F)
  plot(trace)
  plot(Qfunc_new)

  Qfunc_new$set_w(t0)
  lr <- polynomial_schedule(1, 0.01, 100)
  trace_sgd <- SGD(Qfunc_new, Adam_Optimizer(lr, batch_size = 100), stop_crit = 100)
  plot(trace_sgd)
  plot(Qfunc_new)

  Qfunc_new$set_w(t0)
  profvis::profvis(
    GD(Qfunc_new, stop_crit = 1000, gamma0 = 1)
  )

  Qfunc_new$set_w(t0)
  profvis::profvis(
    SGD(Qfunc_new, Adam_Optimizer(lr, batch_size = 100), stop_crit = 1000)
  )

})
