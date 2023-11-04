


test_that("GD works", {


  set.seed(0)
  theta <- list(
    p = c(0.2, 0.4),
    mu = c(-3, 1, 8),
    sigma2 = c(1, 8, 2),
    nu = c(2, 3, 6)
  )
  y <- rtmix(1e4, theta$p, theta$mu, theta$sigma2, theta$nu)
  theta <- theta %>% unlist()

  theta_0 <- list(
    p = c(0.33, 0.33),
    mu = rnorm(3, 0, 9),
    sigma2 = rep(var(y) / 100, 3),
    nu = c(2,3,6)
  ) %>% unlist()
  Qfunc <- MinimalEstep(y, theta_0)

  Qfunc$set_w(0.3*theta + 0.7*theta_0)
  plot(Qfunc)

  trace <- GD(Qfunc, d = 0.7, stop_crit = 50)

  plot(Qfunc)

  plot(trace)

  Qfunc$loglikelihood(theta_0)
  Qfunc$loglikelihood(trace %>% tail())
  Qfunc$loglikelihood(theta)


})
