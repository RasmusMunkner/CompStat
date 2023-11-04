



test_that("rc++ implementation is consistent", {

  set.seed(0)
  theta <- list(
    p = c(0.2, 0.4, 0.1, 0.2),
    mu = c(-3, 1, 8, 4, 5),
    sigma2 = c(1, 8, 2, 3, 9),
    nu = c(2, 3, 6, 1, 4)
  )
  y <- rtmix(1e4, theta$p, theta$mu, theta$sigma2, theta$nu)
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

  A <- Qfunc_new$grad(t)
  B <- Qfunc_new$grad(t, mode = "c")
  C <- Qfunc_new$grad(t, mode = "c2")

  expect_equal(A, B)
  expect_equal(A, C)

  profvis::profvis(
    for (i in 1:100){
      Qfunc_new$grad(t, cache = F, mode ="r")
    },
    simplify = F
  )

  profvis::profvis(
    for (i in 1:100){
      Qfunc_new$grad(t, cache = F, mode ="c")
    },
    simplify = F
  )

  for (i in 1:10){Qfunc_new$grad(t, cache = F, mode ="r")}

  # Shows that there is no magic - Its just caching
  microbenchmark::microbenchmark(
    cond_pC(t_old, y),
    dQC2(t_old, t_old, y),
    environment(Qfunc_new$set_w)$cond_pR(t_old, 1:length(y), y),
    environment(Qfunc_new$set_w)$cond_p(t_old, 1:length(y), y),
    environment(Qfunc_new$set_w)$cond_p_wrapper(t_old, 1:length(y), y, cache = F),
    environment(Qfunc_new$set_w)$cond_p_wrapper(t_old, 1:length(y), y, cache = T)
  )

  clear_cache <- function(mode = "r"){
    environment(Qfunc_new$grad)$cond_p_cache <<- NULL
    Qfunc_new$grad(t, cache = T, mode = mode)
  }

  microbenchmark::microbenchmark(
    clear_cache(mode = "r"),
    clear_cache(mode = "c"),
    Qfunc_new$grad(t, mode = "c"),
    Qfunc_new$grad(t, mode = "c2")
  )


})
