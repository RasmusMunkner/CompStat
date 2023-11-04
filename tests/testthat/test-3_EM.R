



test_that("rc++ implementation is consistent", {


  set.seed(0)
  t <- list(
    be = c(0.2, 0.4, 0.1, 0.2) %>% p_() %>% inverse_softmax() %>% `[`(-length(.)),
    mu = c(-3, 1, 8, 2, -7),
    kp = log(c(1, 8, 2, 2, 4)),
    nu = c(2, 3, 6, 2, 8)
  ) %>% unlist()
  Qfunc <- EstepReparam(0, t)
  y <- rtmix(1e4, Qfunc$get$p(t), Qfunc$get$mu(t),
             Qfunc$get$s2(t), Qfunc$get$nu(t))
  Qfunc <- EstepReparam(y, t) # Wraps the previously shown implementations

  A <- Qfunc$grad(t)
  B <- Qfunc$grad(t, mode = "c")
  C <- Qfunc$grad(t, mode = "c2")

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

  for (i in 1:10){Qfunc$grad(t, cache = F, mode ="r")}

  # Shows that there is no magic - Its just caching
  tc <- c(Qfunc$get$p(t), Qfunc$get$mu(t),Qfunc$get$s2(t), Qfunc$get$nu(t))
  microbenchmark::microbenchmark(
    cond_pC(tc, y),
    environment(Qfunc$set_w)$cond_pR(tc, 1:length(y), y),
    environment(Qfunc$set_w)$cond_p(tc, 1:length(y), y),
    environment(Qfunc$set_w)$cond_p_wrapper(tc, 1:length(y), y, cache = F),
    environment(Qfunc$set_w)$cond_p_wrapper(tc, 1:length(y), y, cache = T)
  )

  clear_cache <- function(mode = "r"){
    environment(Qfunc$grad)$cond_p_cache <<- NULL
    Qfunc$grad(t, cache = T, mode = mode)
  }

  microbenchmark::microbenchmark(
    clear_cache(mode = "r"),
    clear_cache(mode = "c"),
    Qfunc$grad(t, mode = "c"),
    Qfunc$grad(t, mode = "c2")
  )


})
