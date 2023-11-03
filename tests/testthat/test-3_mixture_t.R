


test_that("Derivation of Q-function is correct", {

  set.seed(0)
  y <- rt(30, df = 3)
  N <- 4

  theta <- data.frame(
    p = rnorm(N) %>% exp(.) %>% magrittr::divide_by(.,sum(.)),
    mu = rnorm(N),
    sigma2 = rexp(N, 1),
    nu = 2 + rexp(N, 1)
    )

  theta_prime <- data.frame(
    #p = theta$p,
    p = rnorm(N) %>% exp(.) %>% magrittr::divide_by(.,sum(.)),
    #mu = theta$mu,
    mu = rnorm(N),
    #sigma2 = theta$sigma2,
    sigma2 = rexp(N, 1),
    #nu = theta$nu
    nu = 2 + rexp(N, 1)
  )

  Estep <- Estep_Factory_tmix(y)

  par <- c(unlist(theta), unlist(theta_prime))
  A <- numDeriv::grad(Estep$plain_Q, par) %>%
    `[`(1:(4*N)) %>%
    matrix(ncol = 4)

  B <- Estep$dQ(theta, theta_prime)

  checks <- unlist(A[,1:(ncol(A)-1)] - B[,1:(ncol(A)-1)], use.names = F)

  expect_equal(checks, rep(0, length(checks)))

  set.seed(NULL)

})
