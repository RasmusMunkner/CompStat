

# Note that we are not testing for derivatives wrt. p_K
test_that("Derivation of Q-function is correct", {

  set.seed(0)
  y <- rt(1, df = 3)
  K <- 2

  theta <- list(
    p = rnorm(K) %>% exp(.) %>% magrittr::divide_by(.,sum(.)) %>% `[`(-1),
    mu = rnorm(K),
    sigma2 = rexp(K, 1),
    nu = 2 + rexp(K, 1)
    ) %>% unlist()

  theta_prime <- list(
    p = rnorm(K) %>% exp(.) %>% magrittr::divide_by(.,sum(.)) %>% `[`(-1),
    mu = rnorm(K),
    sigma2 = rexp(K, 1),
    nu = 2 + rexp(K, 1)
  ) %>% unlist()

  Estep <- MinimalEstep(y, init_par = theta_prime)

  #Estep_Old <- Estep_Factory_tmix(y, init_par = theta_prime, K = K)
  #Estep$get_w() - unlist(Estep_Old$get_w())[-K]
  #Estep$objective(theta) - Estep_Old$objective(theta)

  #Estep$grad(theta) - Estep_Old$grad(theta)

  A <- numDeriv::grad(Estep$objective, theta)[1:(3*K-1)]

  B <- Estep$grad(theta)[1:(3*K-1)]

  expect_equal((A-B) %>% magrittr::set_names(NULL), rep(0, length(A)))

  set.seed(NULL)

})

test_that("fisher information is calculated correctly", {

  set.seed(0)
  N <- 2
  theta <- list(
    p = rnorm(N) %>% exp(.) %>% magrittr::divide_by(.,sum(.)) %>% `[`(-1),
    mu = rnorm(N),
    sigma2 = rexp(N, 1),
    nu = 2 + rexp(N, 1)
  )
  theta_prime <- list(
    p = rnorm(N) %>% exp(.) %>% magrittr::divide_by(.,sum(.)) %>% `[`(-1),
    mu = rnorm(N),
    sigma2 = rexp(N, 1),
    nu = 2 + rexp(N, 1)
  )
  n <- 1e3
  y <- rtmix(n, theta$p, theta$mu, theta$sigma2, theta$nu)

  theta <- theta %>% unlist()
  theta_prime <- theta_prime %>% unlist()

  Estep <- Estep_Factory_tmix(y, init_par = theta_prime)

  Estep$set_w(theta)

  A <- Estep$fisher(unlist(theta), method = 1)

  B <- Estep$fisher(unlist(theta), method = 2)

  C <- Estep$fisher(unlist(theta), method = 3)

  D <- Estep$fisher(unlist(theta), method = 4)

  expect_equal(A, B)
  expect_equal(max(abs(A - D)) < 1, TRUE) # If this fails, set the tolerance up

  set.seed(NULL)

})











