

# Note that we are not testing for derivatives wrt. p_K
test_that("Derivation of Q-function is correct", {

  set.seed(0)
  y <- rt(30, df = 3)
  N <- 4

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

  Estep <- Estep_Factory_tmix(y, init_par = theta_prime)

  Estep$objective(theta)
  Estep$grad(theta)

  A <- numDeriv::grad(Estep$objective, unlist(theta))

  B <- Estep$grad(theta)

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

  Estep <- Estep_Factory_tmix(y, init_par = theta_prime, n_components = N)

  Estep$set_w(theta)

  A <- Estep$fisher(unlist(theta), method = 1)

  B <- Estep$fisher(unlist(theta), method = 2)

  C <- Estep$fisher(unlist(theta), method = 3)

  D <- Estep$fisher(unlist(theta), method = 4)

  expect_equal(A, B)
  expect_equal(max(abs(A - D)) < 1, TRUE) # If this fails, set the tolerance up

  set.seed(NULL)

})











