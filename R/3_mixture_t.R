
#' Simulate from a mixture t-distribution
#'
#' @param n The number of simulations
#' @param p A vector of probabilities. Note that it must sum to 1.
#' @param mu A vector of location parameters for the t-components.
#' @param sigma A vector of scale parameters for the t-components.
#' @param nu A vector of shape parameters for the t-components.
#'
#' @return
#' @export
#'
#' @examples
#' rtmix(1e6, c(0.1, 0.4, 0.5), c(10,20,30), c(1,2,3), c(999, 999, 999)) %>% hist()
rtmix <- function(n, p, mu, sigma2, nu){

  if(!all.equal(length(p), length(mu), length(sigma), length(nu))){
    stop("Inputs are not of equal length.")
  }

  Z <- rmultinom(1, size = n, prob = p) %>% as.vector()
  indZ <- Z %>%
    purrr::imap(.f = function(k, n){rep(n,k)}) %>%
    purrr::reduce(.f = c)

  rt(n, df = nu[indZ]) * sqrt(sigma2)[indZ] + mu[indZ]

}

#' Implementation of the E-step for the EM algorithm for mixture-t distributions
#'
#' @param y Observations from the mixture-t distribution of interest
#' @param n_components The number of mixture components
#'
#' @details
#' The object returned is a list of three functions.
#' - Q is the Q-function.
#' - dQ is the gradient of the Q-function.
#' - plain_Q is the Q-function, but modified such that it takes numeric vectors
#' as input. It also implements some safeguards against invalid probabilities.
#'
#'
#' @return An object of class 'CompStatQfunc'.
#' @export
#'
#' @examples
#' Estep <- Estep_Factory_tmix(c(1,2,3))
Estep_Factory_tmix <- function(y, n_components = 4){

  # Conditional marginal density Y
  f_y_by_z <- function(w, index = 1:length(y)){
    1:length(w$p) %>%
      purrr::map(.f = function(k){
        dt((y[index] - w$mu[k])/sqrt(w$sigma2[k]),df=w$nu[k])/sqrt(w$sigma2[k])
      }) %>%
      purrr::reduce(.f = cbind) %>%
      magrittr::set_colnames(NULL)
  }

  # Unconditional marginal density Y
  f_y <- function(w, index = 1:length(y)){
    (f_y_by_z(w, index) %*% w$p) %>% as.vector()
  }

  # Conditional densities for Z
  cond_p <- function(w, index = 1:length(y)){
    outer(1 / f_y(w, index), w$p) * f_y_by_z(w, index)
  }

  dmu <- function(v, w, index = 1:length(y)){
    mudiff <- outer(y[index], 1:length(v$p), function(y, k){
      (v$nu[k]+1)*(y - v$mu[k]) / (v$nu[k] * v$sigma2[k] + (y - v$mu[k])^2)
    })
    colSums(mudiff * cond_p(w, index))
  }

  dsigma2 <- function(v, w, index = 1:length(y)){
    sigmadiff <- outer(y[index], 1:length(v$p), function(y, k){
      (1 - (v$nu[k]+1)*(y - v$mu[k])^2 / (v$nu[k] * v$sigma2[k] + (y - v$mu[k])^2)) / (-2 * v$sigma2[k])
    })
    colSums(sigmadiff * cond_p(w, index))
  }

  # Be aware that the last entry of the resulting vector does not have meaning
  # since there are one fewer probability parameters
  dp <- function(v, w, index = 1:length(y)){
    pSums <- colSums(cond_p(w, index)) / v$p
    pSums <- pSums - pSums[length(v$p)]
    pSums[length(v$p)] <- - sum(pSums) # Ensure gradient does not break p
    pSums
  }

  dQ <- function(v, w, index = 1:length(y)){
    data.frame(
      p = -dp(v, w, index),
      mu = -dmu(v, w, index),
      sigma2 = -dsigma2(v, w, index),
      nu = -rep(0, length(w$nu)) # This is wrong and numDeriv is wrong! I know
      )
  }

  Q <- function(v, w, index = 1:length(y)){
    if(!all.equal(v$p %>% sum(), 1)){
      stop("Probabilities dont add up")
    }
    if(!all.equal(w$p %>% sum(), 1)){
      stop("Probabilities dont add up")
    }
    llz <- outer(y[index], 1:length(v$p), function(y, k){
      -log(v$sigma2[k]) / 2 - (v$nu[k]+1)/2 * log(1 + (y - v$mu[k])^2 / (v$nu[k] * v$sigma2[k])) + log(v$p[k])
    })
    -sum(llz * cond_p(w, index))
  }

  wrap_input <- function(par, input_type = "vw"){
    K <- switch (input_type,
      "vw" = floor(length(par)/8),
      "v" = floor(length(par)/4),
      stop("Unrecognized input type in Estep/wrap_input.")
    )
    if (FALSE){
      stop(paste0(
        "Number of parameters is not a multiple of 8 (or 4 if input type is v).
        Floor(K) = ",
        floor(K))
        )
    }
    v <- data.frame(
      p = c(par[1:(K-1)], 1-sum(par[1:(K-1)])),
      mu = par[(K+1):(2*K)],
      sigma2 = par[(2*K+1):(3*K)],
      nu = par[(3*K+1):(4*K)]
    )
    if (input_type == "vw"){
      w <- data.frame(
        p = c(par[(4*K+1):(5*K - 1)], 1 - sum(par[(4*K+1):(5*K - 1)])),
        mu = par[(5*K+1):(6*K)],
        sigma2 = par[(6*K+1):(7*K)],
        nu = par[(7*K+1):(8*K)]
      )
    } else {
      w <- NULL
    }
    if (input_type == "v" & is.null(memory_w)){ # Ensures theta_prime ok
      memory_w <<- v
    }
    list(v = v, w = w)
  }

  plain_Q <- function(par, index = 1:length(y)){
    inputs <- wrap_input(par)
    Q(inputs$v, inputs$w, index)
  }

  plain_dQ <- function(par, index = 1:length(y)){
    inputs <- wrap_input(par)
    dQ(inputs$v, inputs$w, index)
  }

  # Background weights (used with SGD interface)
  memory_w <- NULL

  # Objective with respect to (theta / v)
  objective <- function(par_v, index = 1:length(y)){
    input <- wrap_input(par_v, input_type = "v")
    Q(input$v, memory_w, index)
  }

  # Gradient with respect to (theta / v)
  grad <- function(par_v, index = 1:length(y)){
    input <- wrap_input(par_v, input_type = "v")
    dQ(input$v, memory_w, index) %>% unlist()
  }

  # Updates the background distributional parameters (theta_prime / w)
  update_w <- function(par_v){
    input <- wrap_input(par_v, input_type = "v")
    memory_w <<- input$v
  }

  # Checks if a set of parameters is valid
  check_par_validity <- function(par_v){
    input <- wrap_input(par_v, input_type = "v")
    out <- list(p = F, mu = T, sigma2 = F, nu = F)
    if (
      isTRUE(all.equal(sum(input$v$p), 1)) &&
      all(0 < input$v$p) &&
      all(input$v$p < 1)
      ){
      out$p <- T
    }
    if (all(input$v$sigma2 > 0)){
      out$sigma2 <- T
    }
    if (all(input$v$sigma2 > 0)){
      out$sigma2 <- T
    }
    if (all(input$v$nu > 0)){
      out$nu <- T
    }
    out
  }

  # Density plot
  plotdens <- function(par_v, xx = seq(-10, 10, 0.01)){
    input <- wrap_input(par_v, input_type = "v")
    yy <- xx %>% purrr::map_dbl(.f = function(x){
      1:n_components %>% purrr::map_dbl(.f = function(k){
        dt((x-input$v$mu[k]) / sqrt(input$v$sigma2[k]), df = input$v$nu[k]) /
          sqrt(input$v$sigma2[k])
      }) %>%
        magrittr::multiply_by(input$v$p) %>%
        sum()
    })
    ggplot2::ggplot(mapping = ggplot2::aes(x = xx, y = yy)) +
      ggplot2::geom_line()
  }

  structure(list(
    Q = Q,
    dQ = dQ,
    plain_Q = plain_Q,
    plain_dQ = plain_dQ,
    objective = objective,
    grad = grad,
    update_w = update_w,
    plotdens = plotdens,
    check_par_validity = check_par_validity,
    par_alloc = c(
      rep("p", n_components),
      rep("mu", n_components),
      rep("sigma2", n_components),
      rep("nu", n_components)
      ),
    n_param = 4 * n_components,
    n_index = length(y)
  ), class = c("CompStatQfunc", "CompStatOptimizable"))
}















