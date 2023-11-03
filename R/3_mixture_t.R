
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
#' rtmix(1e4, c(0.1, 0.4, 0.5), c(10,20,30), c(1,2,3), c(999, 999, 999)) %>% hist()
#' rtmix(1e4, c(0.1, 0.4), c(10,20,30), c(1,2,3), c(999, 999, 999)) %>% hist()
#' rtmix(1e4, c(0.1, 0.4), c(10,20,30), c(1,2,3), c(999, 999)) # Produces error
rtmix <- function(n, p, mu, sigma2, nu){

  if (length(p) == length(mu) - 1){
    p <- c(p, 1-sum(p))
  }

  lengths <- c(length(p), length(mu), length(sigma2), length(nu))
  if(max(lengths) > min(lengths)){
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
Estep_Factory_tmix <- function(y, init_par, n_components = 4){

  # Converts parameters to wide-format
  par_wrap <- function(par){
    par <- par %>% unlist(use.names = F)
    K <- n_components
    data.frame(
      p = c(par[1:(K-1)], 1-sum(par[1:(K-1)])),
      mu = par[(K):(2*K-1)],
      sigma2 = par[(2*K):(3*K-1)],
      nu = par[(3*K):(4*K-1)]
    )
  }
  memory_w <- init_par %>% par_wrap()
  par_wrap_safe <- function(par){
    par <- par %>% unlist(use.names = F)
    v <- par_wrap(par)
    v$nu <- memory_w$nu
    v
  }
  # Updates the background distributional parameters (theta_prime / w)
  set_w <- function(par_v){
    memory_w <<- par_wrap_safe(par_v)
  }
  get_w <- function(){
    memory_w
  }
  # Checks if a set of parameters is valid
  check_par_validity <- function(par_v){
    v <- par_wrap(par_v)
    out <- list(p = F, mu = T, sigma2 = F, nu = F)
    if (
      all(0 < v$p)
    ){
      out$p <- T
    }
    if (all(v$sigma2 > 0)){
      out$sigma2 <- T
    }
    if (all(v$sigma2 > 0)){
      out$sigma2 <- T
    }
    if (all(v$nu > 0)){
      out$nu <- T
    }
    out
  }

  # Conditional marginal density Y
  f_y_by_z <- function(w, index = 1:length(y), y0 = y){
    1:length(w$p) %>%
      purrr::map(.f = function(k){
        dt((y0[index] - w$mu[k])/sqrt(w$sigma2[k]),df=w$nu[k])/sqrt(w$sigma2[k])
      }) %>%
      purrr::reduce(.f = cbind) %>%
      magrittr::set_colnames(NULL)
  }

  # Unconditional marginal density Y
  f_y <- function(w, index = 1:length(y), y0 = y){
    (f_y_by_z(w, index, y0) %*% w$p) %>% as.vector()
  }

  # Conditional densities for Z
  cond_p <- function(w, index = 1:length(y)){
    outer(1 / f_y(w, index), w$p) * f_y_by_z(w, index)
  }

  dmu <- function(v, w, index = 1:length(y), sum = T){
    mudiff <- outer(y[index], 1:length(v$p), function(y, k){
      (v$nu[k]+1)*(y - v$mu[k]) / (v$nu[k] * v$sigma2[k] + (y - v$mu[k])^2)
    })
    if (sum){
      colSums(mudiff * cond_p(w, index))
    } else {
      mudiff * cond_p(w, index)
    }
  }

  dsigma2 <- function(v, w, index = 1:length(y), sum = T){
    sigmadiff <- outer(y[index], 1:length(v$p), function(y, k){
      (1 - (v$nu[k]+1)*(y - v$mu[k])^2 / (v$nu[k] * v$sigma2[k] + (y - v$mu[k])^2)) / (-2 * v$sigma2[k])
    })
    if (sum){
      colSums(sigmadiff * cond_p(w, index))
    } else {
      sigmadiff * cond_p(w, index)
    }
  }

  # Be aware that the last entry of the resulting vector does not have meaning
  # since there are one fewer probability parameters
  dp <- function(v, w, index = 1:length(y), sum = T){
    if (sum){
      pSums <- colSums(cond_p(w, index)) / v$p
      pSums <- pSums - pSums[length(v$p)]
      return(pSums[-n_components])
    } else {
      pfrac <- cond_p(w, index) / matrix(rep(v$p, length(index)), nrow = length(index), byrow = T)
      pfrac <- pfrac - matrix(rep(pfrac[,length(v$p)],length(v$p)), ncol = length(v$p))
      return(pfrac[,-n_components])
    }
  }

  dQ <- function(v, w, index = 1:length(y)){
    list(
      p = -dp(v, w, index),
      mu = -dmu(v, w, index),
      sigma2 = -dsigma2(v, w, index),
      nu = -rep(0, length(w$nu)) # This is wrong and numDeriv is wrong! I know
      ) %>% unlist()
  }

  Q <- function(v, w, index = 1:length(y)){
    llz <- outer(y[index], 1:length(v$p), function(y, k){
      -log(v$sigma2[k]) / 2 - (v$nu[k]+1)/2 * log(1 + (y - v$mu[k])^2 / (v$nu[k] * v$sigma2[k])) + log(v$p[k])
    })
    -sum(llz * cond_p(w, index))
  }

  # Objective with respect to (theta / v)
  objective <- function(par_v, index = 1:length(y)){
    v <- par_wrap_safe(par_v)
    Q(v, memory_w, index)
  }

  # Gradient with respect to (theta / v)
  grad <- function(par_v, index = 1:length(y)){
    v <- par_wrap_safe(par_v)
    dQ(v, memory_w, index)
  }

  ## Likelihood
  loglikelihood <- function(par_v, index = 1:length(y), y0 = y){
    v <- par_wrap_safe(par_v)
    f_y(v, index, y0) %>% log() %>% sum()
  }

  # Fisher information
  fisher <- function(par_v_mle, index = 1:length(y), method = 1){

    v_mle <- par_wrap_safe(par_v_mle)

    # Note that there is a minus built into dQ called by plain_dQ
    SymmGrad <- function(par_v){
      v <- par_wrap_safe(par_v)
      dQ(v, v, index)
    }
    Grad2MLE <- function(par_v){
      v <- par_wrap_safe(par_v)
      dQ(v, v_mle, index)
    }
    AsymGradMLE <- function(par_v){
      v <- par_wrap_safe(par_v)
      dQ(v_mle, v, index)
    }

    if (method == 1){
      return(numDeriv::jacobian(SymmGrad, par_v_mle))
    }

    if (method == 2){
      return(numDeriv::jacobian(Grad2MLE, par_v_mle) +
               numDeriv::jacobian(AsymGradMLE, par_v_mle))
    }

    if (method == 3){
      gradmat <- list(
        dp(v_mle, v_mle, index, sum = F),
        dmu(v_mle, v_mle, index, sum = F),
        dsigma2(v_mle, v_mle, index, sum = F),
        matrix(rep(0, length(index) * n_components), ncol = n_components)
      ) %>%
        purrr::reduce(.f = cbind)
      loglik <- colMeans(gradmat)
      gradmat <- gradmat - matrix(rep(loglik, nrow(gradmat)), ncol = ncol(gradmat), byrow = T)
      return(
        1:nrow(gradmat) %>%
          purrr::map(.f = function(i){
            outer(gradmat[i,], gradmat[i,])
          }) %>%
          purrr::reduce(.f = `+`)
      )
    }

    if (method == 4){
      return(-optimHess(par_v_mle, loglikelihood))
    }

    stop("Invalid method specified.")

  }

  structure(list(
    objective = objective,
    grad = grad,
    set_w = set_w,
    get_w = get_w,
    loglikelihood = loglikelihood,
    fisher = fisher,
    check_par_validity = check_par_validity,
    par_alloc = c(
      rep("p", n_components-1),
      rep("mu", n_components),
      rep("sigma2", n_components),
      rep("nu", n_components)
      ),
    n_param = 4 * n_components - 1,
    n_index = length(y)
  ), class = c("CompStatQfunc", "CompStatOptimizable"))
}

#' Plot function for CompStatQfunc
#'
#' @param Qfunc A CompStatQfunc. Needs to have memory_w set.
#'
#' @return
#' @export
plot.CompStatQfunc <- function(Qfunc){
  env <- environment(Qfunc$set_w)
  w <- env$memory_w
  if (is.null(w)){
    stop("Qfunc does not have a set memory_w. Please call update_w() first.")
  }
  rg <- range(w$mu) + c(-4,4) * max(sqrt(w$sigma2))
  xx <- seq(rg[1], rg[2], length.out = 1e3)
  yy <- env$f_y(w, index = 1:length(xx), y0 = xx)
  yyz <- env$f_y_by_z(w, index = 1:length(xx), y0 = xx) *
    matrix(rep(w$p, length(xx)), ncol = length(w$p), byrow = T)
  yyz <- yyz %>%
    as.data.frame() %>%
    magrittr::set_colnames(
      paste0("t(mu, sigma2, nu)~(",
             round(w$mu,1), ", ", round(w$sigma2,1), ", ", round(w$nu,1),
             ")")) %>%
    dplyr::mutate(Combined = rowSums(.)) %>%
    dplyr::mutate(xx = xx) %>%
    tidyr::pivot_longer(cols = -xx, names_to = "Component", values_to = "value")

  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = yyz$xx, y = yyz$value, color = yyz$Component)) +
    ggplot2::guides(color=ggplot2::guide_legend(title = "Component")) +
    ggplot2::labs(x = "y", y = "f(y)")
}















