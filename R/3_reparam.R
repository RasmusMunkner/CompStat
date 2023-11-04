#' Softmax
#'
#' @param args The numbers to calculate softmax for
#'
#' @return
#' @export
#'
#' @examples
#' softmax(rep(0,3))
#' softmax(c(1,2,3))
#' softmax(c(749, 750, 751))
softmax <- function(args){
  args <- args - mean(args)
  denom <- sum(exp(args))
  exp(args-log(denom))
}

#' Inverse to the softmax under identifiability condition
#'
#' @param p Vector of probabilities summing to 1.
#' @param last_beta The convention used for the final coefficient.
#'
#' @return The vector of coefficients such that softmax(beta) = p.
#' @export
#'
#' @examples
#' beta <- c(rnorm(100),0)
#' p <- softmax(beta)
#' beta_back <- inverse_softmax(p)
#' max(abs(beta-beta_back))
#'
#' p <- c(0.2, 0.4)
#' p %>% p_() %>% inverse_softmax() %>% softmax()
inverse_softmax <- function(p, last_beta = 0){
  logC <- last_beta - log(p[length(p)])
  log(p) + logC
}

#' Extend incomplete vector of probabilities to a complete one.
#'
#' @param p_incomplete
#'
#' @return
#' @export
#'
#' @examples
#' p_(c(0.1, 0.3))
p_ <- function(p_incomplete){
  p_incomplete %>% c(., 1-sum(.))
}

#' Implementation of the E-step for the EM algorithm for mixture-t distributions
#'
#' @param y Observations from the mixture-t distribution of interest
#' @param K The number of mixture components
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
EstepReparam <- function(y, init_par, K = ceiling(length(init_par) / 4), cache_by_default = T){

  # Parameters are parsed as a NumericVector
  # There are (K - 1) probability parameters
  # There are K mean parameters
  # There are K scale parameters
  # There are K shape parameters

  # Parameters should be accessed via accessor functions

  # Updates the background distributional parameters (theta_prime / w)
  memory_w <- init_par
  set_w <- function(v){
    memory_w <<- v
  }
  get_w <- function(){
    memory_w
  }

  # Parameters are accessed through functions
  be <- function(w, k = 1:(K-1)){
    w[k]
  }
  p <- function(w, k = 1:K){
    softmax(c(w[1:(K-1)], 0))[k]
  }
  mu <- function(w, k = 1:K){
    w[K - 1 + k]
  }
  kp <- function(w, k = 1:K){
    w[2*K - 1 + k]
  }
  s2 <- function(w, k = 1:K){
    exp(kp(w, k))
  }
  nu <- function(w, k = 1:K){
    memory_w[3*K - 1 + k]
  }

  # Conditional marginal density Y
  f_y_by_z <- function(w, index = 1:length(y), y0 = y){
    1:K %>%
      purrr::map(.f = function(k){
        dt((y0[index] - mu(w, k))/sqrt(s2(w, k)), df=nu(w,k))/sqrt(s2(w, k))
      }) %>%
      purrr::reduce(.f = cbind) %>%
      magrittr::set_colnames(NULL)
  }

  # Unconditional marginal density Y
  f_y <- function(w, index = 1:length(y), y0 = y){
    (f_y_by_z(w, index, y0) %*% p(w)) %>% as.vector()
  }

  # Conditional densities for Z
  cond_pR <- function(w, index = 1:length(y), y0 = y){
    (outer(1 / f_y(w, index, y0), p(w)) * f_y_by_z(w, index, y0)) %>%
      magrittr::set_colnames(paste0("p", 1:K))
  }

  cond_p <- function(w, index = 1:length(y), y0 = y, mode = "r"){
    if (mode == "c"){
      cond_pC(c(1-sum(p(w)), w), y0[index])
    } else if (mode == "r"){
      cond_pR(w, index, y0)
    }
  }

  cond_p_par_cache_w <- NULL
  cond_p_par_cache_index <- NULL
  cond_p_par_cache_y0 <- NULL
  cond_p_cache <- NULL
  cond_p_wrapper <- function(w, index = 1:length(y), y0 = y, cache = cache_by_default, mode = "r"){
    if (cache){
      if (
        !is.null(cond_p_cache) &
        isTRUE(all.equal(cond_p_par_cache_w, w)) &
        isTRUE(all.equal(cond_p_par_cache_index, index)) &
        isTRUE(all.equal(cond_p_par_cache_y0, y0))
        ){
        return(cond_p_cache)
      } else {
        cond_p_par_cache_w <<- w
        cond_p_par_cache_index <<- index
        cond_p_par_cache_y0 <<- y0
        cond_p_cache <<- cond_p(w, index, y0, mode)
        return(cond_p_cache)
      }
    } else {
      cond_p(w, index, y0, mode)
    }
  }

  dmu <- function(v, w, index = 1:length(y), y0 = y, sum = T, cache = cache_by_default){
    mudiff <- outer(y0[index], 1:K, function(y, k){
      (nu(v,k)+1)*(y - mu(v,k)) / (nu(v,k) * s2(v,k) + (y - mu(v,k))^2)
    })
    if (sum){
      colSums(mudiff * cond_p_wrapper(w, index, y0, cache=cache)) %>% magrittr::set_names(paste0("mu", 1:K))
    } else {
      mudiff * cond_p_wrapper(w, index, y0, cache=cache) %>% magrittr::set_names(paste0("mu", 1:K))
    }
  }

  dkappa <- function(v, w, index = 1:length(y), y0 = y, sum = T, cache = cache_by_default){
    kappadiff <- outer(y0[index], 1:K, function(y, k){
      (1 - (nu(v,k)+1)*(y - mu(v,k))^2 / (nu(v,k) * s2(v,k) + (y - mu(v,k))^2)) / (-2)
    })
    if (sum){
      colSums(kappadiff * cond_p_wrapper(w, index, y0, cache=cache)) %>% magrittr::set_names(paste0("kappa", 1:K))
    } else {
      kappadiff * cond_p_wrapper(w, index, y0, cache=cache) %>% magrittr::set_names(paste0("kappa", 1:K))
    }
  }

  dsigma2 <- function(v, w, index = 1:length(y), y0 = y, sum = T, cache = cache_by_default){
    sigmadiff <- outer(y0[index], 1:K, function(y, k){
      (1 - (nu(v,k)+1)*(y - mu(v,k))^2 / (nu(v,k) * s2(v,k) + (y - mu(v,k))^2)) / (-2 * s2(v,k))
    })
    if (sum){
      colSums(sigmadiff * cond_p_wrapper(w, index, y0, cache=cache)) %>% magrittr::set_names(paste0("sigma2", 1:K))
    } else {
      sigmadiff * cond_p_wrapper(w, index, y0, cache=cache) %>% magrittr::set_names(paste0("sigma2", 1:K))
    }
  }

  dbeta <- function(v, w, index = 1:length(y), y0 = y, sum = T, cache = cache_by_default){
    if (sum){
      (colSums(cond_p_wrapper(w, index, y0, cache=cache)) - length(index) * p(v))[-K]
    } else{
      (cond_p_wrapper(w, index, y0, cache=cache) - matrix(rep(p(v), length(index)), nrow = length(index)))[,-K]
    }
  }

  # Be aware that the last entry of the resulting vector does not have meaning
  # since there are one fewer probability parameters
  dp <- function(v, w, index = 1:length(y), y0 = y, sum = T, cache = cache_by_default){
    if (sum){
      pSums <- colSums(cond_p_wrapper(w, index, y0, cache=cache)) / p(v)
      pSums <- pSums - pSums[K]
      return(pSums[-K])
    } else {
      pfrac <- cond_p_wrapper(w, index, y0, cache=cache) / matrix(rep(p(v), length(index)), nrow = length(index), byrow = T)
      pfrac <- pfrac - matrix(rep(pfrac[,K],K), ncol = K)
      return(pfrac[,-K])
    }
  }

  dQ <- function(v, w, index = 1:length(y), y0 = y, cache = cache_by_default, mode = "r"){
    if (mode == "r"){
      -c(
        dbeta(v, w, index, y0, cache=cache),
        dmu(v, w, index, y0, cache=cache),
        dkappa(v, w, index, y0, cache=cache),
        rep(0, K)
      ) %>%
        magrittr::set_names(names(memory_w))
    } else if (mode == "c"){
      v_cpar <- c(p(v), mu(v), s2(v), nu(v))
      w_cpar <- c(p(w), mu(w), s2(w), nu(w))
      dQC(v_cpar, w_cpar, y0[index]) %>%
        magrittr::multiply_by(-1) %>%
        as.vector() %>% # Output comes as a matrix
        `[`(-K) %>% # CPP version spits out a dB_K, which we don't want
        c(rep(0, K)) %>% # Adds dNu (which we don't consider)
        magrittr::set_names(names(memory_w))
    } else if (mode == "c2"){
      v_cpar <- c(p(v), mu(v), s2(v), nu(v))
      w_cpar <- c(p(w), mu(w), s2(w), nu(w))
      dQC2(v_cpar, w_cpar, y0[index]) %>%
        magrittr::multiply_by(-1) %>%
        as.vector() %>% # Output comes as a matrix
        `[`(-K) %>% # CPP version spits out a dB_K, which we don't want
        c(rep(0, K)) %>% # Adds dNu (which we don't consider)
        magrittr::set_names(names(memory_w))
    }
  }

  Q <- function(v, w, index = 1:length(y), y0 = y, cache = cache_by_default){
    llz <- outer(y0[index], 1:K, function(y, k){
      -log(s2(v,k)) / 2 - (nu(v,k)+1)/2 * log(1 + (y - mu(v,k))^2 / (nu(v,k) * s2(v,k))) + log(p(v,k))
    })
    -sum(llz * cond_p_wrapper(w, index, y0, cache=cache))
  }

  # Objective with respect to (theta / v)
  objective <- function(v, index = 1:length(y), y0 = y, cache = cache_by_default){
    Q(v, memory_w, index, y0, cache=cache)
  }

  # Gradient with respect to (theta / v)
  grad <- function(v, index = 1:length(y), y0 = y, cache = cache_by_default, mode = "r"){
    dQ(v, memory_w, index, y0, cache=cache, mode = mode)
  }

  ## Likelihood
  loglikelihood <- function(v, index = 1:length(y), y0 = y){
    f_y(v, index, y0) %>% log() %>% sum()
  }

  # Fisher information
  fisher <- function(v_mle, index = 1:length(y), method = 1){

    # Note that there is a minus built into dQ called by plain_dQ
    SymmGrad <- function(v){
      dQ(v, v, index)
    }
    Grad2MLE <- function(v){
      dQ(v, v_mle, index)
    }
    AsymGradMLE <- function(v){
      dQ(v_mle, v, index)
    }

    if (method == 1){
      return(numDeriv::jacobian(SymmGrad, v_mle))
    }

    if (method == 2){
      return(numDeriv::jacobian(Grad2MLE, v_mle) +
               numDeriv::jacobian(AsymGradMLE, v_mle))
    }

    if (method == 3){
      gradmat <- list(
        dp(v_mle, v_mle, index, sum = F),
        dmu(v_mle, v_mle, index, sum = F),
        dsigma2(v_mle, v_mle, index, sum = F),
        matrix(rep(0, length(index) * K), ncol = K)
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
      return(-optimHess(v_mle, loglikelihood))
    }

    stop("Invalid method specified.")

  }

  check_par <- function(v){
    rep(TRUE, 4*K-1)
  }

  structure(list(
    objective = objective,
    grad = grad,
    set_w = set_w,
    get_w = get_w,
    loglikelihood = loglikelihood,
    f_y = f_y,
    f_y_by_z = f_y_by_z,
    cond_p = cond_p_wrapper,
    fisher = fisher,
    check_par = check_par,
    get = list(p = p, mu=mu, s2 = s2, nu = nu),
    n_param = 4 * K - 1,
    K = K,
    n_index = length(y)
  ), class = c("CompStatQfunc", "CompStatOptimizable"))
}
