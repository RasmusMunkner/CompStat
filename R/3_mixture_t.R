
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
Estep_Factory_tmix <- function(y, init_par, K = ceiling(length(init_par) / 4)){

  # Converts parameters to wide-format
  par_wrap <- function(par){
    par <- par %>% unlist(use.names = F)
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
    memory_w %>% unlist() %>% `[`(-K)
  }
  # Checks if a set of parameters is valid
  check_par_validity <- function(v){
    #v <- par_wrap(par_v)
    check_p <- function(x){
      0 < x & x < 1 & sum(x) < 1
    }
    check_pos <- function(x){
      0 < x
    }
    c(
      check_p(v[1:(K-1)]),
      rep(TRUE, K),
      check_pos(v[(2*K):(3*K-1)]),
      check_pos(v[(3*K):(4*K-1)])
    )
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
      return(pSums[-K])
    } else {
      pfrac <- cond_p(w, index) / matrix(rep(v$p, length(index)), nrow = length(index), byrow = T)
      pfrac <- pfrac - matrix(rep(pfrac[,length(v$p)],length(v$p)), ncol = length(v$p))
      return(pfrac[,-K])
    }
  }

  dQ <- function(v, w, index = 1:length(y)){
    -c(
      dp(v, w, index),
      dmu(v, w, index),
      dsigma2(v, w, index),
      rep(0, length(w$nu))
    )
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
    check_par = check_par_validity,
    par_alloc = c(
      rep("p", K-1),
      rep("mu", K),
      rep("sigma2", K),
      rep("nu", K)
      ),
    n_param = 4 * K - 1,
    n_index = length(y)
  ), class = c("CompStatQfunc", "CompStatOptimizable"))
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
MinimalEstep <- function(y, init_par, K = ceiling(length(init_par) / 4)){

  # Parameters are parsed as a NumericVector
  # There are (K - 1) probability parameters
  # There are K mean parameters
  # There are K scale parameters
  # There are K shape parameters

  # Updates the background distributional parameters (theta_prime / w)
  memory_w <- init_par
  set_w <- function(v){
    memory_w <<- v
  }
  get_w <- function(){
    memory_w
  }

  # Parameters are accessed through functions
  p <- function(w, k = 1:K){
    w_extend <- c(w[1:(K-1)], 1-sum(w[1:(K-1)]))
    return(w_extend[k])
  }
  mu <- function(w, k = 1:K){
    w[K - 1 + k]
  }
  s2 <- function(w, k = 1:K){
    w[2*K - 1 + k]
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
  cond_p <- function(w, index = 1:length(y), y0 = y){
    (outer(1 / f_y(w, index, y0), p(w)) * f_y_by_z(w, index, y0)) %>%
      magrittr::set_colnames(paste0("p", 1:K))
  }

  cond_p_wrapper <- function(w, index = 1:length(y), y0 = y, mode = "r"){
    if (mode == "c"){
      cond_pC(c(1-sum(p(w)), w), y0[index])
    } else if (mode == "r"){
      cond_p(w, index, y0)
    }
  }

  dmu <- function(v, w, index = 1:length(y), y0 = y, sum = T){
    mudiff <- outer(y0[index], 1:K, function(y, k){
      (nu(v,k)+1)*(y - mu(v,k)) / (nu(v,k) * s2(v,k) + (y - mu(v,k))^2)
    })
    if (sum){
      colSums(mudiff * cond_p(w, index, y0)) %>% magrittr::set_names(paste0("mu", 1:K))
    } else {
      mudiff * cond_p(w, index, y0) %>% magrittr::set_names(paste0("mu", 1:K))
    }
  }

  dsigma2 <- function(v, w, index = 1:length(y), y0 = y, sum = T){
    sigmadiff <- outer(y0[index], 1:K, function(y, k){
      (1 - (nu(v,k)+1)*(y - mu(v,k))^2 / (nu(v,k) * s2(v,k) + (y - mu(v,k))^2)) / (-2 * s2(v,k))
    })
    if (sum){
      colSums(sigmadiff * cond_p(w, index, y0)) %>% magrittr::set_names(paste0("sigma2", 1:K))
    } else {
      sigmadiff * cond_p(w, index, y0) %>% magrittr::set_names(paste0("sigma2", 1:K))
    }
  }

  # Be aware that the last entry of the resulting vector does not have meaning
  # since there are one fewer probability parameters
  dp <- function(v, w, index = 1:length(y), y0 = y, sum = T){
    if (sum){
      pSums <- colSums(cond_p(w, index, y0)) / p(v)
      pSums <- pSums - pSums[K]
      return(pSums[-K])
    } else {
      pfrac <- cond_p(w, index, y0) / matrix(rep(p(v), length(index)), nrow = length(index), byrow = T)
      pfrac <- pfrac - matrix(rep(pfrac[,K],K), ncol = K)
      return(pfrac[,-K])
    }
  }

  dQ <- function(v, w, index = 1:length(y), y0 = y){
    -c(
      dp(v, w, index, y0),
      dmu(v, w, index, y0),
      dsigma2(v, w, index, y0),
      rep(0, K) %>% magrittr::set_names(paste0(paste0("nu", 1:K)))
    )
  }

  Q <- function(v, w, index = 1:length(y), y0 = y){
    llz <- outer(y0[index], 1:K, function(y, k){
      -log(s2(v,k)) / 2 - (nu(v,k)+1)/2 * log(1 + (y - mu(v,k))^2 / (nu(v,k) * s2(v,k))) + log(p(v,k))
    })
    -sum(llz * cond_p(w, index, y0))
  }

  # Objective with respect to (theta / v)
  objective <- function(v, index = 1:length(y), y0 = y){
    if (all(check_par(v))){
      Q(v, memory_w, index, y0)
    } else {
      Inf
    }
  }

  # Gradient with respect to (theta / v)
  grad <- function(v, index = 1:length(y), y0 = y){
    dQ(v, memory_w, index, y0)
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

  # Check and minimally corrects parameters that exceed the allowed boundaries
  validate_par <- function(v, eps_p = 1e-8, eps_s2 = 1e-4, eps_nu = 1){
    vp <- p(v)
    #vmu <- mu(v) This is always fine
    vs2 <- s2(v)
    vnu <- nu(v)

    vp <- ifelse(vp > 1-eps_p, 1 - rep(eps_p,K),
          ifelse(vp < eps_p,   rep(eps_p,K),
                               vp))
    vp[K] <- 1 - vp[-K]

    vs2 <- ifelse(vs2 < eps_s2, rep(eps_s2, K), vs2)
    vnu <- ifelse(vnu < eps_nu, rep(eps_nu, K), vnu)

    c(vp, vmu, vs2, vnu)
  }

  check_par <- function(v){
    check_p <- function(x){
      0 < x & x < 1 & sum(x) < 1
    }
    check_pos <- function(x){
      0 < x
    }
    c(
      check_p(v[1:(K-1)]),
      rep(TRUE, K),
      check_pos(v[(2*K):(3*K-1)]),
      check_pos(v[(3*K):(4*K-1)])
      )
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
    validate_par = validate_par,
    check_par = check_par,
    get = list(p = p, mu=mu, s2 = s2, nu = nu),
    n_param = 4 * K - 1,
    K = K,
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
  w <- Qfunc$get_w()
  rg <- range(Qfunc$get$mu(w)) +
    c(-4,4) * max(sqrt(Qfunc$get$s2(w)))
  xx <- seq(rg[1], rg[2], length.out = 1e3)
  yy <- Qfunc$f_y(w, index = 1:length(xx), y0 = xx)
  yyz <- Qfunc$f_y_by_z(w, index = 1:length(xx), y0 = xx) *
    matrix(rep(Qfunc$get$p(w), length(xx)), ncol = Qfunc$K, byrow = T)
  yyz <- yyz %>%
    magrittr::set_colnames(
      paste0("t", 1:ncol(.), "(p, mu, sigma2, nu)~(",
             round(Qfunc$get$p(w),3), ", ",
             round(Qfunc$get$mu(w),1), ", ",
             round(Qfunc$get$s2(w),1), ", ",
             round(Qfunc$get$nu(w),1),
             ")")) %>%
    as.data.frame() %>%
    dplyr::mutate(Combined = rowSums(.)) %>%
    dplyr::mutate(xx = xx) %>%
    tidyr::pivot_longer(cols = -xx, names_to = "Component", values_to = "value")

  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = yyz$xx, y = yyz$value, color = yyz$Component)) +
    ggplot2::guides(color=ggplot2::guide_legend(title = "Component")) +
    ggplot2::labs(x = "y", y = "f(y)")
}















