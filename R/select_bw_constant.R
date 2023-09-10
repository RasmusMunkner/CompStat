#' Bandwidth selection for kernel smoothing
#'
#' Select an optimal bandwidth for kernel smoothing based on the AMISE.
#' This relies on an initial guess f_0 for the unknown density.
#'
#' @param f_0 An intial guess at the unknown density.
#' @param K The kernel used for smoothing.
#'
#' @return The estimated optimal bandwidth (a positive number).
#'
#' @examples
#' select_bw_constant(dnorm, dnorm)
#' select_bw_constant(dnorm, function(x) ifelse(abs(x) < 1, 3/4 * (1-x^2), 0))
select_bw_constant <- function(
    f_0 = dnorm,
    K = dnorm
){
  K_norm <- l2norm(K)
  f_0_norm <- l2norm(second_derivative(f_0))
  sigma_K <- moment(K, n=2)
  (K_norm / f_0_norm / sigma_K^2)^(1/5)
}

#' Numerical second derivative
#'
#' Calculates the second derivative of the input function using numDeriv.
#'
#' @param f The function that is be differentiated.
#'
#' @return The second derivative of f.
#'
#' @examples
#' second_derivative(log)
#' second_derivative(exp)
#' second_derivative(function(x) x^3 / 6)
second_derivative <- function(f){
  Vectorize(function(z) numDeriv::hessian(f, z)[1,1])
}

#' Numerical L2-Norm calculation
#'
#' @inheritParams stats::integrate
#' @param ... additional arguments passed to integrate or f
#'
#' @return The L2-Norm of f.
#'
#' @examples
#' l2norm(dnorm)
#' l2norm(dexp)
#' l2norm(function(x) 1/sqrt(x), lower = 1, upper = exp(1))
l2norm <- function(f, lower = -Inf, upper = Inf, ...){
  integrate(function(z) f(z)^2, lower=lower, upper=upper, ...)$value
}

#' Moments given a density via numerical integration
#'
#' @param k A univariate density function.
#' @param n A vector of exponents.
#'
#' @return A vector of moments corresponding to the input.
#'
#' @examples
#' moment(dnorm, c(1,2,3,4))
#' moment(dexp, c(1,2))
#' moment(dunif, c(1,2))
moment <- function(k, n = 2){
  n %>%
    purrr::map_dbl(.f=function(n){
      integrate(function(z) z^n * k(z), lower=-Inf, upper=Inf)$value
    })
}


