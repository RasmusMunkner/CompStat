data("poisson")
xz_poisson_prior <- sum(poisson$x * poisson$z)
x_poisson_mean <- mean(poisson$x)

vandermonde <- function(y, n){
  outer(y, 0:n, `^`) %>% t()
}

polycoef <- function(n){
  0:n %>%
    purrr::map_dbl(.f=function(n){
      sum((poisson$x - mean(poisson$x))^n)
    }) / factorial(0:n)
}

#' Approximation of the poisson prior density
#'
#' @param y The points to evaluate the density in.
#' @param maxpow The order of approximation. Should be at least 2.
#' The higher the order, the more accurate (and slower) the approximation.
#'
#' @return The approximate density evaluated at y.
#' @export
#'
#' @examples
#' poisson_prior_log_f_approx(seq(0.001, 1, 0.001))
poisson_prior_log_f_approx <- function(y, maxpow = 4){
  term1 <- y*xz_poisson_prior
  exp_approx <- colSums(vandermonde(y, maxpow) * polycoef(maxpow))
  term2 <- exp(y*x_poisson_mean) * exp_approx
  term1 - term2
}

#' Log-density for the poisson prior (Example Distribution)
#'
#' @param y The points at which to evaluate the log-density.
#'
#' @return The log-density within the given points.
#' @export
#'
#' @examples
#' poisson_prior_log_f(1)
#' seq(1e-4, 1, 1e-4) %>% plot(., poisson_prior_log_f(.), type="l")
poisson_prior_log_f <- function(y){
  term1 <- y*xz_poisson_prior
  term2 <- rowSums(exp(outer(y, poisson$x)))
  result <- term1 - term2
  result[y <= 0] <- -Inf # The above calculations do extrapolate, but the density should not
  return(result)
}

poisson_prior_log_f_prime <- function(y){
  term1 <- xz_poisson_prior
  term2 <- rowSums(matrix(rep(poisson$x, length(y)), nrow=length(y), byrow = T) * exp(outer(y, poisson$x)))
  term1 - term2
}

#' An object mean to represent a distribution.
#'
#' A RandomVariable object should expose a density and a log-density.
#'
#' @param f A density function.
#' @param log_f The logarithm of the density function.
#' @param f_prime The derivative of the density.
#' @param log_f_prime The derivative of the log-density.
#'
#' Note that the log-versions of the density and the differentiated density can
#' can be inferred from one another. This is done only if NULL is passed.
#' Any function that could not be inferred will be replaced with an exception raiser.
#'
#' @return A RandomVariable object.
#' @export
#'
#' @examples
#' X <- RandomVariable(dnorm, function(z) log(dnorm(z)), function(z) function(z){dnorm(z) * (-z)})
#' P1 <- RandomVariable(log_f = log_rejec_dens5,log_f_prime = log_rejec_dens_prime)
#' P_fail <- RandomVariable(f = dnorm)
RandomVariable <- function(f = NULL, log_f = NULL, f_prime = NULL, log_f_prime = NULL){

  if (length(c(f, log_f, f_prime, log_f_prime)) == 0){
    print(c(f, log_f, f_prime, log_f_prime))
    stop("Cannot define a RandomVariable without passing any information.")
  }

  # Handling f
  if (is.null(f)){
    if (is.null(log_f)){
      warning("A density could not be inferred. Use this RV with caution!")
      f <- function(z) stop("The random variable was created without specifying a density.")
    } else {
      f <- function(z) exp(log_f(z))
    }
  }

  # Handling log_f
  if (is.null(log_f)){
    if (is.null(f)){
      warning("A log-density could not be inferred. Use this RV with caution!")
      log_f <- function(z) stop("The random variable was created without specifying a log-density")
    } else {
      log_f <- function(z) log(f(z))
    }
  }

  # Handling f_prime
  if (is.null(f_prime)){
    if (is.null(log_f_prime)){
      warning("A derivative of the density could not be inferred. Use this RV with caution!")
      f_prime <- function(z) stop("The random variable was created without specifying derivative.")
    } else {
      f_prime <- function(z) exp(log_f_prime(z))
    }
  }

  # Handling log_f_prime
  if (is.null(log_f_prime)){
    if (is.null(f_prime) || is.null(f)){
      warning("A density of the log-density could not be inferred. Use this RV with caution!")
      log_f_prime <- function(z) stop("The random variable was created without specifying a differentiated log-density")
    } else {
      log_f_prime <- function(z) f_prime(z) / f(z)
    }
  }

  # Put the thing together
  return (
    structure(list(f = f, log_f = log_f, f_prime = f_prime, log_f_prime = log_f_prime), class = "RandomVariable")
  )
}

#' Convinience function for encoding common distributions into RV's.
#'
#' @param key The name (first letter) of the desired distribution.
#'
#' @return A RandomVariable object corresponding to the distribution specified.
#' @export
#'
#' @examples
#' get_rv("n")
#' get_rv("g")
#' get_rv("p")
get_rv <- function(key = "n"){
  switch(key,
         "n"=RandomVariable(f = dnorm,
                            log_f = function(z) -log(sqrt(2*pi)) - z^2 / 2,
                            f_prime = function(z) dnorm(z) * (-z)),
         "g"=RandomVariable(f = function(z) dgamma(z, 2),
                            log_f = function(z) log(z) - z,
                            f_prime = function(z) -dgamma(z,2) + exp(-z)
                            ),
         "p"=RandomVariable(log_f = poisson_prior_log_f,
                            log_f_prime = poisson_prior_log_f_prime),
         "a"=RandomVariable(log_f = poisson_prior_log_f_approx,
                            log_f_prime = poisson_prior_log_f_prime
                            ),
         "e"=RandomVariable(f = function(z) 3/4 * (1-z^2) * (abs(z) < 1),
                            log_f_prime = function(z) - 2*z / (1 - z^2) * (abs(z) < 1)
                            )
         )
}

#' Density plot for random variables
#'
#' @param rv A RandomVariable
#' @param x A grid of points for which to plot the density
#'
#' @return A plot of the density of the rv
#' @export
#'
#' @examples
plot.RandomVariable <- function(rv, x = seq(-3, 3, 0.01)){
  ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = rv$f(x))) +
    ggplot2::geom_line()
}




