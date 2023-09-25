data("poisson")
xz_poisson_prior <- sum(poisson$x * poisson$z)

#' Log-density for the poisson prior (Example Distribution)
#'
#' @param y The points at which to evaluate the log-density.
#'
#' @return The log-density within the given points.
#' @export
#'
#' @examples
#' poisson_prior_log_f(1)
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
                            log_f_prime = poisson_prior_log_f_prime)
         )
}




