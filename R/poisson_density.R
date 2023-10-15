data("poisson")
xz_poisson_prior <- sum(poisson$x * poisson$z)
x_poisson_mean <- mean(poisson$x)

#' Build the Vandermonde matrix for y
#'
#' @param y A numeric vector.
#' @param n The order of the Vandermonde matrix.
#'
#' @return A matrix containing $y^k$ for $k = 0, \ldots, n$.
#' @export
#'
#' @examples
#' vandermonde(c(1,2,3), 4)
vandermonde <- function(y, n){
  outer(y, 0:n, `^`) %>% t()
}

#' Calculates coefficients for the Taylor polynomial arising from the poisson
#' prior.
#'
#' @param n The order of the Taylor expansion
#' @param x The data used for calculating the coefficients
#'
#' @return A vector of polynomial coefficients for the Taylor polynomial
#' @export
#'
#' @examples
#' polycoef(3)
#' polycoef(3, poisson$x[1:50])
polycoef <- function(n, x = poisson$x){
  x_centralized <- x - mean(x)
  0:n %>%
    purrr::map_dbl(.f=function(n){
      sum(x_centralized^n)
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
poisson_prior_log_f_approx <- function(y, breaks = c(0,1,2,3,4), orders = 4){
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
