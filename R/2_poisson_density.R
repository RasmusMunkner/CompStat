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
vandermonde <- function(y, n, mode = "r"){
  switch(mode,
         "c" = vmC(y,n),
         "c2" = vmC2(y,n),
         "r" = outer(y, 0:n, `^`) %>% t(),
         "poly" = cbind(rep(1,length(y)), poly(y, degree = n, raw = T, simple = T)) %>% t()
         )
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

#' Approximations of the Poisson density
#'
#' @param x X-values for the poisson prior
#' @param z Z-values for the poisson prior
#' @param breaks Breaks used to split X's
#' @param orders The order of Taylor approximation to use. Can be a vector,
#' but should probably be a scalar, since we use broadcasting.
#'
#' @return An approximating function for the poisson prior.
#' @export
#'
#' @examples
#' p1 <- poisson_prior_approximation(K = 8, breaks = c(0,1,2,3,4), vm_method="c2")
#' ppas <- purrr::map(.x = 1:64, .f = function(t) {poisson_prior_approximation(K=t)})
#' y <- purrr::map_dbl(.x = ppas, .f = function(ppa) ppa(1)-poisson_prior_log_f(1))
#'
#' plot(y, z1-z2, type = "l")
#'
poisson_prior_approximation <- function(
    x = poisson$x, z = poisson$z, breaks = c(0,1,2,3,4), K = 4,
    vm_method = "r"){
  xz <- sum(x*z)
  groups <- cut(x, breaks, labels = F)
  x_means <- 1:max(groups) %>%
    purrr::map_dbl(.f = function(k){
      mean(x[groups == k])
    })
  coef <- list(k = 1:max(groups), n = K) %>%
    purrr::pmap(.f = function(k, n){
      polycoef(n, x[groups == k])
    }) %>%
    purrr::reduce(.f = rbind)
  approximation <- function(y){
    VMy <- vandermonde(y, K, vm_method)
    p <- coef %*% VMy
    expmean <- exp(outer(x_means, y))
    y*xz - colSums(p * expmean)
  }
  approximation
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
poisson_prior_log_f_approx <- function(y, breaks = c(0,1,2,3,4), K = 4){
  term1 <- y*xz_poisson_prior
  exp_approx <- colSums(vandermonde(y, K) * polycoef(K))
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

#' Derivative for the log of the poisson prior
#'
#' @param y The points for the derivative to be evaluated in
#'
#' @return The derivative evaluated at the specified points
#' @export
poisson_prior_log_f_prime <- function(y){
  term1 <- xz_poisson_prior
  term2 <- rowSums(matrix(rep(poisson$x, length(y)), nrow=length(y), byrow = T) * exp(outer(y, poisson$x)))
  term1 - term2
}

#' Benchmark for the approximate implementations of the poisson prior
#'
#' @param N_seq A sequence of input sizes (logscale) for which the implementations should be implemented.
#' @param rg The range on which the implementations are tested.
#'
#' @return A data frame containing benchmark results
#' @export
#'
#' @examples
#' bm <- benchmark_poisson_prior_densities()
#' bm %>%
#' ggplot2::ggplot(ggplot2::aes(x = log(N), y = log(Time), color = Method)) +
#' ggplot2::geom_line() +
#' ggplot2::geom_point()
benchmark_poisson_prior_densities <- function(N_seq = 3:12, rg = c(0.0001, 1)){

  ppa_nogroup <- poisson_prior_approximation(breaks = c(0,4), K = 8)
  ppa_base <- poisson_prior_approximation(breaks = c(0, 1, 2, 3, 4), K = 8)
  ppa_loworder <- poisson_prior_approximation(breaks = c(0, 1, 2, 3, 4), K = 4)
  ppa_inaccurate <- poisson_prior_approximation(breaks = c(0, 4), K = 4)

  bm <- N_seq %>%
    purrr::imap_dfr(.f = function(N, i){
      calls <- list(
        call("ppa_nogroup", y = seq(rg[1], rg[2], length.out = 2^N)),
        call("ppa_base", y = seq(rg[1], rg[2], length.out = 2^N)),
        call("ppa_loworder", y = seq(rg[1], rg[2], length.out = 2^N)),
        call("ppa_inaccurate", y = seq(rg[1], rg[2], length.out = 2^N)),
        call("poisson_prior_log_f", y = seq(rg[1], rg[2], length.out = 2^N))
      )
      names(calls) <- c("No Grouping, n = 8", "Grouping, n = 8", "Grouping, n = 4", "No Grouping, n = 4", "Exact")
      microbenchmark::microbenchmark(
        list=calls
      ) %>%
        dplyr::group_by(expr) %>%
        dplyr::summarise(Time = median(time)) %>%
        dplyr::mutate(N = 2^N, Method = expr) %>%
        dplyr::select(-c(expr))
    })
  bm
}
