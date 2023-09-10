
#' Iterative density estimation
#'
#' Kernel-based density estimation with in-built bandwidth selection using AMISE.
#' The method runs until convergence of bandwidth selection and returns all
#' intermediate density and bandwidth estimates.
#'
#' @param x A numerical vector. Data to fit the density to
#' @param maxiter A positive integer. The maximal number of iterations
#' @param f_0 A scalar function. Initial guess at the unknown density.
#' Defaults to maximum-likelihood Gaussian density
#' @param kernel A scalar function. The kernel used for smoothing
#' @param tol A scalar. Tolerance level bandwidth selection convergence
#' @param reltol A scalar. Relative tolerance for bandwidth selection convergence
#'
#' @return A list of bandwidths/estimates for the unknown density.
#' @export
#'
#' @examples
#' iter_dens_est(rnorm(1000))
#' iter_dens_est(4*rbinom(1000, 1, 0.5) + rnorm(1000))
iter_dens_est <- function(x,
                          maxiter = 3L,
                          f_0 = NULL,
                          kernel = get_kernel("gaussian"),
                          tol = 1e-7,
                          reltol = 1e-3
                          ){

  # A reasonable initial guess at a density is Gaussian with
  # the same mean and std.dev. as the observed data
  if (is.null(f_0)){
    f_0 <- function(z) dnorm(z, mean=mean(x), sd = sd(x))
  }

  bw <- vector("numeric", maxiter+1) #Bandwidths
  bw[1] <- NA
  dens_est <- vector("list", maxiter+1) #Density estimation functions
  dens_est[[1]] <- f_0

  for (i in 2:(maxiter+1)){

    bw[i] <- select_bw_constant(dens_est[[i-1]], kernel) * length(x)^(-1/5)

    # Check if convergence is reached
    if (!is.na(bw[i-1])){
      if (abs(bw[i] - bw[i-1]) < tol || abs((bw[i] - bw[i-1])/bw[i-1]) < reltol){
        return(
          list(bw = bw[2:(i-1)], dens_est = dens_est[2:(i-1)])
        )
      }
    }

    dens_est[[i]] <- compile_density(x, bw[i], kernel)
  }

  return(
    list(bw = bw[2:(maxiter+1)], dens_est = dens_est[2:(maxiter+1)])
  )

}

#' Compiles kernel density estimate
#'
#' @param x Data used for density estimation
#' @param bw Bandwidth
#' @param kernel Kernel used for smoothing. Defaults to Gaussian
#'
#' @return The estimated density function
#'
#' @examples
#' set.seed(0)
#' x <- rnorm(1000)
#' dens_hat <- compile_density(x, 0.1)
#' ggplot() +
#'   geom_histogram(aes(x = x, y = after_stat(density))) +
#'   geom_line(aes(x = x, y = dens_hat(x)), color = "red")
#' set.seed(NULL)
#'
compile_density <- function(x, bw, kernel = get_kernel("g")){

  # The arguments are evaluated non-lazily (see 10.2.3 AdvR)
  force(x)
  force(bw)
  force(kernel)

  # Compile the density function
  g <- function(z){
    z %>%
      purrr::map_dbl(
        .f=function(z){
          mean(kernel((z-x)/bw)/bw)
        }
      )
  }
  return(g)
}

#' Access commonly used kernels
#'
#' The valid codes are gaussian/normal, rectangular/uniform, epanechnikov and triangular.
#' These can also be accessed by the corresponding 1-letter codes.
#'
#' @param kernel_code An alias for the desired kernel.
#' @param hessian If true, returns the hessian of the kernel instead.
#'
#' @return The kernel corresponding to the given input.
#' @export
#'
#' @examples
#' get_kernel("gaussian")
#' get_kernel("g")
#' get_kernel("t")
get_kernel <- function(kernel_code, hessian=FALSE){
  if (!hessian){
    switch(kernel_code %>% substr(1,1) %>% tolower(),
           n= ,
           g= dnorm,
           u= ,
           r= function(z){ifelse(abs(z)<1, 1/2, 0)},
           t= function(z){ifelse(abs(z)<1, 1-abs(z), 0)},
           e= function(z){ifelse(abs(z)<1, 3/4 * (1-z^2), 0)}
    )
  } else {
    switch(kernel_code %>% substr(1,1) %>% tolower(),
           n= ,
           g= function(z){z^2 * dnorm(z) - dnorm(z)},
           u= ,
           r= function(z){rep(0, length(x))},
           t= function(z){rep(0, length(x))},
           e= function(z){ifelse(abs(z) < 1, -3/2, 0)}
    )
  }
}

#' Calculate the ISE
#'
#' Calculates Integrated-Square-Error between two densities.
#'
#' @param d1 A density function
#' @param d2 A density function
#' @inheritParams l2norm
#'
#' @return A scalar
#' @export
#'
#' @examples
#' ise(dnorm, dexp)
#' ise(dnorm, function(z) dnorm(z, 1))
ise <- function(d1, d2, lower = -Inf, upper = Inf){
  l2norm(function(z) d1(z) - d2(z), lower = lower, upper = upper)
}

# library(ggplot2)
# n <- 1000
# x <- rnorm(n)
# p <- rbinom(n, prob = 0.5, size = 1)
# x <- p * rnorm(n) + (1-p) * rnorm(n, mean=4)
# grid <- seq(-3,7,0.01)
#
# true_dens <- function(x){
#   dnorm(x) * 0.5 + dnorm(x, 4) * 0.5
# }
#
# density_est <- iter_dens_est(x, 20)
# ggplot() +
#   geom_histogram(aes(x = x, y = after_stat(density))) +
#   geom_density(aes(x = x)) +
#   geom_line(aes(x = grid, y = true_dens(grid)), color = "orange") +
#   geom_line(aes(x = grid, y = density_est$dens_est[[1]](grid)), color = "red") +
#   geom_line(aes(x = grid, y = density_est$dens_est[[3]](grid)), color = "blue") +
#   geom_line(aes(x = grid, y = density_est$dens_est[[10]](grid)), color = "green") +
#   geom_line(aes(x = grid, y = density_est$dens_est[[10]](grid)), color = "yellow")
#
#
# for (i in seq_along(density_est$dens_est)){
#   print(paste0("Estimate ", i, ": ", sum((true_dens(grid) - density_est$dens_est[[i]](grid))^2)))
# }
#
# ise <- seq_along(density_est$dens_est) %>%
#   purrr::map_dbl(.f=function(i){
#     sum((true_dens(grid) - density_est$dens_est[[i]](grid))^2)
#   })

