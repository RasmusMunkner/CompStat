
#' Iterative density estimation
#'
#' Kernel-based density estimation with in-built bandwidth selection using AMISE.
#' The method runs until convergence of bandwidth selection and returns all
#' bandwidth estimates. To access density estimates, use compile_density.
#'
#' @param x A numerical vector. Data to fit the density to
#' @param maxiter A positive integer. The maximal number of iterations
#' @param bw0 A scalar. Initial guess at the optimal bandwidth.
#' Defaults to maximum-likelihood Gaussian density
#' @param kernel_code One of "gaussian", "uniform", "triangular" or "epanechnikov".
#' @param tol A scalar. Tolerance level bandwidth selection convergence
#' @param reltol A scalar. Relative tolerance for bandwidth selection convergence
#'
#' @return A list of bandwidths for estimation of the unknown density.
#' The final bandwidth represents the best-estimate for the optimal bandwidth.
#' @export
#'
#' @examples
#' iter_dens_est(rnorm(1000))
#' iter_dens_est(4*rbinom(1000, 1, 0.5) + rnorm(1000))
iter_dens_est <- function(x,
                          maxiter = 3L,
                          kernel_code = "g",
                          bw0 = NULL,
                          tol = 1e-7,
                          reltol = 1e-3
                          ){

  bw_seq <- rep(NA, maxiter + 1)
  if (is.null(bw0)){
    bw_seq[1] <- 0.9 * sd(x) * length(x)^(1/5) # Silverman
  } else {
    bw_seq[1] <- bw0 # Prespecified guess
  }

  H <- get_kernel(kernel_code)
  for (i in 2:(maxiter+1)){
    l2est <- estimate_l2norm(x, kernel_code, bw_seq[i-1])
    bw_seq[i] <- (H$l2norm / l2est / H$sigma2^2 / length(x))^(1/5)
    if (abs(bw_seq[i] - bw_seq[i-1]) < tol ||
        abs(bw_seq[i]-bw_seq[i-1]) / bw_seq[i-1] < reltol){
      break
    }
  }

  bw_seq[!is.na(bw_seq)]

}

#' Compiles kernel density estimate
#'
#' @param x Data used for density estimation
#' @param bw Bandwidth
#' @param kernel Kernel used for smoothing. Defaults to Gaussian
#'
#' @return The estimated density function
#' @export
#'
#' @examples
#' set.seed(0)
#' x <- rnorm(1000)
#' z <- seq(-3,3,0.01)
#' dens_hat <- compile_density(z, x, 0.1)
#' ggplot() +
#'   geom_histogram(aes(x = x, y = after_stat(density))) +
#'   geom_line(aes(x = z, y = dens_hat), color = "red")
#' set.seed(NULL)
#'
compile_density <- function(z, x, bw, kernel_code = "g"){
  H <- get_kernel(kernel_code)$kernel
  z %>%
    purrr::map_dbl(.f=function(z){
      mean(H((z-x)/bw)/bw)
    })
}

#TODO - Make a nice plotting facility for this (vectorize bw)
plot_density <- function(x, bw, kernel_code = "g"){
  z <- seq(range(x)[1] - 3*bw, range(x)[2] + 3*bw, length.out=200)
  dens_hat <- compile_density(z, x, bw, kernel_code)
  ggplot() +
    geom_histogram(aes(x=x, y=after_stat(density))) +
    geom_line(aes(x=z, y = dens_hat), color="red")
}
