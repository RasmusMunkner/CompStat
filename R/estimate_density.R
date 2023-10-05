
#' Iterative bandwidth estimation
#'
#' Iterative bandwidth selection using AMISE.
#' The method runs until convergence of bandwidth selection and returns all
#' bandwidth estimates. To access density estimates, use compile_density.
#'
#' @param x A numerical vector. Data to fit the density to
#' @param maxiter A positive integer. The maximal number of iterations
#' @param bw0 A scalar. Initial guess at the optimal bandwidth.
#' Defaults to maximum-likelihood Gaussian density
#' @param kernel_code One of "gaussian" or "epanechnikov".
#' @param cv_k Number of folds to use for cross validation.
#' @param tol A scalar. Tolerance level bandwidth selection convergence
#' @param reltol A scalar. Relative tolerance for bandwidth selection convergence
#' @param cvtol A scalar. Absolute tolerance for cv-score improvement
#' @param ... Additional arguments passed to estimate_l2norm.
#'
#' @return A list of bandwidths for estimation of the unknown density.
#' The final bandwidth represents the best-estimate for the optimal bandwidth.
#' @export
#'
#' @examples
#' set.seed(0)
#' iter_bw_est(runif(100), maxiter = 10, kernel="e", cv_k = 3)
#' res <- iter_bw_est(runif(100), maxiter = 30, kernel = "g", cv_k = 5) %>% as.data.frame()
#' ggplot(res, aes(x = 1:nrow(res))) + geom_line(aes(y = bw))
#' ggplot(res, aes(x = 1:nrow(res))) + geom_line(aes(y = cvscore))
#' iter_bw_est(4*rbinom(1000, 1, 0.5) + rnorm(1000))
#' set.seed(NULL)
iter_bw_est <- function(x,
                        maxiter = 3L,
                        kernel = "e",
                        bw0 = NULL,
                        cv_k = 5,
                        tol = 1e-7,
                        reltol = 1e-3,
                        cvtol = -Inf,
                        ...
                        ){
  #Initialize containers -------------------------------------------------------

  if (!("CompStatKernel" %in% class(kernel))){
    tryCatch(
      error = function(cnd){
        stop("Kernel must be specified as either a identifying string or an object of class 'CompStatKernel'")
      },
      kernel <- get_kernel(kernel)
    )
  }

  bw_seq <- rep(NA, maxiter + 1)
  if (is.null(bw0)){
    if (sd(x) == 0){
      stop("Standard deviation of x is 0. Automatic bandwidth selection is not meaningful.")
    }
    bw_seq[1] <- 0.9 * sd(x) * length(x)^(1/5) # Silverman
  } else {
    bw_seq[1] <- bw0 # Pre-specified guess
  }

  fnorm_seq <- rep(NA, maxiter+1)

  cv_seq <- rep(NA, maxiter + 1)
  cv_seq[1] <- cvscore(kernel, x, bw_seq[1], cv_k)

  # Pre-computed quantities
  sigmaK4 <- kernel$sigma2^2
  n <- length(x)

  # Iterate --------------------------------------------------------------------
  for (i in 1:maxiter){

    fnorm_seq[i] <- l2norm(kernel, x, bw_seq[i], ...)

    bw_seq[i+1] <- (kernel$l2norm / fnorm_seq[i] / sigmaK4 / n)^(1/5)
    cv_seq[i+1] <- cvscore(kernel, x, bw_seq[i+1], cv_k)

    if (
      abs(bw_seq[i+1] - bw_seq[i]) < tol ||
      abs(bw_seq[i+1]-bw_seq[i]) / bw_seq[i] < reltol ||
      !is.finite(cv_seq[i+1]) || # Can happen that the epanechnikov kernel misses a point during cv. Then cv score = -Inf
      cv_seq[i+1] - cv_seq[i] < cvtol
        ){
      break
    }

  }

  # Return results -------------------------------------------------------------
  tibble::tibble(
    bw=bw_seq[!is.na(bw_seq)],
    fnorm=fnorm_seq[!is.na(bw_seq)],
    cvscore=cv_seq[!is.na(bw_seq)]
  )

}


#' Plot density estimates
#'
#' @param x Data for which density is computed.
#' @param bw The bandwidth. Can be a vector of bandwidths.
#' @param kernel_code The identifier for which kernel is to be used.
#' @param grid_length The number of evaluation points.
#'
#' @return A ggplot with the density estimates for each bandwidth.
#' @export
#'
#' @examples
#' plot_density(rnorm(1000), 0.1)
plot_density <- function(x, bw, kernel_code = "g", grid_length = 200){
  z <- seq(range(x)[1] - 3*max(bw), range(x)[2] + 3*max(bw), length.out=grid_length)
  bw %>%
    purrr::map_dfr(.f=function(h){
      tibble::tibble(
        bw = h,
        y = compile_density(z, x, h, kernel_code),
        x = z
      ) %>%
        dplyr::mutate(bw = factor(bw))
    }) %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = y, color = bw)) +
      ggplot2::geom_line()
}
