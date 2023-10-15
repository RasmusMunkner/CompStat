
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
#'
#' u <- runif(10000)
#' profvis::profvis(iter_bw_est(u, maxiter = 10, kernel="e", cv_k = 3))
#' profvis::profvis(iter_bw_est(u, maxiter = 10, kernel="e", cv_k = 0))
#'
#' res <- iter_bw_est(runif(100), maxiter = 30, kernel = "g", cv_k = 5)
#' ggplot2::ggplot(res, ggplot2::aes(x = 1:nrow(res))) + ggplot2::geom_line(ggplot2::aes(y = bw))
#' ggplot2::ggplot(res, ggplot2::aes(x = 1:nrow(res))) + ggplot2::geom_line(ggplot2::aes(y = cvscore))
#' iter_bw_est(4*rbinom(1000, 1, 0.5) + rnorm(1000))
#' set.seed(NULL)
iter_bw_est <- function(x,
                        maxiter = 3L,
                        kernel = "e",
                        bw0 = NULL,
                        cv_k = 3,
                        evaluation_method = "c",
                        l2_method = "default",
                        tol = 1e-4,
                        reltol = 1e-2,
                        cvtol = -Inf
                        ){
  # Initialize containers ------------------------------------------------------
  kernel <- CompStatKernel(kernel)
  bw_seq <- rep(NA, maxiter + 1)
  fnorm_seq <- rep(NA, maxiter+1)
  cv_seq <- rep(NA, maxiter + 1)

  # Initial choice of bandwidth-------------------------------------------------
  if (is.null(bw0)){
    if (sd(x) == 0){
      stop("Standard deviation of x is 0.
           Kernel density estimation is not meaningful.")
    }
    bw_seq[1] <- 0.9 * sd(x) * length(x)^(1/5) # Silverman
  } else {
    bw_seq[1] <- bw0 # Pre-specified guess
  }

  # Check if cross-validation should be performed-------------------------------
  cv_enabled <- FALSE
  if (cv_k == as.integer(cv_k) & cv_k > 1){
    cv_enabled <- TRUE
    cv_seq[1] <- cvscore(kernel, x, bw_seq[1], cv_k, method=evaluation_method)
  } else if (cv_k != as.integer(cv_k)){
    message("Number of CV splis (cv_k) is not an integer.
            CV will be skipped. Set cv_k = 0 to supress this message.")
  }

  # Pre-computed quantities-----------------------------------------------------
  sigmaK4 <- kernel$sigma2^2
  n <- length(x)
  exit_code <- "Maximum number iterations reached"

  # Iterate --------------------------------------------------------------------
  for (i in 1:maxiter){

    # Update bandwidth
    fnorm_seq[i] <- l2norm(kernel, x, bw_seq[i], method = l2_method)
    bw_seq[i+1] <- (kernel$l2norm / fnorm_seq[i] / sigmaK4 / n)^(1/5)

    # Cross-validation
    if (cv_k > 0){
      cv_seq[i+1] <- cvscore(
        kernel, x, bw_seq[i+1], cv_k,
        method=evaluation_method
        )
    }

    # Bandwidth stopping criteria
    if (
      abs(bw_seq[i+1] - bw_seq[i]) < tol ||
      abs(bw_seq[i+1]-bw_seq[i]) / bw_seq[i] < reltol
        ){
      exit_code <- "Bandwidth tol reached"
      break
    }

    # CV stopping criteria
    if (cv_enabled){
      if (
        !is.finite(cv_seq[i+1]) ||
        cv_seq[i+1] - cv_seq[i] < cvtol
        ){
        if (!is.finite(cv_seq[i+1])){
          exit_code <- "Infinite CV score"
        } else {
          exit_code <- "CV tol reached"
        }
        break
      }
    }
  }

  # Return results -------------------------------------------------------------
  estimates <- tibble::tibble(
    bw=bw_seq[!is.na(bw_seq)],
    fnorm=fnorm_seq[!is.na(bw_seq)]
  )
  if (cv_k > 0){
    estimates <- estimates %>%
      dplyr::mutate(cvscore=cv_seq[!is.na(bw_seq)])
  }

  structure(
    list(
      estimates = estimates,
      x = x,
      kernel = kernel,
      exit_code = exit_code
    ),
    class = "IterBwEstimate"
  )

}

#' Plotting for IterBwEstimate-object produced by iter_bw_est
#'
#' @param IterBwEst A data frame of results produced by iter_bw_est
#'
#' @return A summary plot for the IterBwEstimate-object
#' @export
#'
#' @examples
#' set.seed(0)
#' res_e <- iter_bw_est(runif(1000), maxiter = 30, kernel = "e", cv_k = 3)
#' res_g <- iter_bw_est(runif(1000), maxiter = 30, kernel = "g", cv_k = 3)
#' plot(res_e)
#' plot(res_g)
#' set.seed(NULL)
plot.IterBwEstimate <- function(IterBwEst){
  p1 <- IterBwEst$estimate %>%
    ggplot2::ggplot(ggplot2::aes(x = 1:nrow(IterBwEst$estimate), y = bw)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iteration", y = "Bandwidth")

  if (!is.null(IterBwEst$estimate[["cvscore"]])){
    plot_rows <- 2
    p2 <- IterBwEst$estimate %>%
      ggplot2::ggplot(ggplot2::aes(x = 1:nrow(IterBwEst$estimate), y = cvscore)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(x = "Iteration", y = "CV Score")

    best_cv_iter <- which.min(IterBwEst$estimate$cvscore)

    p3 <- IterBwEst$estimate %>%
      ggplot2::ggplot(ggplot2::aes(x = bw, y = cvscore)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = IterBwEst$estimate$bw[best_cv_iter], color ="red")
      ggplot2::labs(x = "Bandwidth", y = "CV Score")
  } else {
    plot_rows <- 1
    best_cv_iter <- nrow(IterBwEst$estimate)
    p2 <- NULL
    p3 <- NULL
  }

  xdens <- kernel_density(IterBwEst$kernel, IterBwEst$x, IterBwEst$estimate$bw[best_cv_iter], return_grid = T)
  ydens <- kernel_density(IterBwEst$kernel, IterBwEst$x, IterBwEst$estimate$bw[best_cv_iter])
  p4 <- ggplot2::ggplot(mapping = ggplot2::aes(x = xdens, y = ydens)) +
    ggplot2::geom_line(color = "red") +
    ggplot2::labs(x = "x", y = "Estimated Density")

  plots <- list(p1, p2, p3, p4)
  plots <- plots[!(plots %>% purrr::map_lgl(is.null))]
  gridExtra::grid.arrange(
    grobs = plots,
    nrow = plot_rows, ncol = 2,
    top = grid::textGrob(IterBwEst$kernel$name)
  )

}

#' Print method for IterBwEstimate
#'
#' @param x A result from iter_bw_est
#'
#' @export
print.IterBwEstimate <- function(x){
  print(paste0("Iterative Bandwidth Estimate - ", x$kernel))
  print(paste0("Iteration terminated with exit code - ", x$exit_code))
  print(x$estimates)
}
