
#' L2 Norm through matrix expansion
#'
#' @param x A numeric vector
#' @param r A positive bandwidth
#'
#' @return The estimated L2-norm of the kernel density estimate for the data
#' using the epanechnikov kernel.
#' @export
#'
#' @examples
#' epanechnikov_l2norm_matrix(rnorm(100), 1)
epanechnikov_l2norm_matrix <- function(x, r){
  x <- 9/4 * x / r^6 / length(x)^2
  r <- 9/4 * r / r^6 / length(x)^2
  rep_x <- matrix(rep(x, length(x)), nrow=length(x))
  diff_matrix <- 2*r - abs(rep_x - t(rep_x))
  sum(diff_matrix[diff_matrix > 0])
}

#' L2 Norm for compact data
#'
#' This method only works if the range of the data is shorter than 2 times the
#' bandwidth. An error is thrown if this in not the case.
#'
#' @inheritParams epanechnikov_l2norm_matrix
#' @param do_sort Whether to sort the input. Defaults to true. Unsorted input
#' will produce wrong results.
#'
#' @return The estimated L2-norm of the kernel density estimate for the data
#' using the epanechnikov kernel.
#' @export
#'
#' @examples
#' epanechnikov_l2norm_naive(runif(100), 0.51)
epanechnikov_l2norm_naive <- function(x, r, do_sort = TRUE){
  if (do_sort){
    x <- sort(x) # If x is sorted we can skip this one for speed
  }
  n <- length(x)
  if (x[n] - x[1] > 2*r){
    stop("Naive method only works for data with |x_n-x_1| <= 2r")
  }
  return(
    9/4 * (2*r*n^2 + 2*sum(x*(n-2*seq_along(x)+1))) / r^6 / n^2
  )
}

#' L2 Norm based on running windows
#'
#' Calculates the L2-norm for the kernel density estimate through cleverly
#' exploiting the sorted data to avoid summing over 0's.
#'
#' @inheritParams epanechnikov_l2norm_matrix
#'
#' @return The estimated L2-norm of the kernel density estimate for the data
#' using the epanechnikov kernel.
#' @export
#'
#' @examples
#' epanechnikov_l2norm_running(rnorm(1000), 1)
epanechnikov_l2norm_running <- function(x, r){
  x <- sort(x)
  n <- length(x)
  x[n+1] <- Inf
  results <- 0
  i <- 2
  j <- 1
  while(j < n){
    if(i <= j || x[i+1] - x[j] < 2*r){
      i <- i + 1
    } else if (x[i] - x[j] < 2*r) {
      results <- results + sum(2*r - x[(j+1):i] + x[j])
      j <- j + 1
    } else {
      j <- j+1
    }
  }
  results <- 2*results + 2*r*n
  return(
    results * 9 / 4 / n^2 / r^6
  )
}

epanechnikov_l2norm_running2 <- function(x, r, agg = TRUE, debug = FALSE){
  if (debug){
    browser()
  }
  x <- sort(x)
  n <- length(x)
  x[n+1] <- Inf
  results <- 0#matrix(rep(0, n^2), nrow=n)
  i <- 1
  j <- 1
  k <- 1
  PUSH <- TRUE
  while(j <= n){
    if (PUSH){
      if (x[i+1] - x[j] < 2*r) { # Can we push further? - Push further
        i <- i + 1
      } else if (x[i] - x[j] < 2*r) { # Do we stop pushing?
        for (s in k:i){
          #results[j:s, s] <- 2*r - (x[s] - x[j:s])
          results <- results + sum(2*r - (x[s] - x[j:s]))
        }
        k <- i+1
        j <- j+1
        PUSH <- FALSE
      }
    } else {
      if(x[i+1] - x[j] < 2*r){
        PUSH <- TRUE
      } else {
        j <- j+1
      }
    }
  }
  if (agg){
    return(
      (2 * sum(results) - 2*r*n) * 9/4 / n^2 / r^6
    )
  } else {
    return(results)
  }

}

#' L2 Norm based on exact binning
#'
#' @inheritParams epanechnikov_l2norm_matrix
#' @param helper Another epanechnikov-l2norm function to be used within each bin
#'
#' @return The estimated L2-norm of the kernel density estimate for the data
#' using the epanechnikov kernel.
#' @export
#'
#' @examples
#' epanechnikov_l2norm_binning(rnorm(1000), 1)
epanechnikov_l2norm_binning <- function(
    x, r,
    helper = epanechnikov_l2norm_matrix
    ){
  x <- sort(x)
  n <- length(x)
  if (x[n] - x[1] <= 2*r){
    return(epanechnikov_l2norm_naive(x, r))
  } else {
    bin_index <- floor((x - x[1])/(2*r))
    x_binned <- split(x, bin_index)
    partial_results <- vector(mode="numeric", length = length(x_binned) - 1)
    for (i in 1:(length(x_binned) - 1)){
      partial_results[i] <- helper(
        x=c(x_binned[[i]], x_binned[[i+1]]),
        r=r
      ) * (length(x_binned[[i]]) + length(x_binned[[i+1]]))^2
      if (i > 1){
        partial_results[i] <- partial_results[i] - epanechnikov_l2norm_naive(
          x=x_binned[[i]],
          r=r,
          do_sort=FALSE
        ) * (length(x_binned[[i]]))^2
      }
    }
    return(sum(partial_results) / length(x)^2)
  }
}
