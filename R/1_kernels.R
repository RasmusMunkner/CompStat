
#' Instantiate a CompStatKernel
#'
#' The valid codes are gaussian/normal, rectangular/uniform, epanechnikov and triangular.
#' These can also be accessed by the corresponding 1-letter codes.
#'
#' @param kernel_code An alias for the desired kernel.
#'
#' @return The kernel corresponding to the given input.
#' @export
#'
#' @examples
#' set.seed(0)
#' x <- rnorm(100000)
#' ekern <- CompStatKernel("e")
#' profvis::profvis(kernel_density(ekern, x, 0.1, method = "r"))
CompStatKernel <- function(kernel_code){

  if("CompStatKernel" %in% class(kernel_code)){
    return(kernel_code)
  } else {
    kernel_code <- kernel_code %>% substr(1,1) %>% tolower()
    valid_kernels <- c("n", "g", "e")
    if (!(kernel_code %in% valid_kernels)){
      stop(purrr::reduce(.x = c("Invalid kernel code. Valid codes are:", valid_kernels), .f = paste))
    }
  }

  name <- switch(kernel_code,
    n= ,
    g= "Gaussian",
    u= ,
    r= "Rectangular",
    t= "Triangular",
    e= "Epanechnikov"
  )
  kernel <- switch(kernel_code,
         n= ,
         g= dnorm,
         u= ,
         r= function(z){ifelse(abs(z)<1, 1/2, 0)},
         t= function(z){ifelse(abs(z)<1, 1-abs(z), 0)},
         #e= function(z){ifelse(abs(z)<1, 3/4 * (1-z^2), 0)}, # Old implementation, ~30% slower
         e= epa_kernel
  )
  hessian <- switch(kernel_code,
         n= ,
         g= function(z){(z^2 - 1) * dnorm(z)},
         u= ,
         r= function(z){rep(0, length(z))},
         t= function(z){rep(0, length(z))},
         e= function(z){ifelse(abs(z) < 1, -3/2, 0)}
  )
  l2norm <- switch(kernel_code,
                   n = ,
                   g = 3/8 / sqrt(pi),
                   u = ,
                   r = 1/2,
                   t = 2/3,
                   e = 3/5
                   )
  sigma2 <- switch(kernel_code,
                  n = ,
                  g = 1,
                  u = ,
                  r = 1/3,
                  t = 1/6,
                  e = 1/5
                  )
  return(
    structure(
      list(
        name = name,
        kernel = kernel,
        hessian = hessian,
        l2norm = l2norm,
        sigma2 = sigma2
      ),
      class = "CompStatKernel"
    )
  )
}

#' Plotting function for CompStatKernel's
#'
#' @param kernel A CompStatKernel
#' @param x A numeric vector with the points the kernel is evaluated in
#'
#' @return A ggplot showing the kernel
#' @export
#'
#' @examples
#' plot(CompStatKernel("g"))
#' plot(CompStatKernel("e"))
plot.CompStatKernel <- function(kernel, x = seq(-2,2,0.01)){
  ggplot2::ggplot(NULL, mapping = ggplot2::aes(x = x, y = kernel$kernel(x))) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "x", y = "K(x)")
}

#' Printing method for CompStatKernel's
#'
#' @param kernel
#'
#' @return Prints a tibble of information on the CompStatKernel
#' @export
#'
#' @examples
#' "g" %>% CompStatKernel()
#' "e" %>% CompStatKernel()
print.CompStatKernel <- function(kernel){
  print("CompStatKernel", quote = F)
  print(paste0("Type: ", kernel$name), quote = F)
  print(paste0("L2: ", round(kernel$l2norm, 4)), quote = F)
  print(paste0("Sigma2: ", round(kernel$sigma2, 4)), quote = F)
}



#' Fast implementation of the Epanechnikov kernel
#'
#' It is possible to make a comparable version of this function where x itself is used to store the values.
#' This is slightly faster if the kernel is evaluated many time on small inputs, but otherwise slightly slower
#'
#' @param x Vector for the kernel to be evaluated on
#'
#' @return The values of the kernel for the given input
epa_kernel <- function(x){
  res <- numeric(length(x))
  ind <- abs(x) <= 1
  res[!ind] <- 0
  res[ind] <- 3/4 * (1-x[ind]^2)
  res
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
#' set.seed(NULL)
#' grid <- density(x, 0.2, kernel = "g")$x
#'
#' microbenchmark::microbenchmark(
#' kernel_density("e", x, 0.2, method="c", grid = grid),
#' kernel_density("e", x, 0.2, method="r", grid = grid),
#' density(x, 0.2 * sqrt(CompStatKernel("e")$sigma2), kernel="e")
#' )
#'
#' microbenchmark::microbenchmark(
#' kernel_density("g", x, 0.2, method="c"),
#' kernel_density("g", x, 0.2, method="r"),
#' density(x, 0.2, kernel="g")
#' )
#'
#'
kernel_density <- function(kernel, x, bw,
                           grid = NULL, n = 512, from = NULL, to = NULL, cut = 3, method = "c",
                           return_grid = F
                           ){

  kernel <- CompStatKernel(kernel)

  if (is.null(grid)){
    rg <- range(x)
    if (is.null(from)){
      from <- rg[1] - cut*bw
    }
    if (is.null(to)){
      to <- rg[2] + cut*bw
    }
    grid <- seq(from, to, length.out=n)
  }

  if (return_grid){
    return(grid)
  }

  if (tolower(method) == "r"){
    return(eval_kdensR(kernel$kernel, grid, x, bw))
  }
  if (tolower(method) == "c"){
    if (kernel$name %in% c("Gaussian", "Epanechnikov")){
      return(eval_kdensC(kernel$name, grid, x, bw))
    } else {
      stop("Method 'c' is only implemented for Gaussian and Epanechnikov kernels.")
    }
  }

  stop("Must specify 'method' as either 'r' or 'c'.")

}

#' Benchmark the kernel density evaluation implementations
#'
#' @return A data frame containing benchmark information
#' @export
#'
#' @examples
#' bm <- benchmark_kernel_density()
#' bm %>%
#' bm %>%
#'  ggplot2::ggplot(ggplot2::aes(x = log(N), y = log(Time), color = Method)) +
#'  ggplot2::geom_line() +
#'  ggplot2::geom_point() +
#'  ggplot2::facet_wrap(~Kernel)
benchmark_kernel_density <- function(N_seq = 3:12, bw = 0.2){

  set.seed(0)
  xs <- N_seq %>% purrr::map(.f = function(n){rnorm(2^n)})
  set.seed(NULL)

  bm <- N_seq %>%
    purrr::imap_dfr(.f = function(n, i){
      purrr::map_dfr(.x = c("gaussian", "epanechnikov"), .f = function(kernel){
        calls <- list(
          call("kernel_density", kernel = kernel, x = xs[[i]], bw=bw, method = "c"),
          call("kernel_density", kernel = kernel, x = xs[[i]], bw=bw, method = "r"),
          call("density", kernel = kernel, x = xs[[i]], bw=bw*sqrt(CompStatKernel(kernel)$sigma2))
        )
        names(calls) <- c("cpp", "r", "stats::density")
        microbenchmark::microbenchmark(
          list=calls
        ) %>%
          as.data.frame() %>%
          dplyr::group_by(expr) %>%
          dplyr::summarise(Time = median(time)) %>%
          dplyr::mutate(Method = expr, Kernel = kernel, N = 2^n) %>%
          dplyr::select(-c(expr))
      })
    })
  bm
}

#' Implementation of kernel density evaluation in R.
#'
#' @param k The kernel evaluation function
#' @param grid The grid of values to evaluate the estimated density in
#' @param x The data points used to estimate the density
#' @param bw The bandwidth
#'
#' @return The estimated density evaluated at the grid points
#'
eval_kdensR <- function(k, grid, x, bw){
  grid %>%
    purrr::map_dbl(.f=function(z){
      mean(k((z-x)/bw))/bw
    })
}

l2norm <- function(f, ...){
  UseMethod("l2norm")
}

cvscore <- function(kernel, x, r, ...){
  UseMethod("cvscore")
}

#' Calculate an estimate of the L2-norm for the current kernel approximation
#'
#' @param x The data used for estimation
#' @param kernel_code The name of the kernel used
#' @param r The bandwidth
#' @param method A string. Specifies the method used to calculate the L2-Norm.
#'
#' @return An estimate of the L2-norm for the current kernel approximation
#' @export
#'
#' @examples
#' set.seed(0)
#' x <- rnorm(100)
#' H <- CompStatKernel("g")
#' l2norm(H, x, 1)
#' set.seed(NULL)
l2norm.CompStatKernel <- function(kernel, x, r, method = "default"){
  result <- NULL
  if (kernel$name == "Epanechnikov"){
    result <- switch(tolower(method),
           matrix=epanechnikov_l2norm_matrix(x, r),
           default=,
           c=epanechnikov_l2norm_runningC(x, r),
           running=epanechnikov_l2norm_running(x, r),
           binning=epanechnikov_l2norm_binning(x, r),
           NULL
           )
  } else if (kernel$name == "Gaussian"){
    result <- switch(tolower(method),
                     default=,
                     matrix=gaussian_l2norm_matrix(x, r),
                     integrate=gaussian_l2norm_integrate(x, r),
                     NULL
    )
  }
  if (is.null(result)){
    stop("No implementation for the supplied method/kernel combination is available.")
  }
  return(result)
}

#' Compute cross-validated out-of-sample loglikelihood for kernel smoother
#'
#' @param x A numeric vector of data points
#' @param bw The desired bandwidth
#' @param kernel_code The abbreviation for the kernel to use for smoothing
#' @param k The number of folds for cross validation
#' @param ... Additional arguments passed to kernel_density.
#'
#' @return The negative out-of-sample loglikelihood for the kernel smoother
#' @export
#'
#' @examples
#' x <- rnorm(1000)
#' bw <- seq(0.01, 1, 0.01)
#' g <- CompStatKernel("g")
#' scores <- bw %>% purrr::map_dbl(.f = cvscore, kernel = g, x = x)
#' plot(bw, scores, type = "l")
cvscore.CompStatKernel <- function(kernel, x, bw, k = 5, ...){
  set.seed(0)
  splits <- sample(1:k, length(x), replace = T)
  set.seed(NULL)
  scores <- numeric(k)
  for (i in 1:k){
    scores[i] <- -sum(log(kernel_density(
      kernel,
      x=x[splits != i],
      bw=bw,
      grid=x[splits == i],
      ...
      )))
  }
  mean(scores)
}
