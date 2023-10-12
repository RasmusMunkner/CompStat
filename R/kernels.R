
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
#' CompStatKernel("gaussian")
#' CompStatKernel("g")
#' CompStatKernel("t")
CompStatKernel <- function(kernel_code){

  if("CompStatKernel" %in% class(kernel_code)){
    return(kernel_code)
  }

  kernel_code <- kernel_code %>% substr(1,1) %>% tolower()
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

#' Fast implementation of the Epanechnikov kernel
#'
#' It is possible to make a comparable version of this function where x itself is used to store the values.
#' This is slightly faster if the kernel is evaluated many time on small inputs, but otherwise slightly slower
#'
#' @param x Vector for the kernel to be evaluated on
#'
#' @return The values of the kernel for the given input
#' @export
#'
#' @examples
#' epa_kernel(seq(-2,2,0.1))
epa_kernel <- function(x){
  res <- numeric(length(x))
  ind <- abs(x) <= 1
  res[!ind] <- 0
  res[ind] <- 3/4 * (1-x[ind]^2)
  res
}

gauss <- CompStatKernel("g")

microbenchmark::microbenchmark(
  gauss$kernel(u),
  mock_dnorm(u)
)

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
#' grid <- seq(-3,3,0.01)
#' kern <- CompStatKernel("g")
#' dens_hatR <- density(kern, x, 1, method = "r", grid = grid)
#' dens_hatC <- density(kern, x, 1, method = "c", grid = grid)
#' max(abs(dens_hatR - dens_hatC))
#' ggplot() +
#'   geom_histogram(aes(x = x, y = after_stat(density))) +
#'   geom_line(aes(x = z, y = dens_hat), color = "red")
#' set.seed(NULL)
#'
#' microbenchmark::microbenchmark(
#' kernel_density("e", x, 0.2, method="c"),
#' kernel_density("e", x, 0.2, method="r")
#' )
#'
#' u <- rnorm(100000)
#' profvis::profvis({
#' kernel_density("g", u, 0.2, method="r")
#' kernel_density("e", u, 0.2, method="r")
#' })
#'
kernel_density <- function(kernel, x, bw, grid = NULL, n = 512, from = NULL, to = NULL, cut = 3, method = "c"){

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

  if (method == "r"){
    return(eval_kdensR(kernel$kernel, grid, x, bw))
  }
  if (method == "c"){
    if (kernel$name %in% c("Gaussian", "Epanechnikov")){
      return(eval_kdensC(kernel$name, grid, x, bw))
    } else {
      stop("Method 'c' is only implemented for Gaussian and Epanechnikov kernels.")
    }
  }

  stop("Must specify 'method' as either 'r' or 'c'.")

}

#' Implementation of kernel density evaluation in R.
#'
#' @param k The kernel evaluation function
#' @param grid The grid of values to evaluate the estimated density in
#' @param x The data points used to estimate the density
#' @param bw The bandwidth
#'
#' @return The estimated density evaluated at the grid points
#' @export
#'
#' @examples
#'
#' gkern <- CompStatKernel("g")
#' grid <- rnorm(1000)
#' x <- rnorm(1000)
#' microbenchmark::microbenchmark(
#' eval_kdensR(gkern$kernel, grid, x, 0.1),
#' eval_kdensC("g", grid, x, 0.1)
#' )
#'
#'
eval_kdensR <- function(k, grid, x, bw){
  grid %>%
    purrr::map_dbl(.f=function(z){
      mean(k((z-x)/bw))/bw
    })
}

l2norm <- function(x, ...){
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
#' scores <- bw %>% purrr::map_dbl(.f = density_cv_score, x = x)
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
