
#' Access commonly used kernels
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
#' get_kernel("gaussian")
#' get_kernel("g")
#' get_kernel("t")
get_kernel <- function(kernel_code){
  key <- kernel_code %>% substr(1,1) %>% tolower()
  kernel <- switch(key,
         n= ,
         g= dnorm,
         u= ,
         r= function(z){ifelse(abs(z)<1, 1/2, 0)},
         t= function(z){ifelse(abs(z)<1, 1-abs(z), 0)},
         e= function(z){ifelse(abs(z)<1, 3/4 * (1-z^2), 0)}
  )
  hessian <- switch(key,
         n= ,
         g= function(z){z^2 * dnorm(z) - dnorm(z)},
         u= ,
         r= function(z){rep(0, length(z))},
         t= function(z){rep(0, length(z))},
         e= function(z){ifelse(abs(z) < 1, -3/2, 0)}
  )
  l2norm <- switch(key,
                   n = ,
                   g = 1/(2*sqrt(pi)),
                   u = ,
                   r = 1/2,
                   t = 2/3,
                   e = 16/15
                   )
  sigma2 <- switch(key,
                  n = ,
                  g = 1,
                  u = ,
                  r = 1/3,
                  t = 1/6,
                  e = 4/15
                  )
  return(
    structure(
      list(
        kernel = kernel,
        hessian = hessian,
        l2norm = l2norm,
        sigma2 = sigma2
      ),
      class = "kernel"
    )
  )
}

#' Calculate an estimate of the L2-norm for the current kernel approximation
#'
#' @param x The data used for estimation
#' @param kernel_code The name of the kernel used
#' @param r The bandwidth
#'
#' @return An estimate of the L2-norm for the current kernel approximation
#' @export
#'
#' @examples
#' estimate_l2norm(1, "e", 1)
#' estimate_l2norm(rnorm(100), "r", 0.9 * 100^(1/5))
estimate_l2norm <- function(x, kernel_code, r){
  n <- length(x)
  H <- get_kernel(kernel_code)
  1/(n^2 * r^6) * integrate(
    Vectorize(function(z) sum(H$hessian((z-x)/r))^2),
    lower = -Inf,
    upper = Inf
  )$value
}











