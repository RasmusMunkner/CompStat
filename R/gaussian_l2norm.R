
gaussian_l2norm_matrix <- function(x, r, debug = F){
  if (debug){
    browser()
  }

  delta <- outer(x/2, x/2, FUN = `-`)
  q <- r/sqrt(2)
  a0 <- (delta^2 - r^2)^2
  a2 <- -2*(delta^2 + r^2)
  c <- r / sqrt(4*pi) * exp(-delta^2 / r^2)

  m <- c * (3*q^4 + a2*q^2 + a0) / r^10 / length(x)^2

  sum(m)

}

gaussian_l2norm_integrate <- function(x, r){
  H <- get_kernel("g")
  n <- length(x)
  result <- 1/(n^2 * r^6) * integrate(
    Vectorize(function(z) sum(H$hessian((z-x)/r))^2),
    lower = -Inf,
    upper = Inf,
    subdivisions = 500
  )$value
  result
}
