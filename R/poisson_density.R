data("poisson")


#' Density evaluation for rejection sampling
#'
#' @param y Evaluation point(s).
#'
#' @return The density evaluated in the given points.
#' @export
#'
#' @examples
#' rejec_dens1(0)
rejec_dens1 <- Vectorize(function(y){
  xy <- y * poisson$x
  prod(exp(xy * poisson$z - exp(xy)))
})

#' Density evaluation for rejection sampling
#'
#' @param y Evaluation point(s).
#'
#' @return The density evaluated in the given points.
#' @export
#'
#' @examples
#' rejec_dens2(0)
rejec_dens2 <- Vectorize(function(y){
  xy <- y * poisson$x
  exp(sum(xy * poisson$z - exp(xy)))
})


precomputed_xz <- sum(poisson$x * poisson$z)
#' Density evaluation for rejection sampling
#'
#' @param y Evaluation point(s).
#'
#' @return The density evaluated in the given points.
#' @export
#'
#' @examples
#' rejec_dens3(0)
rejec_dens3 <- Vectorize(function(y){
  exp(y*precomputed_xz - sum(exp(y*poisson$x)))
})

#' Density evaluation for rejection sampling
#'
#' @param y Evaluation point(s).
#'
#' @return The density evaluated in the given points.
#' @export
#'
#' @examples
#' rejec_dens4(0)
rejec_dens4 <- function(y){
  term1 <- y*precomputed_xz
  term2 <- rowSums(exp(outer(y, poisson$x)))
  exp(term1 - term2)
}

#' Density evaluation for rejection sampling
#'
#' @param y Evaluation point(s).
#'
#' @return The density evaluated in the given points.
#' @export
#'
#' @examples
#' rejec_dens5(0)
rejec_dens5 <- function(y){
  term1 <- y*sum(poisson$x * poisson$z)
  term2 <- rowSums(exp(outer(y, poisson$x)))
  exp(term1 - term2)
}






