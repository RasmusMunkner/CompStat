
#' Laplacian Envelope Generator
#'
#' @param f The function for which to construct the envelope.
#' @param rg The range used to compute the fraction of densities.
#'
#' @return A list of type 'Envelope'
#' @export
#'
#' @examples
#' enve <- LaplaceEnvelope(dnorm)
LaplaceEnvelope <- function(f, rg = c(-3,3), sim_method = 1){
  infimum <- optimise(function(x) exp(-abs(x)) / 2 / f(x),
                      interval = rg
                      )$objective

  if (sim_method == 1){
    sim <- function(n) {rexp(n) - rexp(n)}
  }
  else if (sim_method == 2){
    sim <- function(n){
      z <- runif(n)
      (z < 1/2) * (log(2*z)) + (z >= 1/2) * (-log(1-2*(z-1/2)))
    }
  }

  envelope <- list(
    f = f,
    g = function(x) exp(-abs(x)) / 2 / infimum,
    sim = sim
  )

  class(envelope) <- c("Envelope")

  envelope
}








