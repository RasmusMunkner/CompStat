
#' Laplacian Envelope Generator
#'
#' @param f The function for which to construct the envelope. Can also pass a RandomVariable.
#' @param rg The range used to compute the fraction of densities.
#' @param security Scaling used for alpha.
#' Safeguards against poor internal optimization.
#'
#' @return A list of type 'Envelope'
#' @export
#'
#' @examples
#' enve <- LaplaceEnvelope(dnorm)
LaplaceEnvelope <- function(f, rg = c(-3,3), security = 1 + 1e-9, sim_method = 1){

  # Handling the input whether it is a density or a RV
  if ("RandomVariable" %in% class(f)){
    base_rv <- f
    f <- f$f
  } else {
    base_rv <- RandomVariable(f=f)
  }

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

  envelope <- structure(
    list(
      base_rv = base_rv,
      f = function(x) exp(-abs(x)) / 2 / infimum,
      sim = sim,
      alpha = 1 / security,
      log_f = function(x) - abs(x) - log(2) - log(infimum),
      f_prime = function(x) - sign(x) * exp(-abs(x)) / 2 / infimum,
      log_f_prime = function(x) - sign(x)
    ),
    class = c("Envelope", "RandomVariable")
  )

  envelope
}

#' Compute a Gaussian Envelope
#'
#' @inheritParams LaplaceEnvelope
#'
#' @return A Gaussian envelope of f.
#' @export
#'
#' @examples
#' enve <- GaussianEnvelope(dunif, rg=c(0,1), security = 1.1)
GaussianEnvelope <- function(f, rg = c(-3,3), security = 1+1e-9){

  # Handling the input whether it is a density or a RV
  if ("RandomVariable" %in% class(f)){
    base_rv <- f
    f <- f$f
  } else {
    base_rv <- RandomVariable(f=f)
  }

  # Finds the smallest proportion that needs to be used for the gaussian
  # envelope to actually be an envelope (alpha * f <= g)
  get_alpha <- function(param){
    optimise(function(x) dnorm(x, param[["mu"]], param[["sigma"]]) / f(x),
             interval = rg
    )$objective
  }

  best_params <- optim(c(mu=0, sigma = 1),
                       fn = get_alpha,
                       control=list(fnscale=-1)
                       )

  best_params[["value"]] <- best_params[["value"]] / security

  envelope <- structure(
    list(
      base_rv = base_rv,
      f = function(x) dnorm(x, mean=best_params$par[["mu"]], sd=best_params$par[["sigma"]]) / best_params$value,
      sim = function(n){rnorm(n, mean=best_params$par[["mu"]], sd=best_params$par[["sigma"]])},
      alpha = best_params[["value"]],
      log_f = function(x) -log(sqrt(2*pi)) - (x - best_params$par[["mu"]])^2 / 2 / best_params$par[["sigma"]]^2,
      f_prime = function(x) stop("Not implemented gaussian density derivative."),
      log_f_prime = function(x) stop("Not implemented log gaussian density derivative.")
    ),
    class = c("Envelope", "RandomVariable")
  )

  envelope
}




