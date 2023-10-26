
#' Print function for envelopes
#'
#' @param enve An envelope
#'
#' @export
print.Envelope <- function(enve){
  print(paste0("Envelope of type: ", enve$name))
  print(enve$base_rv)
}

#' Plot an Envelope with the function is is an envelope to
#'
#' @param enve An envelope
#' @param grid The grid on which the envelope and function are plotted.
#' Defaults to the provided support of the function.
#' @param logscale Should the envelope be plotted on log-scale?
#' @param compare If true, plot the ratio between the envelope and
#' base density. Note that on log-scale, this is simply the difference.
#'
#' @export
#'
#' @examples
#' enve <- GaussianEnvelope(named_rv("unif"))
#' plot(enve)
#' plot(enve, logscale=TRUE)
#' plot(enve, compare=T)
#' plot(enve, logscale=T, compare=T)
#'
plot.Envelope <- function(enve, grid = NULL, logscale=FALSE, compare = FALSE){

  if (is.null(grid)){
    grid <- seq(enve$support[1], enve$support[2], length.out = 512)
  }

  if (logscale){
    plot_data <- tibble::tibble(
      x = grid,
      f = enve$base_rv$log_f(grid),
      envelope = enve$log_f(grid)
    )
  } else {
    plot_data <- tibble::tibble(
      x = grid,
      f = enve$base_rv$f(grid),
      envelope = enve$f(grid)
    )
  }

  if (!compare){
    plot_data <- plot_data %>%
      tidyr::pivot_longer(-x, names_to = "Function", values_to = "y")

    ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = Function)) +
      ggplot2::geom_line()
  } else {
    if (logscale){
      ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = f-envelope)) +
        ggplot2::geom_line() +
        ggplot2::labs(y="Log(Envelope) - Log(f)")
    } else {
      ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = f / envelope)) +
        ggplot2::geom_line()
    }
  }

}


#' Laplacian Envelope Generator
#'
#' @param rv A random variable that the envelope should dominate
#' @param security Scaling used for alpha.
#' Safeguards against poor internal optimization.
#' @param optimize_log Boolean. Should the optimization be performed on log-scale?
#' @param sim_method Integer. Determines the simulation algorithm for the distribution.
#'
#' @return An envelope of the target RV
#' @export
#'
#' @examples
#' enve <- LaplaceEnvelope(named_rv("normal"))
#' plot(enve)
#'
#' enve <- LaplaceEnvelope(named_rv("poisson prior"))
#' plot(enve)
LaplaceEnvelope <- function(rv, security = 1, optimize_log = T, sim_method = 1){

  lap_rv <- named_rv("laplace")

  if (optimize_log){
    infimum <- optimise(function(x) lap_rv$log_f(x) - rv$log_f(x),
                        interval = rv$support
    )$objective
    infimum <- exp(infimum)
  } else {
    infimum <- optimise(function(x) lap_rv$f(x) / rv$f(x),
                        interval = rv$support
    )$objective
  }


  if (sim_method == 1){
    sim <- function(n) {rexp(n) - rexp(n)}
  }
  else if (sim_method == 2){
    sim <- function(n){
      z <- runif(n)
      (z < 1/2) * (log(2*z)) + (z >= 1/2) * (-log(1-2*(z-1/2)))
    }
  } else {
    stop("Invalid sim method. Only 1 and 2 are implemented for laplacian envelope.")
  }

  envelope <- structure(
    list(
      base_rv = rv,
      sim = sim,
      alpha = 1 / security,
      f = function(x) exp(-abs(x)) / 2 / infimum,
      log_f = function(x) - abs(x) - log(2) - log(infimum),
      f_prime = function(x) - sign(x) * exp(-abs(x)) / 2 / infimum,
      log_f_prime = function(x) - sign(x),
      name="Laplacian Envelope",
      support = c(-2,2)
    ),
    class = c("Envelope", "RandomVariable")
  )

  envelope
}

#' Compute a Gaussian Envelope
#'
#' @inheritParams LaplaceEnvelope
#'
#' @details
#' Note that optimization for alpha breaks if there are multiple local minima
#' This will never be the case for log-concave densities. See the examples for
#' a case where the optimization leads to non-coverage.
#'
#'
#' @return A Gaussian envelope of f.
#' @export
#'
#' @examples
#' # Coverage
#' enve <- GaussianEnvelope(named_rv("uniform"), security = 1)
#' plot(enve)
#'
#' # Poisson prior
#' # Warnings seem good-natured
#' enve <- GaussianEnvelope(named_rv("poisson prior"), optimize_log = T)
#' plot(enve)
#' plot(enve, logscale = T)
#'
#' # Non-coverage
#' enve <- GaussianEnvelope(named_rv("mix"), security = 1)
#' plot(enve)
GaussianEnvelope <- function(rv, security = 1, optimize_log = T, sim_method = 1){

  # Define a function that for a given choice of parameters (mu, sigma) finds
  # the optimal alpha that ensures the gaussian(mu, sigma) is an envelope
  if (optimize_log){
    get_alpha <- function(param){
      exp(optimise(function(x) dnorm(x, param[["mu"]], param[["sigma"]], log = T) - rv$log_f(x),
               interval = rv$support
      )$objective)
    }
  } else {
    get_alpha <- function(param){
      optimise(function(x) dnorm(x, param[["mu"]], param[["sigma"]]) / rv$f(x),
               interval = rv$support
      )$objective
    }
  }

  # Optimize to get the best (mu, sigma) in order to get the largest alpha
  best_params <- optim(c(mu=0, sigma = 1),
                       fn = get_alpha,
                       control=list(fnscale=-1)
                       )

  # Scaling the found alpha with the security in order to ensure envelope dominance
  best_params[["value"]] <- best_params[["value"]] / security

  envelope <- structure(
    list(
      base_rv = rv,
      sim = function(n){rnorm(n, mean=best_params$par[["mu"]], sd=best_params$par[["sigma"]])},
      alpha = best_params[["value"]],
      f = function(x) dnorm(x, mean=best_params$par[["mu"]], sd=best_params$par[["sigma"]]) / best_params$value,
      log_f = function(x) dnorm(x, mean=best_params$par[["mu"]], sd=best_params$par[["sigma"]], log = T) - log(best_params$value),
      f_prime = function(x) stop("Not implemented gaussian density derivative."),
      log_f_prime = function(x) stop("Not implemented log gaussian density derivative."),
      name = paste0("Gaussian(", round(best_params$par[["mu"]], 2), ",", round(best_params$par[["sigma"]]^2,2), ") Envelope"),
      support = best_params$par[["mu"]] + best_params$par[["sigma"]] * qnorm(0.975) * c(-1,1)
    ),
    class = c("Envelope", "RandomVariable")
  )

  envelope
}


