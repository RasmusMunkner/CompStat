#' An object mean to represent a distribution.
#'
#' A RandomVariable object should expose a density and a log-density.
#'
#' @param f A density function.
#' @param log_f The logarithm of the density function.
#' @param f_prime The derivative of the density.
#' @param log_f_prime The derivative of the log-density.
#'
#' Note that the log-versions of the density and the differentiated density can
#' can be inferred from one another. This is done only if NULL is passed.
#' Any function that could not be inferred will be replaced with an exception raiser.
#'
#' @return A RandomVariable object.
#' @export
#'
#' @examples
#' X <- RandomVariable(dnorm, function(z) log(dnorm(z)), function(z) function(z){dnorm(z) * (-z)})
#' P1 <- RandomVariable(log_f = log_rejec_dens5,log_f_prime = log_rejec_dens_prime)
#' P_fail <- RandomVariable(f = dnorm)
RandomVariable <- function(f = NULL, log_f = NULL, f_prime = NULL, log_f_prime = NULL, name = NULL, support = NULL){

  if (length(c(f, log_f, f_prime, log_f_prime)) == 0){
    print(c(f, log_f, f_prime, log_f_prime))
    stop("Cannot define a RandomVariable without passing any information.")
  }

  # Handling f
  if (is.null(f)){
    if (is.null(log_f)){
      warning("A density could not be inferred. Use this RV with caution!")
      f <- function(z) stop("The random variable was created without specifying a density.")
    } else {
      f <- function(z) exp(log_f(z))
    }
  }

  # Handling log_f
  if (is.null(log_f)){
    if (is.null(f)){
      warning("A log-density could not be inferred. Use this RV with caution!")
      log_f <- function(z) stop("The random variable was created without specifying a log-density")
    } else {
      log_f <- function(z) log(f(z))
    }
  }

  # Handling f_prime
  if (is.null(f_prime)){
    if (is.null(log_f_prime)){
      warning("A derivative of the density could not be inferred. Use this RV with caution!")
      f_prime <- function(z) stop("The random variable was created without specifying derivative.")
    } else {
      f_prime <- function(z) exp(log_f_prime(z))
    }
  }

  # Handling log_f_prime
  if (is.null(log_f_prime)){
    if (is.null(f_prime) || is.null(f)){
      warning("A density of the log-density could not be inferred. Use this RV with caution!")
      log_f_prime <- function(z) stop("The random variable was created without specifying a differentiated log-density")
    } else {
      log_f_prime <- function(z) f_prime(z) / f(z)
    }
  }

  # Tagging the RV
  if (is.null(name)){
    message("No name supplied for RV creation. Consider naming custom RV's.")
    name <- "Unnamed custom RV"
  }

  # Determine the support
  if (is.null(support)){
    message("No support was specified. Using the default [-2,2].")
    support <- c(-2,2)
  }

  # Put the thing together
  return (
    structure(list(f = f, log_f = log_f, f_prime = f_prime, log_f_prime = log_f_prime, name = name, support = support), class = "RandomVariable")
  )
}

#' Print functionality for RV
#'
#' @param x A random variable
#'
#' @export
print.RandomVariable <- function(rv){
  print(paste0("RV of type: ", rv$name))
}

#' Density plot for random variables
#'
#' @param rv A RandomVariable
#' @param x A grid of points for which to plot the density
#'
#' @return A plot of the density of the rv
#' @export
#'
#' @examples
#' X1 <- named_rv("n")
#' X2 <- named_rv("g")
#' X3 <- named_rv("p")
#' plot(X1)
#' plot(X2, range = c(0,3))
#' plot(X3, range = c(0, 2/3), logscale = T)
plot.RandomVariable <- function(rv, range = NULL, logscale = F, resolution = 100){
  if (is.null(range)){
    range <- seq(rv$support[1], rv$support[2], length.out = resolution)
  } else {
    range <- seq(range[1], range[2], length.out = resolution)
  }
  if (logscale){
    y <- rv$log_f(range)
    ylabel <- "log(f(x))"
  } else {
    y <- rv$f(range)
    ylabel <- "f(x)"
  }

  ggplot2::ggplot(mapping = ggplot2::aes(x = range, y = y)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "x", y = ylabel) +
    ggplot2::ggtitle(rv$name)
}

#' Convenience function for encoding common distributions into RV's.
#'
#' @param key The name (first letter) of the desired distribution.
#'
#' @return A RandomVariable object corresponding to the distribution specified.
#' @export
#'
#' @examples
#' named_rv("n")
#' named_rv("g")
#' named_rv("p")
named_rv <- function(key = "n"){
  key <- key %>% substr(1,1) %>% tolower()
  rv <- switch(key,
         "n"=RandomVariable(f = dnorm,
                            log_f = function(z) -log(sqrt(2*pi)) - z^2 / 2,
                            f_prime = function(z) dnorm(z) * (-z),
                            name = "Standard Gaussian",
                            support = c(-2,2)
                            ),
         "g"=RandomVariable(f = function(z) dgamma(z, 2),
                            log_f = function(z) log(z) - z,
                            f_prime = function(z) -dgamma(z,2) + exp(-z),
                            name = "Gamma(2,1)",
                            support = c(0,5)
                            ),
         "p"=RandomVariable(log_f = poisson_prior_log_f,
                            log_f_prime = poisson_prior_log_f_prime,
                            name = "Poisson prior",
                            support = c(0.001, 2/3)
                            ),
         "a"=RandomVariable(log_f = poisson_prior_log_f_approx,
                            log_f_prime = poisson_prior_log_f_prime,
                            name = "Approximate poisson prior",
                            support = c(0.001, 2/3)
                            ),
         "e"=RandomVariable(f = function(z) ifelse(abs(z) < 1, 3/4 * (1-z^2), 0),
                            log_f_prime = function(z) ifelse(abs(z) < 1, - 2*z / (1 - z^2), NaN),
                            name = "Epanechnikov",
                            support = c(-1 + 1e-9, 1-1e-9) # Support is the open interval
                            ),
         "l"=RandomVariable(f = function(z) exp(-abs(z))/2,
                            log_f = function(z) -abs(z) - log(2),
                            log_f_prime = function(z) -sign(z),
                            name = "Laplace",
                            support = c(-2,2)
                            ),
         "u"=RandomVariable(f= function(z) ifelse(0 < z & z < 1, 1, 0),
                            log_f= function(z) ifelse(0 < z & z < 1, 0, -Inf),
                            log_f_prime = function(z) 0,
                            name = "Uniform",
                            support = c(0,1)
                            ),
         "m"=RandomVariable(f = function(z) dnorm(z) / 3 + 2 * dnorm(z, 4) / 3,
                            log_f_prime = function(z) (dnorm(z) * (-z/2) / 3 + 2*dnorm(z, 4) * (-(z-4)/2) / 3) / (dnorm(z) / 3 + 2 * dnorm(z, 4) / 3),
                            name = "Gaussian(0,1)/Gaussian(4,1) (1/3)/(2/3)-Mixture",
                            support = c(-3,8)
                            )
         )
  rv
}




