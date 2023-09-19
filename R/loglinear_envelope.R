
#' Envelopes for log-concave distributions
#'
#' @param f A log-concave function.
#' @param f_prime The derivative of f.
#' @param tangent_points The points at which the envelope should be tangent to f.
#'
#' @return A LogLinearEnvelope object.
#' @export
#'
#' @examples
#' LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
LogLinearEnvelope <- function(f, f_prime, tangent_points){

  # Absolutely necessary quantities
  tangent_points <- sort(tangent_points)
  a <- f_prime(tangent_points) / f(tangent_points)
  b <- log(f(tangent_points)) - a * tangent_points
  cutpoints <- diff(-b) / diff(a)

  # Extra quantities that are nice for simulation
  R <- exp(b) / a * (exp(a*c(cutpoints, Inf)) - exp(a*c(-Inf, cutpoints)))
  R <- dplyr::coalesce(R, exp(b) * (c(cutpoints,Inf) - c(-Inf, cutpoints)))
  Q <- cumsum(R)
  c <- Q[length(Q)]

  envelope <- list(a=a, b=b,cutpoints=cutpoints, f = f, R = R, Q = Q, c = c)
  class(envelope) <- "LogLinearEnvelope"
  return(envelope)
}

#' Evalute a LogLinearEnvelope in a sequence of points
#'
#' This implementation is particularly slow.
#'
#' @param x The points to evaluate the envelope in.
#' @param envelope The envelope.
#'
#' @return The values of envelope(x).
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
#' eval_envelope1(seq(-3,3,0.1), enve)
eval_envelope1 <- function(x, envelope, logscale = FALSE){

  if(!("LogLinearEnvelope" %in% class(envelope))){
    stop("The provided envelope is not of class 'LogLinearEnvelope'")
  }

  section <- as.integer(cut(x, breaks=c(-Inf, envelope$cutpoints, Inf)))

  log_y <- list(x = x, i = section) %>%
    purrr::pmap_dbl(.f=function(x,i){
      envelope$a[i]*x + envelope$b[i]
    })

  if (logscale){
    log_y
  } else {
    exp(log_y)
  }

}

#' Evalute a LogLinearEnvelope in a sequence of points
#'
#' This implementation is relatively fast.
#'
#' @inheritParams eval_envelope1
#'
#' @return The values of envelope(x).
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
#' eval_envelope2(seq(-3,3,0.1), enve)
eval_envelope2 <- function(x, envelope, logscale=FALSE){

  if(!("LogLinearEnvelope" %in% class(envelope))){
    stop("The provided envelope is not of class 'LogLinearEnvelope'")
  }

  section <- as.integer(cut(x, breaks=c(-Inf, envelope$cutpoints, Inf)))
  x_part <- split(x, section)

  log_y <- seq_along(x_part) %>%
    purrr::map(.f=function(i){
      envelope$a[i] * x_part[[i]] + envelope$b[i]
    }) %>%
    purrr::reduce(c)

  if (logscale){
    log_y
  } else {
    exp(log_y)
  }

}

#' Plot a LogLinearEnvelope with the function is is an envelope to
#'
#' @param enve The LogLinearEnvelope.
#' @param grid The grid on which the envelope and function are plotted.
#'
#' @return A ggplot-object.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-3, -1, 1, 3))
#' sim_true
plot.LogLinearEnvelope <- function(enve, grid = seq(-6, 6, 0.001), logscale=FALSE){

  if (logscale){
    plot_data <- tibble::tibble(
      x = grid,
      f = log(enve$f(grid)),
      envelope = eval_envelope2(grid, enve, logscale)
    )
  } else {
    plot_data <- tibble::tibble(
      x = grid,
      f = enve$f(grid),
      envelope = eval_envelope2(grid, enve, logscale)
    )
  }
  plot_data <- plot_data %>%
    tidyr::pivot_longer(-x, names_to = "what", values_to = "y")

  ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = what)) +
    ggplot2::geom_line()
}

#' Quantile function for LogLinearEnvelope
#'
#' @param p A vector evaluation points, each with 0 < p < 1.
#' @param enve A LogLinearEnvelope.
#'
#' @return The quantiles of the envelope corresponding to p.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
#' p <- 1:99 / 100
#' q <- qLogLinearEnvelope(p, enve)
#' ggplot(mapping = aes(x = q, y = p)) +
#' geom_line() +
#' geom_line(aes(y = pnorm(x)), color = "red")
qLogLinearEnvelope <- function(p, enve){
  if (!("LogLinearEnvelope" %in% class(enve))){
    stop("The passed envelope is not of class 'LogLinearEnvelope'.")
  }

  Q <- c(0, enve$Q)
  i <- cut(Q[length(Q)]*p, Q)
  Fx <- Q[length(Q)]*p - Q[as.integer(i)]
  x <- log(enve$a[i]*Fx / exp(enve$b[i]) + exp(enve$a[i]*c(-Inf, enve$cutpoints)[i])) / enve$a[i]

  # Handling cases with a[i] = 0
  x %>%
    dplyr::coalesce(Fx / exp(enve$b[i]) + c(-Inf, enve$cutpoints)[i])
}

#' Simulate from a LogLinearEnvelope
#'
#' @param n The number of simulations.
#' @param enve A LogLinearEnvelope.
#'
#' @return A vector of the desired length containing simulations from the envelope.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
#' y <- rLogLinearEnvelope(10000, enve)
#' hist(y)
rLogLinearEnvelope <- function(n, enve){
  u <- runif(n)
  qLogLinearEnvelope(u, enve)
}
