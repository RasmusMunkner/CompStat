
#' Envelopes for log-concave distributions
#'
#' @param f A log-concave function.
#' @param f_prime The derivative of f.
#' @param tangent_points The points at which the envelope should be tanget to f.
#'
#' @return A LogLinearEnvelope object.
#' @export
#'
#' @examples
#' LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
LogLinearEnvelope <- function(f, f_prime, tangent_points){

  tangent_points <- sort(tangent_points)
  a <- f_prime(tangent_points) / f(tangent_points)
  b <- log(f(tangent_points)) - a * tangent_points
  cutpoints <- diff(-b) / diff(a)
  envelope <- list(a=a, b=b,cutpoints=cutpoints, f = f)
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
eval_envelope1 <- function(x, envelope){

  if(!("LogLinearEnvelope" %in% class(envelope))){
    stop("The provided envelope is not of class 'LogLinearEnvelope'")
  }

  section <- as.integer(cut(x, breaks=c(-Inf, envelope$cutpoints, Inf)))

  list(x = x, i = section) %>%
    purrr::pmap_dbl(.f=function(x,i){
      envelope$a[i]*x + envelope$b[i]
    })
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
eval_envelope2 <- function(x, envelope){

  if(!("LogLinearEnvelope" %in% class(envelope))){
    stop("The provided envelope is not of class 'LogLinearEnvelope'")
  }

  section <- as.integer(cut(x, breaks=c(-Inf, envelope$cutpoints, Inf)))
  x_part <- split(x, section)

  seq_along(x_part) %>%
    purrr::map(.f=function(i){
      envelope$a[i] * x_part[[i]] + envelope$b[i]
    }) %>%
    purrr::reduce(c)
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
#' plot(enve)
plot.LogLinearEnvelope <- function(enve, grid = seq(-6, 6, 0.001), logscale=FALSE){

  if (logscale){
    plot_data <- tibble::tibble(
      x = grid,
      f = log(enve$f(grid)),
      envelope = eval_envelope2(grid, enve)
    )
  } else {
    plot_data <- tibble::tibble(
      x = grid,
      f = enve$f(grid),
      envelope = exp(eval_envelope2(grid, enve))
    )
  }
  plot_data <- plot_data %>%
    tidyr::pivot_longer(-x, names_to = "what", values_to = "y")

  ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = what)) +
    ggplot2::geom_line()
}
