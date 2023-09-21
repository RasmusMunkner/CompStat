
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
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
LogLinearEnvelope <- function(f, f_prime, tangent_points){

  # Checks if the tangent points are valid
  # They must all reside within the support of f
  # Could potentially be generalized, but hardly seems worth the effort
  if (sum(f(tangent_points) <= 0) > 0){
    stop("The density f is nonpositive at one or more of the supplied points.")
  }

  # Absolutely necessary quantities
  tangent_points <- sort(tangent_points)
  a <- f_prime(tangent_points) / f(tangent_points)
  b <- log(f(tangent_points)) - a * tangent_points
  z <- diff(-b) / diff(a)

  # Extra quantities that are nice for simulation
  R <- exp(b) / a * (exp(a*c(z, Inf)) - exp(a*c(-Inf, z)))
  R <- dplyr::coalesce(R, exp(b) * (c(z,Inf) - c(-Inf, z)))
  Q <- cumsum(R)
  c <- Q[length(Q)]

  # Evaluation and simulation
  eval_envelope <- function(x){
    section <- cut(x, breaks=c(-Inf, envelope$z, Inf))
    exp(envelope$a[section] * x + envelope$b[section])
  }

  # Quantile function
  quant <- function(p){
    Q <- c(0, Q)
    i <- cut(Q[length(Q)]*p, Q)
    Fx <- Q[length(Q)]*p - Q[as.integer(i)]
    x <- log(a[i]*Fx / exp(b[i]) + exp(a[i]*c(-Inf, z)[i])) / a[i]

    # Handling cases with a[i] = 0
    x %>%
      dplyr::coalesce(Fx / exp(b[i]) + c(-Inf, z)[i])
  }

  # Simulation
  sim <- function(n){
    u <- runif(n)
    quant(u)
  }

  envelope <- list(
    a=a, b=b, z=z, R = R, Q = Q,
    g = eval_envelope,
    f = f,
    quant = quant,
    sim = sim,
    alpha = 1/c
  )
  class(envelope) <- c("LogLinearEnvelope", "Envelope")
  return(envelope)
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
#' plot(enve, logscale=TRUE)
plot.Envelope <- function(enve, grid = seq(-6, 6, 0.001), logscale=FALSE){

  if (logscale){
    plot_data <- tibble::tibble(
      x = grid,
      f = log(enve$f(grid)),
      envelope = log(enve$g(grid))
    )
  } else {
    plot_data <- tibble::tibble(
      x = grid,
      f = enve$f(grid),
      envelope = enve$g(grid)
    )
  }
  plot_data <- plot_data %>%
    tidyr::pivot_longer(-x, names_to = "what", values_to = "y")

  ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = what)) +
    ggplot2::geom_line()
}
