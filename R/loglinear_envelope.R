
#' Envelopes for log-concave distributions
#'
#' @param rv An object created by RandomVariable.
#' @param tangent_points The points at which the envelope should be tangent to f.
#'
#' @return A LogLinearEnvelope object.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(get_rv(), c(-2,0,1))
#' enve <- LogLinearEnvelope(get_rv("n"))
#' enve <- LogLinearEnvelope(get_rv("g"))
LogLinearEnvelope <- function(rv, tangent_points = NULL){

  if (is.null(tangent_points)){
    #Tries to find the
    #browser()
    tangent_points <- uniroot(rv$log_f_prime, interval = c(-5, 5), extendInt = "y")$root + c(-1,1)
  }

  # Checks if the tangent points are valid
  # They must all reside within the support of f
  # Could potentially be generalized, but hardly seems worth the effort
  if (sum(rv$f(tangent_points) <= 0) > 0){
    stop("The density f is nonpositive at one or more of the supplied points.")
  }

  # Absolutely necessary quantities
  tangent_points <- sort(tangent_points)
  a <- rv$log_f_prime(tangent_points)
  b <- rv$log_f(tangent_points) - a * tangent_points
  z <- diff(-b) / diff(a)

  # Extra quantities that are nice for simulation
  R <- exp(b) / a * (exp(a*c(z, Inf)) - exp(a*c(-Inf, z)))
  R <- dplyr::coalesce(R, exp(b) * (c(z,Inf) - c(-Inf, z)))
  Q <- cumsum(R)
  c <- Q[length(Q)]
  eb <- exp(b)
  eaz <- exp(a * c(-Inf, z))

  # Evaluation and simulation
  eval_envelope <- function(x){
    section <- cut(x, label = FALSE, breaks=c(-Inf, z, Inf))
    a[section] * x + b[section]
  }

  # Quantile function (This could potentially be improved a lot)
  quant <- function(p){
    Q <- c(0, Q)
    i <- cut(Q[length(Q)]*p, Q, labels=FALSE)
    Fx <- Q[length(Q)]*p - Q[i]
    aF <- a[i]*Fx
    aFb <- aF / exp(b[i])
    aFb <- aF / eb[i]
    zTrunc <- c(-Inf, z)[i]
    aZ <- exp(a[i]*zTrunc)
    aZ <- eaz[i]
    almost <-  log(aFb + aZ)
    x <- almost / a[i]

    # Handling cases with a[i] = 0
    x %>%
      dplyr::coalesce(Fx / exp(b[i]) + c(-Inf, z)[i])
  }

  # Simulation
  sim <- function(n){
    u <- runif(n)
    quant(u)
  }

  envelope_rv <- RandomVariable(log_f = eval_envelope)

  envelope <- structure(
    list(
      a=a, b=b, z=z, R = R, Q = Q, #LogLinearEnvelope Fields
      tangent_points = tangent_points,
      quant = quant,
      base_rv = rv, # Envelope fields
      sim = sim,
      alpha = 1/c,
      f = envelope_rv$f, # RandomVariable fields
      log_f = envelope_rv$log_f,
      f_prime = envelope_rv$f_prime,
      log_f_prime = envelope_rv$log_f_prime
      ),
    class = c("LogLinearEnvelope", "Envelope", "RandomVariable")
  )
  return(envelope)
}

#' Re-estimate a LogLinearEnvelope
#'
#' Mostly useful for adding additional tangent points.
#'
#' @param enve The previous envelope.
#' @param tangent_points The tangent points the new envelope should be fitted to.
#'
#' @return A LogLinearEnvelope.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(get_rv("n"))
#' plot(enve)
#' enve <- update.LogLinearEnvelope(enve, c(enve$tangent_points, 0))
#' plot(enve)
update.LogLinearEnvelope <- function(enve, tangent_points){
  enve <- LogLinearEnvelope(enve$base_rv, tangent_points)
  return(enve)
}

# This is black magic needed by update.LogLinearEnvelope
# Truth be told, the function just needs to pretend to have a formula, like y~x
# Which it does not, but we can hack it (be careful though)
getCall.LogLinearEnvelope <- function(enve){
  NULL
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
#' enve <- LogLinearEnvelope(get_rv(), c(-2,0,1))
#' plot(enve)
#' plot(enve, logscale=TRUE)
plot.Envelope <- function(enve, grid = seq(-6, 6, 0.001), logscale=FALSE){

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
  plot_data <- plot_data %>%
    tidyr::pivot_longer(-x, names_to = "what", values_to = "y")

  ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = what)) +
    ggplot2::geom_line()
}
