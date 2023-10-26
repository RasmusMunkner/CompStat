
#' Envelopes for log-concave distributions
#'
#' @param rv A RandomVariable.
#' @param tangent_points The points at which the envelope should be tangent to
#' the RandomVariable log-density.
#' @param autoselection_msg Turn this off to suppress reporting of automatically
#' selected tangent points.
#' @param cache_eval Should the envelope use precomputed quantities for
#' simulation? (This option is only included to show an example of
#' poor optimization when it is turned off).
#'
#' @return A LogLinearEnvelope object.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(named_rv(), c(-2,0,1))
#' enve <- LogLinearEnvelope(named_rv("n"))
#' enve <- LogLinearEnvelope(named_rv("g"))
LogLinearEnvelope <- function(rv, tangent_points = NULL, autoselection_msg = T, precompute = T){

  if (is.null(tangent_points)){
    #Tries to find good tangent points.
    maxpoint <- uniroot(rv$log_f_prime, interval = rv$support, extendInt = "y")$root
    tangent_points <- c((maxpoint - rv$support[1])/2 + rv$support[1], rv$support[2]-(rv$support[2] - maxpoint)/2)
    if (autoselection_msg){
      message(paste0("Tangent points were not specified. Automatic selection (", round(tangent_points[1],2), ",", round(tangent_points[2],2), ") is used."))
    }
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

  # Safeguard against flat tangents (they slow down computation)
  tangent_points <- tangent_points[a != 0]
  a <- rv$log_f_prime(tangent_points) # Recompute a

  # Check that the envelope is well-defined
  if (a[1] <= 0 || a[length(a)] >= 0){
    stop("Envelope is not well-defined. Usually this is due to inadequate tangent point spread.")
  }

  b <- rv$log_f(tangent_points) - a * tangent_points
  z <- diff(-b) / diff(a)

  # Extra quantities that are nice for simulation
  R <- sign(a) * (exp(a*c(z, Inf) + b - log(abs(a))) - exp(a*c(-Inf, z) + b - log(abs(a))))
  Q <- cumsum(R)
  c <- Q[length(Q)]

  # Caching calculations for fast simulation
  eb_cache <- exp(b)
  aeb_cache <- a / eb_cache
  eaz_cache <- exp(a * c(-Inf, z))

  # Evaluation and simulation
  eval_envelope <- function(x){
    section <- .bincode(x, breaks=c(-Inf, z, Inf))
    a[section] * x + b[section]
  }

  # Quantile function (This could potentially be improved a lot)
  if (precompute){
    quant <- function(p){
      Q <- c(0, Q)
      i <- .bincode(Q[length(Q)]*p, Q)
      Fx <- Q[length(Q)]*p - Q[i]

      # Term 1
      aFb <- aeb_cache[i] * Fx

      # Term 2
      aZ <- eaz_cache[i]

      # Finally
      x <-  log(aFb + aZ) / a[i]
      x
    }
  } else {
    quant <- function(p){
      Q <- c(0, Q)
      i <- .bincode(Q[length(Q)]*p, Q)
      Fx <- Q[length(Q)]*p - Q[i]

      # Term 1
      aF <- a[i]*Fx
      aFb <- aF / exp(b[i])

      # Term 2
      aZ <- exp(a[i]*c(-Inf, z)[i])

      # Finally
      x <-  log(aFb + aZ) / a[i]
      x
    }
  }

  # Simulation
  sim <- function(n){
    u <- runif(n)
    quant(u)
  }

  envelope <- structure(
    list(
      a=a, b=b, z=z, R = R, Q = Q, #LogLinearEnvelope Fields
      tangent_points = tangent_points,
      quant = quant,
      base_rv = rv, # Envelope fields
      sim = sim,
      alpha = 1, #Niels says this should be 1/c, but we are already on unnormalized scale
      f = function(z) exp(eval_envelope(z)), # RandomVariable fields
      log_f = eval_envelope,
      f_prime = function(z) stop("No implementation of derivatives for log-linear-envelopes."),
      log_f_prime = function(z) stop("No implementation of derivatives for log-linear-envelopes."),
      name = "Log-Linear Envelope",
      support = rv$support
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
update.LogLinearEnvelope <- function(enve, tangent_points, project_on_support = T){

  if (project_on_support){
    supp_len <- diff(enve$base_rv$support) / 2
    for (i in 1:length(tangent_points)){
      if (tangent_points[i] <= enve$base_rv$support[1]){
        tangent_points[i] <- enve$base_rv$support[1] + 1e-14 + 1/(1/supp_len + enve$base_rv$support[1] - tangent_points[i])
      } else if (tangent_points[i] >= enve$base_rv$support[2]){
        tangent_points[i] <- enve$base_rv$support[2] - 1e-14 - 1/(1/supp_len + tangent_points[i]-enve$base_rv$support[2])
      }
    }
  }

  enve <- LogLinearEnvelope(enve$base_rv, tangent_points)
  return(enve)
}

# This is black magic needed by update.LogLinearEnvelope
# Truth be told, the function just needs to pretend to have a formula, like y~x
# Which it does not, but we can hack it (be careful though)
getCall.LogLinearEnvelope <- function(enve){
  NULL
}
