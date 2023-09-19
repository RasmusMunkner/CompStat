

#' Slow rejection sampling for LogLinearEnvelopes
#'
#' This method is abysmally slow and was only implemented to provide a baseline.
#'
#' #' Currently, something is clearly wrong with this implementation.
#' WORK IN PROGRESS.
#'
#' @param n The number of simulations.
#' @param enve A LogLinearEnvelope.
#' @param alpha
#'
#' @return Simulations from the distribution that the envelope is based on.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
#' rejection_sampler_naive(10, enve)
rejection_sample_naive <- function(n, enve){
  alpa <- 1/enve$c
  y_vec <- vector("numeric")
  n_sim <- 1
  while (n_sim <= n){
    U <- runif(1)
    Y <- rLogLinearEnvelope(1, enve)
    if (U <= alpha * enve$f(Y) / eval_envelope2(Y, enve)){
      y_vec[n_sim] <- Y
      n_sim <- n_sim + 1
    }
  }
  y_vec
}

#' Compile rejection samplers for LogLinearEnvelopes
#'
#' Currently, something is clearly wrong with this implementation.
#' WORK IN PROGRESS.
#'
#' @param enve A logLinearEnvelope \eqn{g}.
#' @param alpha A scalar value ensuring \eqn{\alpha f \leq g}.
#'
#' @return A rejection sampler for the given envelope.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
#' sampler <- rejection_sample_factory(enve, 1/enve$c)
#' sim <- sampler(5000)
#' hist(sim)
#' qqnorm(sim)
#' qqline(sim)
#' environment(sampler)$p
#' environment(sampler)$credibility
rejection_sample_factory <- function(enve, alpha = NULL){
  if (is.null(alpha)){
    alpha <- 1/enve$c
  }
  p <- alpha
  credibility <- 20 # Arbitrary
  rejection_sampler <- function(n, train = TRUE){
    sim <- vector("numeric", n)
    n_sim <- 0
    while (n_sim < n){

      # Calculate RV's
      U <- runif(ceiling((n - n_sim)/p))
      Y <- rLogLinearEnvelope(ceiling((n - n_sim)/p), enve)
      Y_accept <- Y[U <= alpha * enve$f(Y) / eval_envelope2(Y, enve)]

      # Store simulations
      if (length(Y_accept) >= 1){
        new_sim_length <- length(Y_accept) - max(0, n_sim + length(Y_accept) - n)
        sim[(n_sim+1):(n_sim+new_sim_length)] <- Y_accept[1:new_sim_length]
        n_sim <- n_sim + length(Y_accept)
      }


      # Update simulation parameters
      if (train){
        p <<- (p * credibility + length(Y_accept) / length(Y))/(credibility + 1)
        credibility <<- credibility + 1
      }
    }
    sim
  }
}
