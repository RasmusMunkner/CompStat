

#' Slow rejection sampling for arbitrary envelope
#'
#' This method is abysmally slow and was only implemented to provide a baseline.
#'
#' @param n The number of simulations.
#' @param enve An object of class 'Envelope'.
#' @param alpha
#'
#' @return Simulations from the distribution that the envelope is based on.
#' @export
#'
#' @examples
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
#' rejection_sampler_naive(10, enve)
rejection_sampler_naive <- function(n, enve){
  alpha <- 1/enve$c
  y_vec <- vector("numeric")
  n_sim <- 1
  while (n_sim <= n){
    U <- runif(1)
    Y <- enve$sim(1)
    if (U <= alpha * enve$f(Y) / enve$g(Y)){
      y_vec[n_sim] <- Y
      n_sim <- n_sim + 1
    }
  }
  y_vec
}

#' Compile rejection samplers for arbitrary envelope
#'
#' @param enve An object of class 'Envelope'.
#' @param alpha A scalar value ensuring \eqn{\alpha f \leq g}.
#'
#' @return A rejection sampler for the given envelope.
#' @export
#'
#' @examples
#' #Normal Distribution
#' enve <- LogLinearEnvelope(dnorm, function(z){dnorm(z) * (-z)}, c(-2,0,1))
#' sampler <- rejection_sampler_factory(enve)
#' sim <- sampler(50000)
#' hist(sim, prob=TRUE)
#' curve(enve$f(x), -3, 3, col="blue", add=TRUE)
#' qqnorm(sim)
#' qqline(sim)
#' shapiro.test(sim[1:5000])
#' environment(sampler)$p
#' environment(sampler)$credibility
#'
#' #Gamma Distribution
#' enve <- LogLinearEnvelope(
#'   function(z){dgamma(z,2)},
#'   function(z){dgamma(z,2) / z - dgamma(z, 2)},
#'   c(0.5, 2, 4.5)
#' )
#' sampler <- rejection_sampler_factory(enve)
#' sim <- sampler(50000)
#' hist(sim, prob=TRUE)
#' curve(enve$f(x), -3, 3, col="blue", add=TRUE)
rejection_sampler_factory <- function(enve, alpha = NULL){
  if (is.null(alpha) & "LogLinearEnvelope" %in% class(enve)){
    alpha <- 1/enve$c
  } else if (is.null(alpha)) {
    stop("For a general envelope, you must specify an appropriate alpha.")
  }
  p <- alpha
  credibility <- 20 # Arbitrary
  rejection_sampler <- function(n, train = TRUE){
    sim <- vector("numeric", n)
    n_sim <- 0
    while (n_sim < n){

      # Calculate RV's
      U <- runif(ceiling((n - n_sim)/p))
      Y <- enve$sim(ceiling((n - n_sim)/p))
      Y_accept <- Y[U <= alpha * enve$f(Y) / enve$g(Y)]

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


