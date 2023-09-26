
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
#' enve <- LogLinearEnvelope(get_rv(), c(-2,0,1))
#' sampler <- rejection_sampler_factory(enve)
#' new_sampler <- rejection_sampler_factory(enve)
#' sim <- sampler(50000)
#' plot(enve)
#' sampler(5000000, adapt_enve = T)
#' sampler %>% environment %>% `[[`("enve") %>% plot()
#' hist(sim, prob=TRUE)
#' curve(enve$base_rv$f(x), -3, 3, col="blue", add=TRUE)
#' qqnorm(sim)
#' qqline(sim)
#' shapiro.test(sim[1:5000])
#' environment(sampler)$p
#' environment(sampler)$credibility
#'
#' #Gamma Distribution
#' enve <- LogLinearEnvelope(get_rv("g"))
#' sampler <- rejection_sampler_factory(enve)
#' sim <- sampler(50000)
#' hist(sim, prob=TRUE, ylim = c(0,0.5))
#' curve(enve$base_rv$f(x), 0, 8, col="blue", add=TRUE)
#'
#' # Laplacian Envelope
#' enve <- LaplaceEnvelope(dnorm, sim_method = 2)
#' sampler <- rejection_sampler_factory(enve)
#' sim <- sampler(50000)
#' hist(sim, prob=TRUE, ylim=c(0,0.6))
#' curve(enve$base_rv$f(x), -3, 3, col="blue", add=TRUE)
rejection_sampler_factory <- function(enve, evalmode = 2){
  p <- min(1 - 1 / (2 + enve$alpha), 0.9) # An initial guess at the rejection probability
  credibility <- 1#20 # Arbitrary
  rejection_sampler <- function(n, env = NULL, train = TRUE, adapt_enve = FALSE){

    browser()

    if (is.null(env)){
      env <- enve
    }
    sim <- vector("numeric", n)
    n_sim <- 0
    while (n_sim < n){

      # Calculate RV's
      Y <- env$sim(ceiling((n - n_sim)/p))
      if (0 %in% evalmode){
        U <- runif(ceiling((n - n_sim)/p))
        x1 <- env$alpha
        x2 <- env$base_rv$f(Y)
        x3 <- env$f(Y)
        filter <- U <= x1 * x2 / x3
      }
      if (1 %in% evalmode) {
        U <- runif(ceiling((n - n_sim)/p))
        x1 <- log(env$alpha)
        x2 <- env$base_rv$log_f(Y)
        x3 <- env$log_f(Y)
        filter <- log(U) <= x1 + x2 - x3
      }
      if (2 %in% evalmode){
        E <- rexp(ceiling((n - n_sim)/p))
        x1 <- log(env$alpha)
        x2 <- env$base_rv$log_f(Y)
        x3 <- env$log_f(Y)
        filter <- -E <= x1 + x2 - x3
      }
      Y_accept <- Y[filter]

      if ("LogLinearEnvelope" %in% class(enve) & adapt_enve & sum(!filter) > 15){
        Y_reject_dens <- Y[!filter] %>% sample(min(1e5, length(.))) %>% density(bw="SJ")
        update_point <- sample(Y_reject_dens$x, size = 1, prob = (Y_reject_dens$y)^3 / sum((Y_reject_dens$y)^3))
        enve <<- update(enve, tangent_points = c(enve$tangent_points, update_point))
        env <- enve
      }


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


