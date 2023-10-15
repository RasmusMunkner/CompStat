
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
rejection_sampler <- function(enve, evalmode = 2){
  p <- min(1 - 1 / (2 + enve$alpha), 0.9) # An initial guess at the rejection probability
  credibility <- 1#20 # Arbitrary
  sampler <- function(n, env = NULL, train = TRUE, adapt_enve = FALSE){

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
  sampler
}

#' Basic benchmark of different envelope types
#'
#' @param N_seq The implementations are tested for 2^n for each n in this sequence.
#' @param rv A RandomVariable describing the target distribution
#'
#' @return The benchmark data
#' @export
#'
#' @examples
#' bm <- benchmark_envelopes()
#' bm %>%
#' ggplot2::ggplot(ggplot2::aes(x = log(N), y = log(Time), color = Method)) +
#' ggplot2::geom_line() +
#' ggplot2::geom_point()
benchmark_envelopes <- function(N_seq = 4:12, rv = named_rv("epanechnikov")){
  set.seed(0)
  gauss_enve <- GaussianEnvelope(rv, security = 1)
  laplace_enve <- LaplaceEnvelope(rv, security = 1)
  loglinear_enve <- LogLinearEnvelope(rv, autoselection_msg = F, precompute = F)
  adapted_enve <- LogLinearEnvelope(rv, autoselection_msg = F, precompute = F)

  gauss_sampler <- rejection_sampler(gauss_enve, evalmode = 0)
  laplace_sampler <- rejection_sampler(laplace_enve, evalmode = 0)
  loglinear_sampler <- rejection_sampler(loglinear_enve, evalmode = 0)
  adapted_sampler <- rejection_sampler(adapted_enve, evalmode = 0)

  for (i in 1:50){
    adapted_sampler(10000, train = T, adapt_enve = T)
  }

  bm <- N_seq %>%
    purrr::imap_dfr(.f = function(N, i){
      calls <- list(
        call("gauss_sampler", n = 2^N),
        call("laplace_sampler", n = 2^N),
        call("loglinear_sampler", n = 2^N, train = F),
        call("adapted_sampler", n = 2^N, train = F)
      )
      names(calls) <- c("Gaussian", "Laplace", "LogLinear", "Adapted LogLinear")
      microbenchmark::microbenchmark(
        list=calls
      ) %>%
        dplyr::group_by(expr) %>%
        dplyr::summarise(Time = median(time)) %>%
        dplyr::mutate(N = 2^N, Method = expr) %>%
        dplyr::select(-c(expr))
    })

  set.seed(NULL)
  bm

}


