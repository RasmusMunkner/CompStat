
#' Gradient descent using Wolfie condition and backtracing
#'
#' @param optimizable
#' @param init_par
#' @param stop_crit
#' @param d
#' @param gamma0
#' @param c
#' @param tracing
#' @param em_update_tol
#'
#' @return
#' @export
#'
#' @examples
#' See tests.
GD <- function(
    optimizable,
    init_par = NULL,
    stop_crit = 50,
    d = 0.98,
    gamma = 1,
    tracing = T,
    em_update_tol = 1,
    debug = F,
    ...
    ){

  if (debug){
    browser()
  }

  EM_flag <- FALSE
  if ("CompStatQfunc" %in% class(optimizable)){
    EM_flag <- TRUE
  }

  # Ensure stopping criterion is valid
  stop_crit <- stopping_criterion(stop_crit)

  # Initialize parameters
  if (is.null(init_par)){
    init_par <- rep(0, optimizable$n_param)
  }
  init_par <- init_par %>% dplyr::coalesce(0)
  if (EM_flag){
    init_par <- optimizable$get_w()
  }

  par_next <- init_par
  obj_next <- optimizable$objective(init_par)

  trace <- matrix(
    NA, nrow = (stop_crit$maxiter+1), ncol = optimizable$n_param + 1)
  trace[1,] <- c(par_next,
                 ifelse(EM_flag,
                        optimizable$loglikelihood(par_next),
                        obj_next)
    )

  for (epoch in 1:stop_crit$maxiter){

    # Remember the previous parameters
    par_before <- par_next
    obj_before <- obj_next

    grad <- optimizable$grad(par_before, ...) / optimizable$n_index
    while (gamma > 1e-9){
      par_next <- par_before - gamma * grad
      obj_next <- optimizable$objective(par_next)
      if (obj_next < obj_before){
        break
      } else {
        gamma <- gamma * d
      }
    }

    if (tracing == T){
      trace[epoch+1,] <- c(
        par_next,ifelse(EM_flag, optimizable$loglikelihood(par_next), obj_next))
    }

    if (EM_flag){
      if (obj_next <
          obj_before * ifelse(obj_before > 0, em_update_tol, 1/em_update_tol)){
        optimizable$set_w(par_next)
      }
    }

    # Check if stopping criterion is fulfilled
    if (stop_crit$check(
      epoch,
      param = par_next,
      param_old = par_before,
      obj = ifelse(EM_flag,
                   optimizable$loglikelihood(par_next), obj_next),
      obj_old = ifelse(EM_flag,
                       optimizable$loglikelihood(par_before), obj_before)
    )
    ){

      if (EM_flag){
        names <- c(names(par_next), "obj")
      } else {
        names <- c(paste0("p", 1:optimizable$n_param), "obj")
      }

      if (tracing){
        return(trace %>%
                 magrittr::set_colnames(names) %>%
                 as.data.frame() %>%
                 magrittr::set_class(c("CompStatTrace", class(.))))
      } else {
        return(c(par_next, obj_next) %>%
                 matrix(nrow = 1) %>%
                 magrittr::set_colnames(names) %>%
                 as.data.frame())
      }
    }

  }

  stop("An error with maxiter occured. GD should call return from within the iteration loop.")
}
