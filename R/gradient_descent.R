

#' Gradient descent method for CompStatOptimizable
#'
#' @param optimizable A CompStatOptimizable
#' @param init_param Initial parameters for optimization.
#' @param lr Learning rate schedule
#' @param stop_crit A CompStatStoppingCriterion. Alternativly a number of epochs.
#'
#' @return
#' @export
#'
#' @examples
#' parabola_optim <- optimizable_parabola(c(0,1,-2))
#' GD(parabola_optim)
GD <- function(
    optimizable,
    init_param = NULL,
    lr = 1e-3,
    stop_crit = 50
    ){
  browser()

  # Ensure stopping criterion is valid
  if (!(class(stop_crit) %in% c("CompStatStoppingCriterion"))){
    stop_crit <- stopping_criterion(maxiter = stop_crit)
  }

  # Ensure the learning rate is callable
  if (!(class(lr) %in% c("CompStatDecaySchedule"))){
    lr <- constant_schedule(lr)
  }

  # Initialize parameters
  param <- vector(mode = "list", length = attr(stop_crit, "maxiter") + 1)
  if (!is.null(init_param)){
    param[[1]] <- init_param
  } else {
    param[[1]] <- rnorm(optimizable$n_param)
  }

  for (epoch in 2:length(param)){
    param[epoch] <- param[[epoch-1]] - lr$lr(epoch) * optimizable$grad(param[[epoch-1]])
    if (stop_crit(epoch, param = param[[epoch-1]], old_param = param[[epoch]])){
      return(param[1:epoch])
    }
  }

  return(param)

}

#' Test case for optimization algorithms
#'
#' @param minima A series of target parameter values
#'
#' @return A CompStatOptimizable where the minima are the optimal solution
#' @export
#'
#' @examples
#' parabola_optim <- optimizable_parabola(c(0,1,-2))
optimizable_parabola <- function(minima){
  f <- function(coef){
    sum((coef - minima)^2)
  }
  grad <- function(coef){
    sum(2*coef)
  }
  structure(list(
    objective = f,
    grad = grad,
    n_param = length(minima)
  ),
  class = "CompStatOptimizable"
  )
}




















