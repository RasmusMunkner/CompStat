
CompStatOptimizer <- function(key, lr, update_param){

  # Ensure the learning rate is callable
  if (!is.function(lr)){
    lrate <- function(epoch) {lr}
  } else {
    lrate <- lr
  }

  structure(list(
    key = key,
    lr = lrate,
    update_param = update_param
  ),
  class = "CompStatOptimizer")
}

#' Adam Optimizer Class
#'
#' @param lr Learning rate
#' @param beta_1 Momentum parameter
#' @param beta_2 Second-moment parameter
#' @param eps Stability parameter. Only used to avoid division by 0
#'
#' @return A CompStatOptimizer with key 'Adam'
#' @export
Adam_Optimizer <- function(lr = 1e-3, beta_1 = 0.98, beta_2 = 0.999, eps = 1e-8, ...){

  rho <- 0
  nu <- 0

  update_param <- function(grad){
    rho <<- beta_1 * rho + (1-beta_1) * grad
    nu <<- beta_2 * nu + (1-beta_2) * grad^2
    rho / (sqrt(nu) + eps)
  }

  CompStatOptimizer("Adam", lr, update_param)

}

#' Some funny optimizer with more quickly vanishing momentum
#'
#' @param lr
#' @param beta_1
#' @param beta_2
#' @param vanish
#' @param eps
#' @param ...
#'
#' @return
#' @export
Vanishing_Adam_Optimizer <- function(lr = 1e-3, beta_1 = 0.98, beta_2 = 0.999, vanish = 1e-1, eps = 1e-8, ...){

  rho <- 0
  nu <- 0

  update_param <- function(grad){
    rho <<- (grad^2 / (grad^2 + vanish)) * beta_1 * rho + (1-beta_1) * grad
    nu <<- beta_2 * nu + (1-beta_2) * grad^2
    rho / (sqrt(nu) + eps)
  }

  CompStatOptimizer("Vanish", lr, update_param)

}

#' Vanilla/Identity optimizer
#'
#' @param lr A CompStatDecaySchedule or a scalar
#'
#' @return A CompStatOptimizer object applying pure gradient upgrades
#' @export
Vanilla_Optimizer <- function(lr = 1e-3, ...){

  update_param <- function(grad){
    grad
  }

  CompStatOptimizer("Vanilla", lr, update_param)
}













