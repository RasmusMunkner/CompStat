
#' Wrapper CompStatOptimizer construction
#'
#' @param lr Learning rate schedule. Should be a function or a scalar
#' @param update_param The update function. Applied to gradients before those
#' are used in SGD
#'
#' @return A CompStatOptimizer
#' @export
CompStatOptimizer <- function(lr, batch_size, update_param, reset, par){

  # Ensure the learning rate is callable
  if (is.function(lr)){
    lrate <- lr
  } else {
    lrate <- function(epoch, ...) {lr[pmin(epoch, length(lr))]}
  }
  # Ensure batch_size is callable
  if (is.function(batch_size)){
    bsize <- batch_size
  } else {
    bsize <- function(epoch, ...) {batch_size[pmin(epoch, length(batch_size))]}
  }

  structure(list(
    lr = lrate,
    batch_size = bsize,
    update_param = update_param,
    reset = reset,
    par = par
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
#' @return A CompStatOptimizer
#' @export
Adam_Optimizer <- function(
    lr = 1e-3, batch_size = 32, beta_1 = 0.95, beta_2 = 0.97, eps = 1e-8, amsgrad = T
    ){

  nu <- 0
  rho <- 0

  update_param <- function(grad){
    rho <<- beta_1 * rho + (1-beta_1) * grad
    nu_proposal <- beta_2 * nu + (1-beta_2) * grad^2
    if (amsgrad){
      nu <<- pmax(nu_proposal, nu)
    } else {
      nu <<- nu_proposal
    }
    rho / (sqrt(nu) + eps)
  }

  reset <- function(alpha_rho = 0.9, alpha_nu = 0.9, partial = F){
    if (partial){
      rho <<- rho * alpha_rho
      nu <<- nu * alpha_nu
    } else {
      rho <<- 0
      nu <<- 0
    }

  }

  CompStatOptimizer(
    lr, batch_size, update_param, reset,
    par = list(beta_1 = beta_1, beta_2 = beta_2, eps = eps, amsgrad = amsgrad)
    )
}

#' Momentum optimizer constructor
#'
#' @inheritParams Adam_Optimizer
#'
#' @return A CompStatOptimizer
#' @export
Momentum_Optimizer <- function(lr = 1e-3, batch_size = 32, beta_1 = 0.95){
  Adam_Optimizer(lr, batch_size, beta_1 = beta_1, beta_2 = 1, eps = 1)
}

#' Vanilla/Identity optimizer
#'
#' @inheritParams Adam_Optimizer
#'
#' @return A CompStatOptimizer
#' @export
Vanilla_Optimizer <- function(lr = 1e-3, batch_size = 32){
  Adam_Optimizer(lr, batch_size, beta_1 = 0, beta_2 = 1, eps = 1)
}













