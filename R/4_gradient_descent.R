

#' Gradient descent method for CompStatOptimizable
#'
#' @param optimizable A CompStatOptimizable
#' @param init_param Initial parameters for optimization.
#' @param lr Learning rate schedule
#' @param stop_crit A CompStatStoppingCriterion. Alternatively a number of epochs.
#'
#' @return
#' @export
#'
#' @examples
#' parabola_optim <- optimizable_parabola(c(0,1,-2))
#' SGD(parabola_optim, lrate = 0.33) %>% plot()
SGD <- function(
    optimizable,
    optimizer = "vanilla",
    init_param = NULL,
    lr = 1e-3,
    stop_crit = 50,
    shuffle = T,
    batch_size = 1,
    ...
    ){

  # Ensure stopping criterion is valid
  if (!(class(stop_crit) %in% c("CompStatStoppingCriterion"))){
    stop_crit <- stopping_criterion(maxiter = stop_crit)
  }

  # Determine optimizer
  if (!(class(optimizer) %in% c("CompStatOptimizer"))){
    opt <-
      switch(optimizer,
           "vanilla"= Vanilla_Optimizer(lr),
           "adam" = Adam_Optimizer(lr, ...)
           )
  } else {
    opt <- optimizer
  }

  # Initialize parameters
  param <- vector(mode = "list", length = (stop_crit$maxiter + 1) * ceiling(optimizable$n_index / batch_size))
  if (!is.null(init_param)){
    param[[1]] <- init_param
  } else {
    param[[1]] <- rnorm(optimizable$n_param)
  }

  for (epoch in 2:length(param)){

    # Reshuffle observations
    if (shuffle){
      index_permutation <- sample(optimizable$n_index, optimizable$n_index, replace = F)
    } else {
      index_permutation <- 1:optimizable$n_index
    }

    # Apply minibatch gradient updates
    for (b in 1:ceiling(optimizable$n_index / batch_size)){

      grad <- optimizable$grad(
        param[[epoch-1]],
        index_permutation[1+(b-1)*batch_size, min(1+b*batch_size, optimizable$n_index)]
      )

      update <- opt$lr(epoch) * opt$update_param(grad)

      param[[epoch]] <- param[[epoch-1]] - update

      if (stop_crit$check(epoch, param = param[[epoch-1]], param_old = param[[epoch]])){
        return(param[1:epoch] %>% parameter_trace())
      }

    }

  }

  return(param %>% parameter_trace())
}




















