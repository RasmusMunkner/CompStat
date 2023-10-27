

#' Gradient descent method for CompStatOptimizable
#'
#' @param optimizable A CompStatOptimizable
#' @param optimizer A CompStatOptimizer or a key corresponding to one
#' @param init_param Initial parameters for optimization.
#' @param lr Learning rate schedule
#' @param stop_crit A CompStatStoppingCriterion. Alternatively a number of epochs.
#' @param ... Additional arguments passed to the optimizer if it was specified via a key
#'
#' @return
#' @export
#'
#' @examples
#' parabola_optim <- optimizable_parabola(c(0,1,-2))
#' trace <- SGD(parabola_optim, lr = 0.33)
#' trace %>% plot()
#' trace %>% plot(type = "o")
SGD <- function(
    optimizable,
    optimizer = "vanilla",
    init_param = NULL,
    stop_crit = 50,
    shuffle = T,
    batch_size = 1,
    ...
    ){

  browser()

  # Ensure stopping criterion is valid
  if (!(class(stop_crit) %in% c("CompStatStoppingCriterion"))){
    stop_crit <- stopping_criterion(maxiter = stop_crit)
  }

  # Determine optimizer
  if (!(class(optimizer) %in% c("CompStatOptimizer"))){
    opt <-
      switch(optimizer,
           "vanilla"= Vanilla_Optimizer(...),
           "adam" = Adam_Optimizer(...)
           )
  } else {
    opt <- optimizer
  }

  # Useful control quantities
  batches_per_epoch <- ceiling(optimizable$n_index / batch_size)
  max_total_updates <- batches_per_epoch * stop_crit$maxiter

  # Initialize parameters
  trace <- CompStatTrace(optimizable)
  if (is.null(init_param)){
    init_param <- rep(NA, optimizable$n_param)
  }
  init_param <- init_param %>% dplyr::coalesce(rnorm(optimizable$n_param))
  trace <- extend(trace, init_param)

  for (epoch in 1:stop_crit$maxiter){

    # Reshuffle observations
    if (shuffle){
      index_permutation <- sample(optimizable$n_index, optimizable$n_index, replace = F)
    } else {
      index_permutation <- 1:optimizable$n_index
    }

    # Apply minibatch gradient updates
    for (b in 1:batches_per_epoch){
      timestep <- epoch * batches_per_epoch + b

      grad <- optimizable$grad(
        trace %>% tail(),
        index_permutation[1+(b-1)*batch_size, min(1+b*batch_size, optimizable$n_index)]
      )

      update <- opt$lr(epoch) * opt$update_param(grad)

      trace <- trace %>% extend(tail(trace) - update)

    }

    if (stop_crit$check(epoch,
                        param = trace %>% tail(1), param_old = trace %>% tail(2),
                        obj = trace %>% tail(1, type = "o"), trace %>% tail(2, type = "o"))
        ){
      return(trace)
    }

  }

  return(trace)
}




















