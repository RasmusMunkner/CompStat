

#' Gradient descent method for CompStatOptimizable
#'
#' @param optimizable A CompStatOptimizable
#' @param optimizer A CompStatOptimizer or a key corresponding to one
#' @param init_param Initial parameters for optimization.
#' @param lr Learning rate schedule
#' @param stop_crit A CompStatStoppingCriterion. Alternatively a number of epochs.
#' @param trace_precision One of 'batch', 'epoch'. Other values indicate no tracing.
#' @param ... Additional arguments passed to the optimizer if it was specified via a key
#'
#' @return
#' @export
#'
#' @examples
#' parabola_optim <- optimizable_parabola(c(0,1,-2))
#' negloglik <- logistic_opfun <- make_logistic_loglikelihood()
#' trace_batch <- SGD(negloglik, lr = 1e-4, batch_size = 100, trace_precision = "batch")
#' trace_epoch <- SGD(negloglik, lr = 1e-4, batch_size = 100, trace_precision = "epoch")
#' trace_batch %>% plot()
#' trace_epoch %>% plot()
SGD <- function(
    optimizable,
    optimizer = "vanilla",
    init_param = NULL,
    stop_crit = 50,
    shuffle = T,
    batch_size = 1,
    trace_precision = "batch",
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
  par_next <- init_param
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

      par_now <- par_next

      grad <- optimizable$grad(
        par_now,
        index_permutation[(1+(b-1)*batch_size):min(1+b*batch_size, optimizable$n_index)]
      )

      update <- opt$lr(epoch) * opt$update_param(grad)

      par_next <- par_now - update

      if (trace_precision == "batch"){
        trace <- trace %>% extend(par_next)
        if (stop_crit$check(epoch,
                            param = trace %>% tail(1), param_old = trace %>% tail(2),
                            obj = trace %>% tail(1, type = "o"), trace %>% tail(2, type = "o"))
        ){
          return(trace)
        }
      }

    }

    if (trace_precision == "epoch"){
      trace <- trace %>% extend(par_next)
      if (stop_crit$check(epoch,
                          param = trace %>% tail(1), param_old = trace %>% tail(2),
                          obj = trace %>% tail(1, type = "o"), trace %>% tail(2, type = "o"))
      ){
        return(trace)
      }
    }

  }

  trace <- trace %>% extend(par_next)

  return(trace)
}




















