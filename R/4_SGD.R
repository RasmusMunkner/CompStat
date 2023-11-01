

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
    optimizer = Vanilla_Optimizer(),
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
           "momentum"= Momentum_Optimizer(...),
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
      index_permutation <- sample(
        optimizable$n_index, optimizable$n_index, replace = F)
    } else {
      index_permutation <- 1:optimizable$n_index
    }

    # Apply minibatch gradient updates
    for (b in 1:batches_per_epoch){

      par_now <- par_next

      grad <- optimizable$grad(
        par_now,
        index_permutation[(1+(b-1)*batch_size):
                            min(1+b*batch_size, optimizable$n_index)]
      )

      update <- opt$lr(epoch) * opt$update_param(grad)

      par_next <- par_now - update

      if (trace_precision == "batch"){
        trace <- trace %>% extend(par_next)
        if (stop_crit$check(
          epoch,
          param = trace %>% tail(1),
          param_old = trace %>% tail(2),
          obj = trace %>% tail(1, type = "o"),
          obj_old = trace %>% tail(2, type = "o")
          )
        ){
          return(trace)
        }
      }

    }

    if (trace_precision == "epoch"){
      trace <- trace %>% extend(par_next)
      if (stop_crit$check(
        epoch,
        param = trace %>% tail(1),
        param_old = trace %>% tail(2),
        obj = trace %>% tail(1, type = "o"),
        obj_old = trace %>% tail(2, type = "o")
        )
      ){
        return(trace)
      }
    }

  }

  trace <- trace %>% extend(par_next)
  return(trace)
}

#' Wrapper function for cpp sgd for the logistic loglikelihood
#'
#' @param design
#' @param init_coef
#' @param y
#' @param lr
#' @param beta1
#' @param beta2
#' @param eps
#' @param batch_size
#' @param stop_crit
#' @param disable_adam
#'
#' @return
#' @export
#'
#' @examples
#' # Demonstration of the functionality
#' n <- 1000
#' p <- 3
#' sll <- simple_logistic_loglikelihood(n, p)
#' init_coef <- runif(p)
#' lr <-  1e-2
#' batch_size <- 100
#' maxiter <- 1000
#' trace <- SGC_CPP_Wrapper(
#' design = sll$X,
#' init_coef = init_coef,
#' y = sll$y,
#' lr = lr,
#' stop_crit = maxiter,
#' batch_size = batch_size,
#' disable_adam = F
#' )
#'
#' trace %>% dplyr::mutate(iter = dplyr::row_number()) %>%
#' tidyr::pivot_longer(cols = -iter, names_to = "coef", values_to = "value") %>%
#' ggplot2::ggplot(ggplot2::aes(x = iter, y = value, color = coef)) +
#' ggplot2::geom_line()
SGC_CPP_Wrapper <- function(
    design, init_coef, y,
    pen_matrix = matrix(0, nrow = ncol(design), ncol = ncol(design)),
    lambda = 0.001,
    lr = 1e-3, beta1 = 0.9, beta2 = 0.95, eps = 1e-8, batch_size = 32,
    stop_crit = 50,
    disable_adam = T
    ){

  stop_crit <- stopping_criterion(stop_crit)
  if (!is.function(lr)){
    if(length(lr) != stop_crit$maxiter){
      if (length(lr) == 1){
        lrate <- rep(lr, stop_crit$maxiter)
      } else {
        stop("lr should be a scalar, a vector of length maxiter or a function.")
      }
    } else {
      lrate <- lr
    }
  } else {
    lrate <- lr(1:stop_crit$maxiter)
  }


  if (disable_adam){
    beta1 <- 0
    beta2 <- 1
    eps <- 1
  }

  coef_trace <- SGD_CPP(
    design = design,
    coef = init_coef,
    y = y,
    lr = lrate,
    maxiter = stop_crit$maxiter,
    batch_size = batch_size,
    adam_beta1 = beta1,
    adam_beta2 = beta2,
    adam_eps = eps
    ) %>%
    purrr::map_dfr(.f = function(x) x %>% as.vector() %>% setNames(paste0("p", seq_along(init_coef))))

  return(coef_trace)

}



















