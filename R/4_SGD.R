

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
#' negloglik <- logistic_loglikelihood(response = horses$dead, design = horses$Temperature %>% ExpandBspline())
#' trace_simple <- SGD(negloglik, Vanilla_Optimizer(0.1), batch_size = 100, tracing = F, stop_crit = 150, seed = 0)
#' trace_epoch <- SGD(negloglik, Vanilla_Optimizer(0.1), stop_crit = 1e2, seed = 0)
#' trace_epoch %>% plot()
SGD <- function(
    optimizable,
    optimizer = Vanilla_Optimizer(),
    init_param = NULL,
    stop_crit = 50,
    shuffle = T,
    batch_size = 1,
    tracing = T,
    objtarget = -Inf,
    seed = NULL,
    ...
    ){

  dqrng::dqRNGkind("Xoroshiro128+")
  dqrng::dqset.seed(seed)

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
  opt$reset()

  # Useful control quantities
  batches_per_epoch <- ceiling(optimizable$n_index / batch_size)
  max_total_updates <- batches_per_epoch * stop_crit$maxiter

  # Initialize parameters

  if (is.null(init_param)){
    init_param <- rep(0, optimizable$n_param)
  }
  init_param <- init_param %>% dplyr::coalesce(0)
  par_next <- init_param

  trace <- matrix(
    NA, nrow = (stop_crit$maxiter+1), ncol = optimizable$n_param + 1)
  trace[1,] <- c(par_next, optimizable$objective(par_next))

  for (epoch in 1:stop_crit$maxiter){

    # Reshuffle observations
    if (shuffle){
      index_permutation <- dqrng::dqsample.int(
        optimizable$n_index, optimizable$n_index, replace = F)
    } else {
      index_permutation <- 1:optimizable$n_index
    }

    par_before <- par_next
    # Apply minibatch gradient updates
    for (b in 1:batches_per_epoch){

      par_now <- par_next

      grad <- optimizable$grad(
        par_now,
        index_permutation[(1+(b-1)*batch_size):
                            min(b*batch_size, optimizable$n_index)]
      )

      update <- opt$lr(epoch) * opt$update_param(grad)

      par_next <- par_now - update
    }

    if (tracing == T){
      trace[epoch+1,] <- c(par_next, optimizable$objective(par_next))
    }

    if (stop_crit$check(
      epoch,
      param = par_next,
      param_old = par_before,
      obj = optimizable$objective(par_next),
      obj_old = optimizable$objective(par_before)
    )
    ){
      return(
        if (tracing){
          trace %>%
            magrittr::set_colnames(c(paste0("p", 1:optimizable$n_param), "obj")) %>%
            magrittr::set_class(c("CompStatTrace", class(.))) %>%
            return()
        } else {
          c(par_next, optimizable$objective(par_next)) %>%
            matrix(nrow = 1) %>%
            magrittr::set_colnames(c(paste0("p", 1:optimizable$n_param), "obj")) %>%
            return()
        }

      )
    }
  }

  stop("An error with maxiter occured. SGD should call return from within the iteration loop.")
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
SGD_CPP <- function(
    lll, init_coef = NULL,
    lr = 1e-3, beta1 = 0.9, beta2 = 0.95, eps = 1e-8, batch_size = 32,
    stop_crit = 50, objtarget = 0,
    disable_adam = F, amsgrad = T,
    seed = NULL
    ){

  if (!("CompStatLogisticLogLikelihood" %in% class(lll))){
    stop(paste0("Input lll must be of class 'CompStatLogisticLogLikelihood'."))
  } else {
    design <- lll$design
    pen_matrix <- lll$penalty_matrix
    lambda <- lll$lambda
    y <- lll$response
  }

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

  if (is.null(init_coef)){
    init_coef <- rnorm(ncol(design))
  }


  if (disable_adam){
    beta1 <- 0
    beta2 <- 1
    eps <- 1
  }

  trace <- SGD_CPP_PRIMITIVE(
    design = design,
    coef = init_coef,
    y = y,
    lr = lrate,
    pen_matrix = pen_matrix,
    lambda = lambda,
    maxiter = stop_crit$maxiter,
    batch_size = batch_size,
    adam_beta1 = beta1,
    adam_beta2 = beta2,
    adam_eps = eps,
    amsgrad = amsgrad,
    seed = seed,
    objtarget = objtarget
    )

  trace[[1]] %>%
    purrr::map_dfr(.f = function(x) x %>% as.vector() %>% setNames(paste0("p", seq_along(init_coef)))) %>%
    cbind(data.frame(obj = trace[[2]] %>% unlist(use.names = F)))

}



















