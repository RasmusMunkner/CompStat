

#' Gradient descent method for CompStatOptimizable
#'
#' @param optimizable A CompStatOptimizable
#' @param optimizer A CompStatOptimizer or a key corresponding to one
#' @param init_par Initial parameters for optimization.
#' @param lr Learning rate schedule
#' @param stop_crit A CompStatStoppingCriterion. Alternatively a number of epochs.
#' @param trace_precision One of 'batch', 'epoch'. Other values indicate no tracing.
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
    init_par = NULL,
    stop_crit = 50,
    shuffle = T,
    tracing = T,
    seed = NULL
    ){

  # Enable compatibility with EM algorithm
  EM_flag <- ("CompStatQfunc" %in% class(optimizable))

  dqrng::dqRNGkind("Xoroshiro128+")
  dqrng::dqset.seed(seed)

  # Ensure stopping criterion is valid
  stop_crit <- stopping_criterion(stop_crit)

  # Determine optimizer
  if (!(class(optimizer) %in% c("CompStatOptimizer"))){
    opt <-
      switch(optimizer,
           "vanilla"= Vanilla_Optimizer(),
           "momentum"= Momentum_Optimizer(),
           "adam" = Adam_Optimizer()
           )
  } else {
    opt <- optimizer
  }
  opt$reset()

  # Initialize parameters
  if (is.null(init_par)){
    init_par <- rep(0, optimizable$n_param)
  }
  init_par <- init_par %>% dplyr::coalesce(0)

  # Reasonable defaults for EM algorithm
  # Note they are random to ensure algorithm is not stuck
  if (EM_flag){
    valid_flags <- optimizable$check_par_validity(init_par)
    if (!valid_flags$p){
      where_p <- optimizable$par_alloc == "p"
      init_par[where_p] <- 1 / (1+sum(where_p)) # Uniform
    }
    if (!valid_flags$sigma2){
      where_sigma2 <- optimizable$par_alloc == "sigma2"
      init_par[where_sigma2] <- 1 + rexp(sum(where_sigma2), 1)
    }
    if (!valid_flags$nu){
      where_nu <- optimizable$par_alloc == "nu"
      init_par[where_nu] <- 1 + rexp(sum(where_nu), 1)
    }
    optimizable$set_w(init_par)
  }

  par_next <- init_par
  obj_next <- optimizable$objective(init_par)

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

    # Remember the previous parameters
    par_before <- par_next
    obj_before <- obj_next

    # Determine batch size
    batch_size <- optimizer$batch_size(epoch, obj = obj_before, n = optimizable$n_index)
    batches_per_epoch <- ceiling(optimizable$n_index / batch_size)

    # Apply minibatch gradient updates
    for (b in 1:batches_per_epoch){

      grad <- optimizable$grad(
        par_before,
        index_permutation[(1+(b-1)*batch_size):
                            min(b*batch_size, optimizable$n_index)]
      )

      update <- opt$lr(epoch) * opt$update_param(grad)

      if (EM_flag){ # Checks that parameter values are within the allowed limits
        tmp_par_next <- par_before - update
        for (attempt in 1:20){
          check <- optimizable$check_par_validity(tmp_par_next) %>% unlist()
          if (all(check)){
            par_next <- tmp_par_next
            break
          }
          tmp_par_next <- (tmp_par_next + par_before) / 2
        }
      } else {

        par_next <- par_before - update

      }
    }

    # Tracing and keep track of objective function
    obj_next <- optimizable$objective(par_next)
    if (tracing == T){
      trace[epoch+1,] <- c(par_next, obj_next)
    }

    # If the Q-function improved over the epoch, update the underlying par
    if (EM_flag){
      if (obj_next < obj_before){
        optimizable$set_w(par_next)
      }
    }

    # Check if stopping criterion is fulfilled
    if (stop_crit$check(
      epoch,
      param = par_next,
      param_old = par_before,
      obj = obj_next,
      obj_old = obj_before
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
#'
SGD_CPP <- function(
    lll,
    optimizer,
    init_par = NULL,
    stop_crit = 50,
    seed = NULL
    ){

  if (is.null(seed)){
    stop("Seed must be set for CPP implementation.")
  }

  # Check that the objective is a logistic loglikelihood
  if (!("CompStatLogisticLogLikelihood" %in% class(lll))){
    stop(paste0("Input lll must be of class 'CompStatLogisticLogLikelihood'."))
  } else {
    design <- lll$design
    pen_matrix <- lll$penalty_matrix
    lambda <- lll$lambda
    y <- lll$response
  }

  # Ensure the stopping criterion is a CompStatStoppingCriterion
  stop_crit <- stopping_criterion(stop_crit)

  # Call learning rate and batch size
  lr <- optimizer$lr(1:stop_crit$maxiter)
  batch_size <- optimizer$batch_size(1:stop_crit$maxiter, n = lll$n_index)

  # If nothing else is specified, initialize all parameters to 0
  if (is.null(init_par)){
    init_par <- rep(0,ncol(design))
  }

  # Call the cpp implementation
  trace <- SGD_CPP_PRIMITIVE(
    design = design,
    coef = init_par,
    y = y,
    lr = lr,
    pen_matrix = pen_matrix,
    lambda = lambda,
    maxiter = stop_crit$maxiter,
    batch_size = batch_size,
    beta_1 = optimizer$par$beta_1,
    beta_2 = optimizer$par$beta_2,
    eps = optimizer$par$eps,
    amsgrad = optimizer$par$amsgrad,
    seed = seed,
    objtarget = stop_crit$threshhold_obj
    )

  trace[[1]] %>%
    purrr::keep(.p = function(x) !is.null(x)) %>%
    purrr::map_dfr(.f = function(x) x %>% as.vector() %>% setNames(paste0("p", seq_along(init_par)))) %>%
    cbind(data.frame(obj = trace[[2]] %>% unlist(use.names = F))) %>%
    magrittr::set_class(c("CompStatTrace", class(.)))

}



















