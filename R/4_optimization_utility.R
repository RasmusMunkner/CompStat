

#' Stopping criterions for optimization algorithms
#'
#' @param tol_obj The tolerance level for objective value change.
#' @param norm_obj The norm used to evaluate the change in objective values.
#' @param tol_param The tolerance level for parameter value change.
#' @param norm_param The norm used to evaluate the change in parameter values.
#' @param maxiter The maximal number of iterations allowed.
#'
#' @return True if the algorithm should halt, false if it should continue.
#' @export
#'
#' @examples
#' sc <- stopping_criterion(maxiter = 75, tol_param = 0.01)
#' sc(50)
#' sc(100)
stopping_criterion <- function(
    maxiter = 50,
    tol_obj = NULL,
    threshhold_obj = -Inf,
    norm_obj = function(x) sum(x),
    tol_param = NULL,
    norm_param = function(x) sqrt(sum(x^2))
   ){

  # If stopping_criterion is called on a CompStatStoppingCriterion, do nothing
  if (class(maxiter) %in% c("CompStatStoppingCriterion")){
    return(maxiter)
  }

  stopper <- function(
    epoch,
    param = NULL, param_old = NULL,
    obj = NULL, obj_old = NULL){

    # Check maxiter criterion
    if (maxiter <= epoch){
      return(TRUE)
    }

    # Check criteria on parameter
    if (!is.null(param) & !is.null(param_old)){
      if (!is.null(tol_param)){
        if (norm_param(param - param_old) <= tol_param * (norm_param(param_old) + tol_param)){
          return(TRUE)
        }
      }
    }

    # Check criterion on objective
    if (!is.null(obj) & !is.null(obj_old)){
      if (!is.null(tol_obj)){
        if (norm_param(obj_old - obj) <= tol_obj * (norm_param(obj_old) + tol_obj)){
          return(TRUE)
        }
      }
    }

    # Check if objective is small enough
    if (!is.null(obj)){
      if (obj < threshhold_obj){
        return(TRUE)
      }
    }

    return(FALSE)
  }

  structure(list(
    check = stopper,
    maxiter = maxiter,
    tol_obj = tol_obj,
    tol_param = tol_param,
    threshhold_obj = threshhold_obj
    ),
    class = "CompStatStoppingCriterion"
  )
}

#' Polynomial learning rate schedule
#'
#' @param lr_start The initial learning rate
#' @param lr_later The learning rate after K epochs
#' @param later A number of epochs
#' @param p The power of the decay
#'
#' @return An object of class CompStatDecaySchedule.
#' @export
#'
#' @examples
#' poly_decay <- polynomial_schedule(1, 0.01)
polynomial_schedule <- function(lr_start, lr_later, later = 100, p = 1){
  K <- later^p * lr_later / (lr_start - lr_later)
  lr_schedule <- function(epoch){
    lr_start / (1 + epoch^p / K)
  }
  lr_schedule
}

#' Constant learning rate schedule
#'
#' @param lr A constant used as a learning rate
#'
#' @return An object of type CompStatDecaySchedule
#' @export
#'
#' @examples
#' const_decay <- constant_schedule(1)
constant_schedule <- function(lr){
  lr_schedule <- function(epoch) lr
  lr_schedule
}

#' Plotting function for CompStatTrace's
#'
#' @param trace A matrix with columns p1-pp and obj.
#'
#' @return
#' @export
plot.CompStatTrace <- function(trace, what = "po"){
  plotdata <-
  trace %>%
    dplyr::mutate(Iteration = dplyr::row_number())
  p1 <- plotdata %>%
    dplyr::select(-obj) %>%
    tidyr::pivot_longer(cols = -Iteration, names_to = "Parameter", values_to = "Value") %>%
    ggplot2::ggplot(ggplot2::aes(x = Iteration, y = Value, color = Parameter)) +
    ggplot2::geom_line() +
    ggplot2::ggtitle("Parameters")
  p2 <- plotdata %>%
    ggplot2::ggplot(ggplot2::aes(x = Iteration, y = obj)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iteration", y = "Objective") +
    ggplot2::ggtitle("Objective Function")
  switch(what,
         "p" = p1,
         "o" = p2,
         gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)
         )
}

CompareTraces <- function(traces, what = "po"){
  plotdata <- traces %>%
    purrr::imap_dfr(.f = function(trace, name){
      trace %>%
        dplyr::mutate(Iteration = dplyr::row_number(),
                      Algorithm = name)
    })
  p1 <- plotdata %>%
    dplyr::select(-obj) %>%
    tidyr::pivot_longer(cols = -c(Iteration, Algorithm), names_to = "Parameter", values_to = "Value") %>%
    ggplot2::ggplot(ggplot2::aes(x = Iteration, y = Value, color = Parameter)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~Algorithm) +
    ggplot2::ggtitle("Parameters")
  p2 <- plotdata %>%
    ggplot2::ggplot(ggplot2::aes(x = Iteration, y = obj, color = Algorithm)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iteration", y = "Objective") +
    ggplot2::ggtitle("Objective Function")
  switch(what,
         "p" = p1,
         "o" = p2,
         gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)
  )
}

#' Extract the last row from a CompStatTrace
#'
#' @param trace A CompStatTrace
#'
#' @return
#' @export
#'
#' @examples
tail.CompStatTrace <- function(trace, type = "p"){
  switch(type,
         "p" = trace[nrow(trace), 1:(ncol(trace)-1)],
         "o" = trace[nrow(trace), ncol(trace)],
         trace[nrow(trace), ]
         ) %>% unlist()
}

#' Generic wrapper around CompStatOptimizable creation
#'
#' @param objective The objective function to be minimized
#' @param grad The gradient of the objective function
#' @param n_param The dimensionality of the parameter input for the objective
#'
#' The objective function and gradient functions should take two arguments,
#' - A numeric vector of parameters of length n_param
#' - A numeric vector of indicies of any length
#'
#' @return A 'CompStatOptimizable' object
#' @export
CompStatOptimizable <- function(objective, grad, n_param, n_index){
  structure(list(
    objective = objective,
    grad = grad,
    n_param = n_param,
    n_index = n_index
  ),
  class = "CompStatOptimizable")
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
  f <- function(param, index = NULL){
    sum((param - minima)^2)
  }
  grad <- function(param, index = NULL){
    2 * (param - minima)
  }
  n_param <- length(minima)
  CompStatOptimizable(f, grad, n_param, n_index = 1)
}

#' Create a logistic regression setup with known coefficients
#'
#' @param n The number of observations
#' @param p The number of parameters
#' @param beta The coefficients. Defaults to 1:p
#'
#' @return A list containing the regression problem and the true coefficients
#' @export
#'
#' @examples
#' sll <- simple_logistic_loglikelihood(n = 1000)
#' data <- cbind(sll$X, sll$y) %>% magrittr::set_colnames(c(paste0("x", 1:ncol(sll$X)), "y")) %>% data.frame()
#' glmfit <- glm(y ~ .-1, family = binomial(), data = data)
#' glmfit$coef
simple_logistic_loglikelihood <- function(n = 10, p = 2, beta = 1:p){
  X <- matrix(rnorm(n * p), ncol = p)
  y <- rbinom(n, 1, p = exp(X %*% beta)/(1+exp(X %*% beta)))
  return(list(X = X, y = y, beta = beta))
}

#' Create a regression with known targets
#'
#' @param n The number of observations
#' @param p The number of coefficients
#'
#' @return A list containing a model matrix, the response and the true coefficients
#' @export
#'
#' @examples
#' #Showing an example of SGD optimization for this loglikelihood
#' sll <- simple_logistic_loglikelihood_optimizable(n = 1000)
#' trace <- SGD(sll, init_param = c(0.8, 2.1), lr = 0.05, batch_size = 100, stop_crit = stopping_criterion(maxiter = 50))
#' trace %>% plot(type = "o")
#' trace %>% plot(type = "p")
simple_logistic_loglikelihood_optimizable <- function(n = 10, p = 2, beta = 1:p){
  sll <- simple_logistic_likelihood(n,p,beta)
  logistic_loglikelihood(
    design = sll$X,
    response = sll$y,
    penalty_matrix = matrix(0, nrow = ncol(sll$X), ncol = ncol(sll$X)),
    lambda = 0
  )
}

#' Builds a univariate logistic regression dataset along with the
#' regression function.
#'
#' @param n The number of observations
#'
#' @return A list containing a covariate x, a binary response y and a function
#' f_true(x) = E[Y|X=x].
#' @export
#'
#' @examples
#' reg <- univariate_logistic_regression(10)
univariate_logistic_regression <- function(
    n = 1000, depth = 6, width = 10, mean = 0, sd = 2/3
    ){

  x <- runif(n, -1, 1)

  NN <- function(x, layers){
    x <- t(x)
    for (i in seq_along(layers)){
      linmap <- layers[[i]] %*% x
      if (i == length(layers)){
        x <-  exp(linmap) / (1+exp(linmap))
      } else {
        x <- tanh(linmap)
      }
    }
    x
  }

  make_nn_layers <- function(
    depth, width, mean, sd
    ){
    1:depth %>% purrr::map(.f = function(d){
      col <- ifelse(d == 1, 1, width)
      row <- ifelse(d == depth, 1, width)
      matrix(rnorm(row * col, mean, sd), nrow = row, ncol = col)
    })
  }

  layers <- make_nn_layers(depth, width, mean, sd)
  p <- NN(x, layers)
  attributes(p) <- NULL
  y <- rbinom(n, 1, p)

  return(structure(list(
    x = x,
    y = y,
    f_true = function(x) {NN(x,layers) %>% t()}
  ), class = "CompStatRegression"))
}

#' Plotting function for CompStatRegression
#'
#' @param reg A CompStatRegression
#'
#' @export
#'
#' @examples
#' reg <- univariate_logistic_regression(100, sd = 2/3, depth = 6)
#' plot(reg)
plot.CompStatRegression <- function(reg){
  x <- seq(-1,1,0.001)
  y <- reg$f_true(x)
  data.frame(x = x, y = y) %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line() +
    ggplot2::lims(x = c(-1,1)) +
    ggplot2::labs(x = "x", y = "p(y=1|X=x)")
}



#' Build a random trace
#'
#' @param optimizable A CompStatOptimizable to build the trace for
#' @param maxiter A number of maximal iterations
#'
#' @return
#' @export
#'
#' @examples
#' parabola <- optimizable_parabola(c(-3,2,1,0))
#' trace <- just_screw_around(parabola)
just_screw_around <- function(optimizable, maxiter = 50){
  trace <- CompStatTrace(optimizable)
  for (epoch in 1:maxiter){
    trace <- trace %>% extend(rnorm(optimizable$n_param))
  }
  trace
}




