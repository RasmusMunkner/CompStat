

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
    norm_obj = function(x) sum(x),
    tol_param = NULL,
    norm_param = function(x) sqrt(sum(x^2))
   ){

  # If stopping_criterion is called on a CompStatStoppingCriterion, do nothing
  if (class(maxiter) %in% c("CompStatStoppingCriterion")){
    return(maxiter)
  }

  stopper <- function(epoch, param = NULL, param_old = NULL, obj = NULL, obj_old = NULL){

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
    return(FALSE)
  }

  structure(list(
    check = stopper,
    maxiter = maxiter
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

#' S3 object interface for tracing parameter optimization
#'
#' @param trace A list of parameter values obtained during the fitting procedure
#'
#' @return An object of type CompStatParameterTrace
#' @export
#'
#' @examples
#' parabola_optim <- optimizable_parabola(c(0,1,-2))
#' trace <- GD(parabola_optim, lrate = 0.33)
CompStatTrace <- function(optimizable){
  if (!("CompStatOptimizable" %in% class(optimizable))){
    stop("Can only build CompStateTrace for CompStatOptimizable.")
  }
  structure(list(
    parameter_trace = data.frame(matrix(vector(), 0, optimizable$n_param)) %>%
      magrittr::set_colnames(paste0("p", 1:optimizable$n_param)),
    objective_trace = numeric(),
    grad_trace = data.frame(matrix(vector(), 0, optimizable$n_param)) %>%
      magrittr::set_colnames(paste0("g", 1:optimizable$n_param)),
    optimizable = optimizable
  ), class = "CompStatTrace")
}

#' Generic method for extending objects
#'
#' @param extendable An object implementing an extension method
#' @param extension A valid extension for the extendable
#'
#' @return The extended version of the extendable
#' @export
extend <- function(extendable, extension){
  UseMethod("extend")
}

#' Convenience function to log new data to trace
#'
#' @param trace A CompStatTrace
#' @param param A new vector of parameters to add to the trace
#'
#' @return An updates CompStatTrace
#' @export
#'
#' @examples
#' parabola <- optimizable_parabola(c(-2,0,1))
#' trace <- CompStatTrace(parabola)
#' trace <- trace %>% extend(c(3,2,1)) %>% extend(c(0,2,1)) %>% extend(c(-2,0,1))
extend.CompStatTrace <- function(trace, param){
  if (length(param) != trace$optimizable$n_param){
    stop(paste0("New parameter values must be of the correct length", trace$optimizable$n_param, ")"))
  }
  pnames <- colnames(trace$parameter_trace)
  gnames <- colnames(trace$grad_trace)
  trace$parameter_trace <- trace$parameter_trace %>%
    rbind(param, make.row.names = F) %>%
    magrittr::set_colnames(pnames)
  trace$objective_trace <- trace$objective_trace %>%
    append(param %>% trace$optimizable$objective())
  trace$grad_trace <- trace$grad_trace %>%
    rbind(param %>% trace$optimizable$grad(), make.row.names = F) %>%
    magrittr::set_colnames(gnames)
  return(trace)
}

#' Plotting method for CompStatParameterTrace
#'
#' @param x An object of class 'CompStatParameterTrace'
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot of the parameter traces
#' @export
#'
#' @examples
#' parabola_optim <- optimizable_parabola(c(0,1,-2))
#' trace <- just_screw_around(parabola_optim)
#' trace %>% plot(type = "p")
#'
plot.CompStatTrace <- function(trace, type = "p"){
  if (type == "p"){
    trace$parameter_trace %>%
      dplyr::mutate(Iteration = dplyr::row_number()) %>%
      tidyr::pivot_longer(cols = -Iteration, names_to = "Parameter", values_to = "Value", names_ptypes = factor()) %>%
      ggplot2::ggplot(ggplot2::aes(x = Iteration, y = Value, color = Parameter)) +
      ggplot2::geom_line()
  } else if (type == "o") {
    trace$objective_trace %>%
      tibble::tibble() %>%
      magrittr::set_colnames("Objective") %>%
      dplyr::mutate(Iteration = dplyr::row_number()) %>%
      tidyr::pivot_longer(cols = -Iteration, names_to = "Parameter", values_to = "Value") %>%
      ggplot2::ggplot(ggplot2::aes(x = Iteration, y = Value, color = Parameter)) +
      ggplot2::geom_line()
  }
}

#' Print method for CompStatParameterTrace
#'
#' @param trace A CompStatTrace
#' @param ... Additional arguments passed to print.default
#'
#' @return Prints a summary of the parameter trace
#' @export
#'
#' @examples
print.CompStatTrace <- function(trace, ...){
  print("---CompStatTrace---", quote=F)
  print("Final Parameters:", quote=F)
  print(trace$parameter_trace %>% tail(1), ...)
  print("Final Objective:", quote=F)
  print(trace$objective_trace %>% tail(1), ...)
}

#' Extract the most recent parameters for a CompStatTrace
#'
#' @param trace A CompStatTrace
#'
#' @return
#' @export
#'
#' @examples
#' parabola_optim <- optimizable_parabola(c(0,1,-2))
#' trace <- just_screw_around(parabola_optim, maxiter = 5)
#' trace %>% tail()
#' trace %>% tail(2)
#' trace %>% tail(type = "o")
tail.CompStatTrace <- function(trace, n = 1, type = "p"){
  if (type == "p"){
    trace$parameter_trace %>% `[`(nrow(.) - n + 1,) %>% unlist()
  } else if (type == "o"){
    trace$objective_trace %>% .subset2(length(.) - n + 1)
  }
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
  make_logistic_loglikelihood(
    design = sll$X,
    response = sll$y,
    penalty_matrix = matrix(0, nrow = ncol(sll$X), ncol = ncol(sll$X)),
    lambda = 0
  )
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




