
#' Rescale covariates
#'
#' @param x Covariate vector to be rescaled
#' @param method A string specifying the rescaling method
#'
#' @return
#' @export
#'
#' @examples
#' rescale_covariate(horses$Temperature, method = "minmax")
#' rescale_covariate(horses$Temperature, method = "standard")
rescale_covariate <- function(x, method = "minmax"){
  switch(method,
         minmax = (x-min(x)) / (max(x) - min(x)),
         standard = (x - mean(x)) / sd(x),
         quantile = stop("Not implemented quantile"),
         default = x
         )
}

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
stopping_criterion <- function(tol_obj = NULL,
                               norm_obj = function(x) sum(x),
                               tol_param = NULL,
                               norm_param = function(x) sqrt(sum(x^2)),
                               maxiter = 50
                               ){
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
#' @param K A number of epochs
#' @param p The power of the decay
#'
#' @return An object of class CompStatDecaySchedule.
#' @export
#'
#' @examples
#' poly_decay <- polynomial_schedule(1, 0.01)
polynomial_schedule <- function(lr_start, lr_later, K = 100, p = 1){
  lr_schedule <- function(epoch){
    lr_start / (1 + epoch^p / K)
  }
  structure(list(
    lr = lr_schedule
    ),
    class = "CompStatDecaySchedule"
  )
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
  structure(list(
    lr = lr_schedule
  ),
  class = "CompStatDecaySchedule"
  )
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
parameter_trace <- function(trace){
  structure(list(
    trace = trace,
    final = trace[[length(trace)]]
  ),
  class = "CompStatParameterTrace"
  )
}

#' Plotting method for CompStatParameterTrace
#'
#' @param x An object of class 'CompStatParameterTrace'
#' @param ... Additional arguments passed to ggplot
#'
#' @return A ggplot of the parameter traces
#' @export
plot.CompStatParameterTrace <- function(x, ...){
  x$trace %>%
    purrr::reduce(.f = rbind) %>%
    as.data.frame() %>%
    dplyr::mutate(Epoch = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = -Epoch, names_to = "Parameter", values_to = "Value") %>%
    ggplot2::ggplot(ggplot2::aes(x = Epoch, y = Value, color = Parameter), ...) +
    ggplot2::geom_line()
}

#' Print method for CompStatParameterTrace
#'
#' @param x A CompStatParameter Trace
#' @param ... Additional arguments passed to print.default
#'
#' @return Prints a summary of the parameter trace
#' @export
#'
#' @examples
print.CompStatParameterTrace <- function(x, ...){
  print("---CompStatParameterTrace---", quote=F)
  print("Final Parameters:", quote=F)
  print(x$final, ...)
}

#' Generic wrapper around CompStatOptimizable creation
#'
#' @param objective The objective function to be minimized
#' @param grad The gradient of the objective function
#' @param n_param The dimensionality of the parameter input for the objective
#'
#' @return A 'CompStatOptimizable' object
#' @export
CompStatOptimizable <- function(objective, grad, n_param, n_index){
  structure(list(
    objective = objective,
    grad = grad,
    n_param = n_param,
    n_index = n_index
  ), class = "CompStatOptimizable")
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





