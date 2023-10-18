
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

  attr(stopper, "maxiter") <- maxiter

  structure(
    stopper,
    class = "CompStatStoppingCriterion"
  )
}

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

constant_schedule <- function(lr){
  lr_schedule <- function(epoch) lr
  structure(list(
    lr = lr_schedule
  ),
  class = "CompStatDecaySchedule"
  )
}





