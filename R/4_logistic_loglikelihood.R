

#' Build a CompStatOptimizable for the logistic loglikelihood
#'
#' @param response The response.
#' @param design The design matrix. Also accepts a CompStatBasisExpansion
#' @param penalty_matrix The penalty matrix corresponding to the design.
#' @param lambda The strengt of the penalization.
#'
#' @return A CompStatOptimizable exposing the negative, penalized loglikelihood
#' and the gradient.
#' @export
#'
#' @examples
#' logistic_opfun <- logistic_loglikelihood()
#' trace <- just_screw_around(logistic_opfun)
#' trace <- SGD(logistic_opfun, lr = 0.0001, stop_crit = 20)
#' trace <- SGD(logistic_opfun, lr = 0.0001, batch_size = 32)
#' logistic_opfun$objective(rep(0,9))
#' logistic_opfun$grad(rep(0,9))
#' logistic_opfun$objective(rep(1,9), batch = 1:100)
#' logistic_opfun$grad(rep(1,9), batch = 1:100)
logistic_loglikelihood <- function(
    response,
    design,
    penalty_matrix = NULL,
    lambda = 0.001
){

  if (c("CompStatBasisExpansion") %in% class(design)){
    penalty_matrix <- design$Omega
    design <- design$X
  } else if (is.null(penalty_matrix)){
    penalty_matrix <- diagmat(rep(1,ncol(design))) # Ordinary L2 Penalty
  }

  loglikelihood <- function(coef, batch = 1:nrow(design)){
    eta <- exp(design[batch,] %*% coef)
    p <- eta / (1 + eta)
    -mean(response[batch] * log(p) + (1-response[batch]) * log(1-p)) +
      lambda * t(coef) %*% penalty_matrix %*% coef
  }

  grad <- function(coef, batch = 1:nrow(design)){
    X <- design[batch,]
    y <- response[batch]

    eta <- exp(X %*% coef)
    p <- eta / (1 + eta)
    dp <- t(X) %*% diagmat(eta / (1 + eta)^2)
    dg <-  - (y/p - (1-y)/(1-p)) / length(batch)
    grad <- dp %*% dg %>% as.vector() + (2 * lambda * penalty_matrix %*% coef)
    grad %>% as.vector()
  }

  structure(list(
    objective = loglikelihood,
    grad = grad,
    n_param = nrow(penalty_matrix),
    n_index = length(response),
    design = design,
    response = response,
    penalty_matrix = penalty_matrix,
    lambda = lambda
  ),
  class = "CompStatOptimizable"
  )
}

#' Builds a diagonal matrix from a vector
#'
#' @param v A vector
#'
#' @return A diagonal matrix with v on the diagonal
#' @export
#'
#' @examples
#' diagmat(c(1,2,3))
diagmat <- function(v){
  m <- matrix(0, nrow = length(v), ncol = length(v))
  diag(m) <- v
  m
}
