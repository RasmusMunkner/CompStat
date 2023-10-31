

#' Build a CompStatOptimizable for the logistic loglikelihood
#'
#' @param design The design matrix.
#' @param response The response.
#' @param penalty_matrix The penalty matrix corresponding to the design.
#' @param lambda The strengt of the penalization.
#'
#' @return A CompStatOptimizable exposing the negative, penalized loglikelihood
#' and the gradient.
#' @export
#'
#' @examples
#' logistic_opfun <- make_logistic_loglikelihood()
#' trace <- just_screw_around(logistic_opfun)
#' trace <- SGD(logistic_opfun, lr = 0.0001, stop_crit = 20)
#' trace <- SGD(logistic_opfun, lr = 0.0001, batch_size = 32)
#' logistic_opfun$objective(rep(0,9))
#' logistic_opfun$grad(rep(0,9))
#' logistic_opfun$objective(rep(1,9), batch = 1:100)
#' logistic_opfun$grad(rep(1,9), batch = 1:100)
make_logistic_loglikelihood <- function(design = horses$Temperature %>%
                                       rescale_covariate() %>%
                                       splines::splineDesign(
                                         knots = quantile(., seq(0.2, 0.8, 0.1)) %>% extend_and_sort(),
                                         x = ., outer.ok = T),
                                     response = horses$dead,
                                     penalty_matrix = horses$Temperature %>%
                                       rescale_covariate() %>%
                                       quantile(seq(0.2, 0.8, 0.1)) %>%
                                       spline_pen_mat(),
                                     lambda = 0.001
){

  loglikelihood <- function(coef, batch = 1:nrow(design)){
    eta <- exp(design[batch,] %*% coef)
    p <- eta / (1 + eta)
    -mean(response[batch] * log(p) + (1-response[batch]) * log(1-p)) + lambda * t(coef) %*% penalty_matrix %*% coef
  }

  grad <- function(coef, batch = 1:nrow(design)){
    X <- design[batch,]
    y <- response[batch]

    eta <- exp(X %*% coef)
    p <- eta / (1 + eta)
    dp <- t(X) %*% diagmat(eta / (1 + eta)^2)
    dg <-  - (y/p - (1-y)/(1-p)) / length(batch)
    grad <- (dp %*% dg %>% as.vector()) + ((2 * lambda * penalty_matrix %*% coef) %>% as.vector())
    grad

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

# spline_pen_mat(c(1,2,3,4,5,6))
#
# horses$Temperature %>% rescale_covariate() %>%
#   cbind(poly(.,6) %>% as.data.frame()) %>%
#   `names<-`(c("Temperature", paste0("Degree", 1:6))) %>%
#   pivot_longer(cols = -Temperature, names_to="Basis", values_to="splineval") %>%
#   ggplot(aes(x = Temperature, y = splineval, color = Basis)) +
#   geom_line()
#
# horses$Temperature %>% rescale_covariate() %>% hist()
#
# horses$Temperature %>% rescale_covariate() %>% poly(10) %>% cov()
#
# htemp <- horses$Temperature %>% rescale_covariate()
# spline_htemp <- htemp %>% splines::splineDesign(knots = seq(0.1, 0.9, length.out = 8), x = ., outer.ok = T)
#
#
# tibble(htemp = htemp) %>%
#   cbind(spline_htemp %>% setNames(paste0("S", 1:4))) %>%
#   as_tibble() %>%
#   pivot_longer(cols = -c(htemp)) %>%
#   ggplot(aes(x = htemp, y = value, color = name)) +
#   geom_line()
#
#
