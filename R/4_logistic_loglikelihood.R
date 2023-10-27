

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
    -mean(response[batch] * log(p) + (1-response[batch]) * log(1-p)) + 2 * lambda * t(coef) %*% penalty_matrix %*% coef
  }

  grad <- function(coef, batch = 1:nrow(design)){
    eta <- exp(design[batch,] %*% coef)
    -t(design[batch,]) %*% (eta / (1 + eta) + response[batch] * (1-eta) / (1+eta)) / length(batch) + 2 * lambda * penalty_matrix %*% coef
  }

  structure(list(
    objective = loglikelihood,
    grad = grad,
    n_param = nrow(penalty_matrix),
    n_index = length(response)
  ),
  class = "CompStatOptimizable"
  )
}

#' Convenience function for sorting a vector and adding three copies of each endpoint to the vector
#'
#' @param x A vector.
#'
#' @return The same vector, but sorted and with 3 copies of the endpoints added to either end.
#' @export
#'
#' @examples
#' 1:6 %>% extend_and_sort()
extend_and_sort <- function(x){
  sort(c(rep(range(x), 3), x))
}

#' Calculate the penalty matrix using spline-design for basis splines
#'
#' The code is borrowed from chp.2 in the CompStat book
#'
#' @param inner_knots The inner knots used for the basis
#'
#' @return The penalty matrix corresponding to the calculated basis.
#' @export
#'
#' @examples
#' horses$Temperature %>%
#' rescale_covariate() %>%
#'   quantile(seq(0.2, 0.8, 0.1)) %>%
#'   spline_pen_mat()
spline_pen_mat <- function(inner_knots) {
  knots <- extend_and_sort(inner_knots)
  d <- diff(inner_knots)  # The vector of knot differences; b - a
  g_ab <- splines::splineDesign(knots, inner_knots, derivs = 2)
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splines::splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) +
      4 * crossprod(d * g_ab_mid, g_ab_mid) +
      crossprod(d * g_b, g_b)) / 6
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
