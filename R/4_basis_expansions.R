
#' Template for constructing CompStatBasisExpansions
#'
#' @param X The model matrix
#' @param Omega The penalty matrix
#' @param expand_new A function which computes the expansion for a new observation
#'
#' @return A CompStatBasisExpansion
#' @export
#'
#' @examples
#' one_hot <- CompStatBasisExpansion(
#' diagmat(rep(1,9)),
#' 1:9,
#' diagmat(rep(1,9)),
#' function(x){
#' purrr::map(.x = x, .f = function(x){c(rep(0, x-1), 1, rep(0, 9-x))}) %>%
#' purrr::reduce(.f = rbind) %>% magrittr::set_rownames(NULL)
#' }
#' )
#' one_hot$expand_new(c(1,1,1,2,3,8))
CompStatBasisExpansion <- function(X, x, Omega, expand_new){
  structure(list(
    X = X,
    x = x,
    Omega = Omega,
    expand_new = expand_new
  ), class = "CompStatBasisExpansion")
}

#' Plotting method for CompStatBasisExpansions
#'
#' @param BasisExpansion A CompStatBasisExpansion
#' @param beta A matching vector of coefficients
#'
#' @return A ggplot of the basis expansion. The x-values used for creating the
#' the basis expansion are used for x-axis in the plot
#' @export
#'
#' @examples
#' # One-hot encoding example
#' one_hot <- CompStatBasisExpansion(
#' diagmat(rep(1,9)),
#' 1:9,
#' diagmat(rep(1,9)),
#' function(x){
#' purrr::map(.x = x, .f = function(x){c(rep(0, x-1), 1, rep(0, 9-x))}) %>%
#' purrr::reduce(.f = rbind) %>% magrittr::set_rownames(NULL)
#' }
#' )
#' plot(one_hot, rnorm(9))
plot.CompStatBasisExpansion <- function(BasisExpansion, beta, post_transform = function(x) x){
  y <- post_transform(BasisExpansion$X %*% beta)
  data.frame(x = BasisExpansion$x, y = y) %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "x", y = "f(x)")
}

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
rescale_covariate <- function(x, mean_shift = NULL, sd_shift = NULL, method = "minmax"){
  if (!is.null(mean_shift) & !is.null(sd_shift)){
    return((x - mean_shift) / sd_shift)
  } else {
    switch(method,
           minmax = (x-min(x)) / (max(x) - min(x)),
           standard = (x - mean(x)) / sd(x),
           quantile = stop("Not implemented quantile"),
           default = x
    )
  }
}

#' Build a B-spline CompStatBasisExpansion
#'
#' @param x A vector for which the B-splines should be evaluated over
#' @param n_knots The number of knots
#'
#' @return
#' @export
#'
#' @examples
#' bspline_basis <- ExpandBspline(horses$Temperature)
#' plot(bspline_basis, rnorm(bspline_basis$X %>% ncol()))
ExpandBspline <- function(x, n_knots = 8){

  x_min <- min(x)
  x_max <- max(x)
  mean_shift <- x_min
  sd_shift <- x_max - x_min

  x_minmax <- x %>%
    rescale_covariate(mean_shift, sd_shift)

  knots <- x_minmax %>%
    quantile(seq(0,1,length.out=n_knots))

  X <- x_minmax %>%
    splines::splineDesign(knots = knots %>% extend_and_sort())

  Omega <- spline_pen_mat(knots)

  expand_new <- function(x){
    x %>%
      rescale_covariate(mean_shift, sd_shift) %>%
      splines::splineDesign(knots = knots %>% extend_and_sort())
  }

  CompStatBasisExpansion(X, x, Omega, expand_new)

}

#' Convenience function for sorting a vector and adding three copies of each endpoint to the vector
#'
#' @param x A vector.
#'
#' @details
#' Why do we need this? Its useful for construction of B-spline bases. Not sure why.
#'
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
#' The code is borrowed from chp.2 in the CompStat book. It relies on the splines being of order 3.
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

#' Build a Polynomial CompStatBasisExpansion
#'
#' @param x The data within which the basis expansion should be performed
#' @param degree The largest degree of polynomial used
#'
#' @return A CompStatBasisExpansion
#' @export
#'
#' @examples
#' PolyBasis <- ExpandPoly(rnorm(500), degree = 8)
#' plot(PolyBasis, rnorm(ncol(PolyBasis$X), sd = 20))
ExpandPoly <- function(x, degree = 4){
  x_min <- min(x)
  x_max <- max(x)
  mean_shift <- x_min
  sd_shift <- x_max - x_min

  x_minmax <- x %>%
    rescale_covariate(mean_shift, sd_shift)

  X <- x_minmax %>%
    vandermonde(n = degree) %>% t()

  Omega <- outer(0:degree, 0:degree, FUN = function(i,j){
    ifelse(i > 1 & j > 1, i*(i-1)*j*(j-1)/(i + j - 3), 0)
  })

  expand_new <- function(x){
    vandermonde(x, n = degree) %>% t()
  }

  CompStatBasisExpansion(X, x, Omega, expand_new)

}

