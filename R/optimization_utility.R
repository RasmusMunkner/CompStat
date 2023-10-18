
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





