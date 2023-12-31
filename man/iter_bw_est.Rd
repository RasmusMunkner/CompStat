% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_density.R
\name{iter_bw_est}
\alias{iter_bw_est}
\title{Iterative bandwidth estimation}
\usage{
iter_bw_est(
  x,
  maxiter = 3L,
  kernel_code = "e",
  bw0 = NULL,
  tol = 1e-07,
  reltol = 0.001,
  ...
)
}
\arguments{
\item{x}{A numerical vector. Data to fit the density to}

\item{maxiter}{A positive integer. The maximal number of iterations}

\item{kernel_code}{One of "gaussian", "uniform", "triangular" or "epanechnikov".}

\item{bw0}{A scalar. Initial guess at the optimal bandwidth.
Defaults to maximum-likelihood Gaussian density}

\item{tol}{A scalar. Tolerance level bandwidth selection convergence}

\item{reltol}{A scalar. Relative tolerance for bandwidth selection convergence}

\item{...}{Additional arguments passed to estimate_l2norm.}
}
\value{
A list of bandwidths for estimation of the unknown density.
The final bandwidth represents the best-estimate for the optimal bandwidth.
}
\description{
Iterative bandwidth selection using AMISE.
The method runs until convergence of bandwidth selection and returns all
bandwidth estimates. To access density estimates, use compile_density.
}
\examples{
iter_bw_est(rnorm(1000), method = "matrix")
iter_bw_est(4*rbinom(1000, 1, 0.5) + rnorm(1000))
}
