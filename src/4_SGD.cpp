#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

// The following function is supposed to calculate the gradient for the logistic likelihood
// The inputs are the design matrix, the model coefficients and the responses
// [[Rcpp::export]]
Rcpp::List batch_gradient(arma::mat design, arma::vec coef, arma::vec y){
  double N = y.n_elem;
  arma::vec eta = design * coef;
  arma::vec p = exp(eta) / (1 + exp(eta));
  arma::mat dp = arma::diagmat(eta / pow((1 + eta), 2)) * design;
  double g = - arma::accu(y % log(p) + (1 - y) % log(1-p)) / N;
  arma::vec dg = - (y / p + (1 - y) / (1 - p)) / N;
  arma::vec grad = dp.t() * dg;

  Rcpp::List result;
  result["obj"] = g;
  result["grad"] = grad;

  return(result);
}

// The follow function is a full cpp implementation of SGD for the logistic regression
// [[Rcpp::export]]
Rcpp::List SGD_CPP(
    arma::mat design, arma::vec coef, arma::vec y,
    double lr, int maxiter, int batch_size
){
  NumericVector objs (maxiter);
  Rcpp::List grads (maxiter);
  Rcpp::List coef_list (maxiter);
  Rcpp::List h_eval;
  batch_size = std::min<int>(y.n_elem, batch_size);
  int batches_per_epoch = y.n_elem / batch_size; // Integer division rounds down, its intentional
  uvec indicies_all, indicies;

  for (int i = 0; i < maxiter; i++){
    indicies_all = randperm(y.n_elem);

    for (int b = 0; b < batches_per_epoch; b++){
      indicies = indicies_all.subvec(b * batch_size, (b+1) * batch_size - 1);
      h_eval = batch_gradient(
        design.rows(indicies),
        coef,
        y.elem(indicies)
        );
    }

    coef_list[i] = as<NumericVector>(wrap(coef));
    objs[i] = as<double>(h_eval[0]);
    grads[i] = as<NumericVector>(h_eval[1]);

    coef = coef - lr * as<arma::vec>(h_eval[1]);
  }

  Rcpp::List result;
  result["obj"] = objs;
  result["grad"] = grads;
  result["coef"] = coef_list;

  return(result);
}









