#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

// The following function is supposed to calculate the gradient for the logistic likelihood
// The inputs are the design matrix, the model coefficients and the responses
// [[Rcpp::export]]
arma::vec batch_gradient(arma::mat design, arma::vec coef, arma::vec y, arma::mat pen_matrix, double lambda){
  double N = y.n_elem;
  arma::vec eta = exp(design * coef);
  arma::vec p = eta / (1 + eta);
  arma::mat dp = arma::diagmat(eta / pow((1 + eta), 2)) * design;
  //double g = - arma::accu(y % log(p) + (1 - y) % log(1-p)) / N;
  arma::vec dg = - (y / p - (1 - y) / (1 - p)) / N;
  arma::vec grad = dp.t() * dg + 2 * lambda * coef.t() * pen_matrix;

  //Rcpp::List result;
  //result["obj"] = g;
  //result["grad"] = grad;

  return(grad);
}

// The follow function is a full cpp implementation of SGD for the logistic regression
// [[Rcpp::export]]
Rcpp::List SGD_CPP(
    arma::mat design, arma::vec coef, arma::vec y,
    arma::mat pen_matrix, double lambda,
    NumericVector lr, int maxiter, int batch_size,
    double adam_beta1, double adam_beta2, double adam_eps
){
  NumericVector objs (maxiter);
  Rcpp::List grads (maxiter);
  Rcpp::List coef_list (maxiter);
  //Rcpp::List h_eval;
  arma::vec grad;
  arma::vec rho = vec(coef.n_elem); rho.zeros();
  arma::vec nu = vec(coef.n_elem); nu.zeros();
  batch_size = std::min<int>(y.n_elem, batch_size);
  int batches_per_epoch = y.n_elem / batch_size; // Integer division rounds down, its intentional
  uvec indicies_all, indicies;

  for (int i = 0; i < maxiter; i++){
    indicies_all = randperm(y.n_elem);

    for (int b = 0; b < batches_per_epoch; b++){
      indicies = indicies_all.subvec(b * batch_size, (b+1) * batch_size - 1);
      grad = batch_gradient(
        design.rows(indicies),
        coef,
        y.elem(indicies),
        pen_matrix,
        lambda
        );
      rho = rho * adam_beta1 + grad * (1 - adam_beta1);
      nu = nu * adam_beta2 + pow(grad, 2) * (1 - adam_beta2);
      coef = coef - lr[i] * rho / (sqrt(nu) + adam_eps);
    }

    coef_list[i] = as<NumericVector>(wrap(coef)); // Type conversion to break mutability
    //objs[i] = as<double>(h_eval[0]);
    //grads[i] = as<NumericVector>(h_eval[1]);
  }

  // Rcpp::List result;
  // result["obj"] = objs;
  // result["grad"] = grads;
  //result["coef"] = coef_list;

  return(coef_list);
}









