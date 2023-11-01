#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
#include <dqrng.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace dqrng;

// The following function is supposed to calculate the gradient for the logistic likelihood
// The inputs are the design matrix, the model coefficients and the responses
// [[Rcpp::export]]
arma::vec lll_gradC(
    const arma::mat &design, const arma::vec &coef, const arma::vec &y,
    const arma::mat &pen_matrix, const double &lambda
){
  double N = y.n_elem;
  arma::vec eta = exp(design * coef);
  arma::vec p = eta / (1 + eta);
  arma::mat dp = arma::diagmat(eta / pow((1 + eta), 2)) * design;
  //double g = - arma::accu(y % log(p) + (1 - y) % log(1-p)) / N;
  arma::vec dg = - (y / p - (1 - y) / (1 - p)) / N;
  arma::vec grad = dp.t() * dg + 2 * lambda * pen_matrix * coef;
  return(grad);
}

// The follow function is a full cpp implementation of SGD for the logistic regression
// [[Rcpp::export]]
Rcpp::List SGD_CPP_PRIMITIVE(
    const arma::mat &design, arma::vec coef, const arma::vec &y,
    const arma::mat &pen_matrix, const double &lambda,
    const NumericVector &lr, const int &maxiter, int &batch_size,
    const double &adam_beta1, const double &adam_beta2, const double &adam_eps,
    const bool &amsgrad,
    const int &seed
){
  NumericVector objs (maxiter);
  Rcpp::List grads (maxiter);
  Rcpp::List coef_list (maxiter);
  arma::vec grad;
  arma::vec rho = vec(coef.n_elem); rho.zeros();
  arma::vec nu = vec(coef.n_elem); nu.zeros();
  batch_size = std::min<int>(y.n_elem, batch_size);
  int batches_per_epoch = std::ceil(double(y.n_elem) / batch_size);
  uvec indicies_all, indicies;
  dqrng::dqRNGkind("Xoroshiro128+");
  dqrng::dqset_seed(IntegerVector::create(seed)); // Sets the shuffling seed

  for (int i = 0; i < maxiter; i++){
    indicies_all = as<uvec>(dqrng::dqsample_int(y.n_elem, y.n_elem));

    for (int b = 0; b < batches_per_epoch; b++){
      indicies = indicies_all.subvec(
        b * batch_size,
        std::min<int>((b+1) * batch_size - 1, y.n_elem - 1)
      );

      grad = lll_gradC(
        design.rows(indicies),
        coef,
        y.elem(indicies),
        pen_matrix,
        lambda
        );

      rho = rho * adam_beta1 + grad * (1 - adam_beta1);
      if (amsgrad){
        nu = arma::max(
          nu * adam_beta2 + pow(grad, 2) * (1 - adam_beta2), nu
        );
      } else {
        nu = nu * adam_beta2 + pow(grad, 2) * (1 - adam_beta2);
      }

      coef = coef - lr[i] * rho / (sqrt(nu) + adam_eps);
    }

    coef_list[i] = as<NumericVector>(wrap(coef)); // Breaks mutability of coef
  }

  return(coef_list);
}








