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
  arma::vec dg = - (y / p - (1 - y) / (1 - p)) / N;
  arma::vec grad = dp.t() * dg + 2 * lambda * pen_matrix * coef;
  return(grad);
}

//[[Rcpp::export]]
double lllC(
    const arma::mat &design, const arma::vec &coef, const arma::vec &y,
    const arma::mat &pen_matrix, const double &lambda
){
  double N = y.n_elem;
  arma::vec eta = exp(design * coef);
  arma::vec p = eta / (1 + eta);
  arma::mat g = - arma::accu(y % log(p) + (1 - y) % log(1-p)) / N + lambda * coef.t() * pen_matrix * coef;
  return(arma::accu(g)); // Call to accu is just for type conversion
}

// The follow function is a full cpp implementation of SGD for the logistic regression
// [[Rcpp::export]]
Rcpp::List SGD_CPP_PRIMITIVE(
    const arma::mat &design, arma::vec coef, const arma::vec &y,
    const arma::mat &pen_matrix, const double &lambda,
    const NumericVector &lr, const IntegerVector &batch_size,
    const int &maxiter, const double &objtarget,
    const double &beta_1, const double &beta_2, const double &eps,
    const bool &amsgrad,
    const int &seed
){
  Rcpp::List coef_list (maxiter+1);
  Rcpp::List obj_list (maxiter+1);
  double obj;
  coef_list[0] = as<NumericVector>(wrap(coef)); // Breaks mutability of coef
  obj = lllC(
    design,
    coef,
    y,
    pen_matrix,
    lambda
  );
  obj_list[0] = 1 * obj;
  arma::vec grad;
  arma::vec rho = vec(coef.n_elem); rho.zeros();
  arma::vec nu = vec(coef.n_elem); nu.zeros();
  int batch_size_now;
  int batches_per_epoch;
  uvec indicies_all, indicies;
  dqrng::dqRNGkind("Xoroshiro128+");
  dqrng::dqset_seed(IntegerVector::create(seed)); // Sets the shuffling seed

  for (int i = 0; i < maxiter; i++){
    indicies_all = as<uvec>(dqrng::dqsample_int(y.n_elem, y.n_elem));

    batch_size_now = std::min<int>(y.n_elem, batch_size[i]);
    batches_per_epoch = std::ceil(double(y.n_elem) / batch_size_now);

    for (int b = 0; b < batches_per_epoch; b++){
      indicies = indicies_all.subvec(
        b * batch_size_now,
        std::min<int>((b+1) * batch_size_now - 1, y.n_elem - 1)
      );

      grad = lll_gradC(
        design.rows(indicies),
        coef,
        y.elem(indicies),
        pen_matrix,
        lambda
        );

      rho = rho * beta_1 + grad * (1 - beta_1);
      if (amsgrad){
        nu = arma::max(
          nu * beta_2 + pow(grad, 2) * (1 - beta_2), nu
        );
      } else {
        nu = nu * beta_2 + pow(grad, 2) * (1 - beta_2);
      }

      coef = coef - lr[i] * rho / (sqrt(nu) + eps);
    }

    coef_list[i+1] = as<NumericVector>(wrap(coef)); // Breaks mutability of coef
    obj = lllC(
      design,
      coef,
      y,
      pen_matrix,
      lambda
    );
    obj_list[i+1] = 1 * obj;

    if (obj < objtarget){
      break;
    }
  }

  Rcpp::List results (2);
  results[0] = coef_list;
  results[1] = obj_list;

  return(results);
}






