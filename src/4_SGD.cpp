#include <Rcpp.h>
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
arma::vec batch_gradient(arma::mat design, arma::vec coef, arma::vec y){
  arma::vec eta = design * coef;
  arma::vec zeta = eta/(1+eta) + y * (1-eta) / (1+eta);


}
