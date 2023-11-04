#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

//[[Rcpp::export]]
arma::mat f_y_by_zC(
    const arma::vec &w, const arma::vec &y
){
  int K = w.n_elem / 4;
  int N = y.n_elem;
  arma::mat y_by_z (y.n_elem, K);

  for (int k = 0; k < K; k++){
    for (int n = 0; n < N; n++){
      y_by_z(n,k) = R::dt((y(n)-w(K+k))/ sqrt(w(2*K+k)), w(3*K+k), 0) / sqrt(w(2*K+k));
    }
  }
  return(y_by_z);
}

//[[Rcpp::export]]
arma::mat cond_pC(
    const arma::vec &w, const arma::vec &y
){
  uword K = w.n_elem / 4;
  uword N = y.n_elem;
  arma::mat fy_by_z (y.n_elem, K);

  for (uword k = 0; k < K; k++){
    for (uword n = 0; n < N; n++){
      fy_by_z(n,k) = R::dt((y(n)-w(K+k))/ sqrt(w(2*K+k)), w(3*K+k), 0) / sqrt(w(2*K+k));
    }
  }

  arma::vec fy = f_y_by_zC(w, y) * w.subvec(0, K - 1);

  return(fy_by_z % ((1/fy) * w.subvec(0, K - 1).t()));

}

//[[Rcpp::export]]
arma::vec dQC(arma::vec v, arma::vec w, arma::vec y){

  uword K = v.n_elem / 4;
  uword N = y.n_elem;

  arma::mat condp = cond_pC(w, y);

  arma::vec dBeta = arma::sum(condp).t() - N * w.subvec(0, K - 1);

  arma::mat mu_coef (N, K);
  for (uword k = 0; k < K; k++){
    for (uword n = 0; n < N; n++){
      mu_coef(n,k) = (w(3*K+k)+1)*(y(n)-w(K+k)) / (w(3*K+k)*w(2*K+k)+pow((y(n)-w(K+k)),2));
    }
  }
  arma::vec dMu = arma::sum(mu_coef % condp).t();

  arma::mat kappa_coef (N, K);
  for (uword k = 0; k < K; k++){
    for (uword n = 0; n < N; n++){
      kappa_coef(n,k) = -0.5 + 0.5 * (w(3*K+k)+1)*pow((y(n)-w(K+k)),2) / (w(3*K+k)*w(2*K+k)+pow((y(n)-w(K+k)),2));
    }
  }
  arma::vec dKappa = arma::sum(kappa_coef % condp).t();

  return(arma::join_cols(arma::join_cols(dBeta, dMu), dKappa));
}

//[[Rcpp::export]]
arma::mat dQC2(arma::vec v, arma::vec w, arma::vec y){

  uword K = v.n_elem / 4;
  uword N = y.n_elem;

  arma::mat condp = cond_pC(w, y); // This is where all the computational burden is

  arma::vec dBeta = arma::sum(condp).t() - N * w.subvec(0, K - 1);

  arma::mat mumat (N,K);
  arma::mat ymat (N,K);
  for (uword k = 0; k < K; k++){
    mumat.col(k).fill(w[K+k]);
    }
  for (uword n = 0; n < N; n++){
    ymat.row(n).fill(y[n]);
    }
  arma::mat y_minus_mu = ymat - mumat;
  //
  arma::mat nu_sigma (N,K);
  arma::mat nu_plus (N,K);
  for (uword k = 0; k < K; k++){
    nu_sigma.col(k).fill(w[3*K+k]*w[2*K+k]);
    nu_plus.col(k).fill(w[3*K+k]+1);
  }
  arma::mat mu_coef = nu_plus % y_minus_mu / (nu_sigma + pow(y_minus_mu, 2));
  arma::vec dMu = arma::sum(mu_coef % condp).t();
  //
  arma::mat kappa_coef = -0.5 + 0.5 * y_minus_mu % mu_coef;
  arma::vec dKappa = arma::sum(kappa_coef % condp).t();
  //
  return(arma::join_cols(arma::join_cols(dBeta, dMu), dKappa));
}
