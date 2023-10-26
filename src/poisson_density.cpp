#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

//' @importFrom Rcpp evalCpp
//' @useDynLib CompStat


// [[Rcpp::export]]
NumericMatrix vmC(NumericVector y, int n){
  NumericMatrix VMy(n + 1, y.size());

  for (int r = 0; r <= n; ++r){
    for (int c = 0; c < y.size(); ++c){
      VMy(r,c) = pow(y[c], r);
    }
  }

  return(VMy);

}

// [[Rcpp::export]]
NumericMatrix vmC2(NumericVector y, int n){
  NumericMatrix VMy(n + 1, y.size());

  for (int c = 0; c < y.size(); ++c){
    VMy(0,c) = 1;
  }

  for (int c = 0; c < y.size(); ++c){
    for (int r = 1; r <= n; ++r){
      VMy(r,c) = VMy(r-1,c) * y[c];
    }
  }

  return(VMy);

}
