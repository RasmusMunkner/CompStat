#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double getOne(){
  double n = 1;
  return exp(n);
}

// [[Rcpp::export]]
NumericVector eval_kdensC(String kcode, NumericVector grid, NumericVector x, double bw){

  int m = grid.size();
  int n = x.size();

  NumericVector out(m);

  if (kcode == "Gaussian"){

    for(int i = 0; i < m; ++i){
      for(int j = 0; j < n; ++j){
        out[i] += exp(-pow((grid[i] - x[j])/bw, 2)/2);
      }
      out[i] /= bw * n * sqrt(2*M_PI);
    }

  } else { //Assuming Epanechnikov Kernel

    for(int i = 0; i < m; ++i){
      for(int j = 0; j < n; ++j){
        NumericVector tmp(2);
        tmp[0] = 0;
        tmp[1] = 1 - pow((grid[i] - x[j])/bw,2);
        out[i] += max(tmp);
      }
      out[i] /= bw * n * 4 / 3;
    }

  }

  return out;

}
