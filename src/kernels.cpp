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
    double tmp;

    for(int i = 0; i < m; ++i){
      for(int j = 0; j < n; ++j){
        tmp = pow((grid[i] - x[j])/bw,2);
        if (tmp < 1){
          out[i] += 1 - tmp;
        }
      }
      out[i] /= bw * n * 4 / 3;
    }

  }

  return out;

}
