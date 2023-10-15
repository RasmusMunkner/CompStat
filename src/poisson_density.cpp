#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
NumericMatrix vandermondeC(NumericVector y, int n){
  NumericMatrix VMy(n + 1, y.size());

  for (int r = 0; r <= n; ++r){
    for (int c = 0; c < y.size(); ++c){
      VMy(r,c) = pow(y[c], c);
    }
  }

  return(VMy);

}
