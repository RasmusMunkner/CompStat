#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double epanechnikov_l2norm_runningC(NumericVector x, double r){
  std::sort(x.begin(), x.end());
  int n = x.size();
  x.push_back(R_PosInf);
  double results = 0;
  int i = 1;
  int j = 0;
  while(j < n-1){
    if(i <= j || x[i+1] - x[j] < 2*r){
      ++i;
    } else if (x[i] - x[j] < 2*r) {
      for(int k = j+1; k <= i; ++k){
        results -= x[k];
      }
      results = results + (2*r + x[j]) * (i - j);
      ++j;
    } else {
      ++j;
    }
  }
  results = 2*results + 2*r*n;
  return(
    results * 9 / 4 / pow(n,2) / pow(r, 6)
  );
}
