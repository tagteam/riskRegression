#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
NumericMatrix colCumSum(NumericMatrix m){
  for (int i = 1; i < m.nrow(); ++i) {
    for (int j = 0; j < m.ncol(); ++j) {
      m(i, j) += m(i-1, j);
    }
  }
  return m;
}
