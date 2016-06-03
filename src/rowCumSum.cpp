#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
NumericMatrix rowCumSum(NumericMatrix m){
  for (int i = 0; i < m.nrow(); ++i) {
    for (int j = 1; j < m.ncol(); ++j) {
      m(i, j) += m(i, j-1);
    }
  }
  return m;
}
