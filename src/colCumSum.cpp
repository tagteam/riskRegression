// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//' Apply cumsum in each column 
//'
//' @description Fast computation of apply(x,2,cumsum)
//' @param x A matrix.
//' @return A matrix of same size as x.
//' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
//' @examples
//' x <- matrix(1:8,ncol=2)
//' colCumSum(x)
//' @export
// [[Rcpp::export]]
NumericMatrix colCumSum(NumericMatrix x){
  arma::mat m(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat result;
  result=cumsum(m,0);
  return wrap(result);
}
