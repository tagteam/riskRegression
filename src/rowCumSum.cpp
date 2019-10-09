// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Apply cumsum in each row 
//'
//' @description Fast computation of t(apply(x,1,cumsum))
//' @param x A matrix.
//' @return A matrix of same size as x.
//' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
//' @examples
//' x <- matrix(1:8,ncol=2)
//' rowCumSum(x)
//' @export
// [[Rcpp::export]]
NumericMatrix rowCumSum(NumericMatrix x){
  arma::mat m(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat result;
  result=cumsum(m,1);
  return wrap(result);
}

//' Apply cumprod in each row 
//'
//' @description Fast computation of t(apply(x,1,cumprod))
//' @param x A matrix.
//' @return A matrix of same size as x.
//' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
//' @examples
//' x <- matrix(1:8,ncol=2)
//' rowCumProd(x)
//' @export
// [[Rcpp::export]]
NumericMatrix rowCumProd(NumericMatrix x){
  arma::mat m(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat result;
  result=cumprod(m,1);
  return wrap(result);
}
