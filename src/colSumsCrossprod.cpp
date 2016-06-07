// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Apply crossprod and colSums 
//'
//' @description Fast computation of crossprod(colSums(X),Y) 
//' @param X A matrix with dimensions k*n. Hence the result of \code{colSums(X)} has length n.
//' @param Y A matrix with dimenions n*m. Can be a matrix with dimension m*n but then \code{transposeY} should be \code{TRUE}.
//' @param transposeY Logical. If \code{TRUE} transpose Y before matrix multiplication.
//' @return A vector of length m.
//' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
//' @examples
//' x <- matrix(1:8,ncol=2)
//' y <- matrix(1:16,ncol=8)
//' colSumsCrossprod(x,y,0)
//' 
//' x <- matrix(1:8,ncol=2)
//' y <- matrix(1:16,ncol=2)
//' colSumsCrossprod(x,y,1)
//' @export
// [[Rcpp::export]]
NumericMatrix colSumsCrossprod(NumericMatrix X, NumericMatrix Y, bool transposeY){
  arma::mat A(X.begin(), X.nrow(), X.ncol(), false);
  arma::mat B(Y.begin(), Y.nrow(), Y.ncol(), false);
  arma::rowvec result;
  // result of colSums(A) has to be a matrix
  // with one row and as many columns as B has rows
  if (transposeY)
    result = arma::sum(A,0)*B.t();
  else
    result = arma::sum(A,0)*B;
  return wrap(result); 
}



 
