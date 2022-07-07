
// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"
using namespace Rcpp;
//' Apply crossprod and rowSums
//'
//' @description Fast computation of crossprod(rowSums(X),Y)
//' @param X A matrix with dimensions n*k. Hence the result of \code{rowSums(X)} has length n.
//' @param Y A matrix with dimenions n*m. Can be a matrix with dimension m*n but then \code{transposeY} should be \code{TRUE}.
//' @param transposeY Logical. If \code{TRUE} transpose Y before matrix multiplication.
//' @return A vector of length m.
//' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
//' @examples
//' x <- matrix(1:10,nrow=5)
//' y <- matrix(1:20,ncol=4)
//' rowSumsCrossprod(x,y,0)
//'
//' x <- matrix(1:10,nrow=5)
//' y <- matrix(1:20,ncol=5)
//' rowSumsCrossprod(x,y,1)
//' @export
// [[Rcpp::export]]
NumericMatrix rowSumsCrossprod(NumericMatrix X, NumericMatrix Y, bool transposeY){
  arma::mat A(X.begin(), X.nrow(), X.ncol(), false);
  arma::mat B(Y.begin(), Y.nrow(), Y.ncol(), false);
  arma::rowvec result;
  // result of colSums(A) has to be a matrix
  // with one row and as many columns as B has rows
  // since sum(A,1) is a matrix with one column
  // we transpose before multiplication
  if (transposeY)
    result = arma::sum(A,1).t()*B.t();
  else
    result = arma::sum(A,1).t()*B;
  return wrap(result); 
}
