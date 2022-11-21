// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// * colCenter (documentation)
//' @title Apply - by column
//' @description Fast computation of sweep(X, MARGIN = 1, FUN = "-", STATS = center)
//' @name colCenter_cpp
//' 
//' 
//' @param X A matrix.
//' @param center a numeric vector of length equal to the number of rows of \code{x}
//' 
//' @return A matrix of same size as X.
//' 
//' @author Brice Ozenne <broz@@sund.ku.dk>
//' @examples
//' x <- matrix(1,6,5)
//' sweep(x, MARGIN = 1, FUN = "-", STATS = 1:6)
//' colCenter_cpp(x, 1:6 )

// * colCenter (code)
//' @rdname colCenter_cpp
//' @export
// [[Rcpp::export]]
arma::mat colCenter_cpp(arma::mat X, const arma::colvec& center){
  X.each_col() -= center;
  return(X);
}
