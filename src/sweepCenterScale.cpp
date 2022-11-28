// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"

using namespace Rcpp;
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

// * rowCenter (documentation)
//' @title Apply - by row
//' @description Fast computation of sweep(X, MARGIN = 2, FUN = "-", STATS = center)
//' @name rowCenter_cpp
//' 
//' @param X A matrix.
//' @param center a numeric vector of length equal to the number of rows of \code{x}
//' 
//' @return A matrix of same size as X.
//' @author Brice Ozenne <broz@@sund.ku.dk>
//' @examples
//' x <- matrix(1,6,5)
//' sweep(x, MARGIN = 2, FUN = "-", STATS = 1:5)
//' rowCenter_cpp(x, 1:5 )
//' 
//' rowCenter_cpp(x, colMeans(x) )

// * rowCenter (code)
//' @rdname rowCenter_cpp
//' @export
// [[Rcpp::export]]
arma::mat rowCenter_cpp(arma::mat X, const arma::rowvec& center){
  X.each_row() -= center;
  return(X);
}


// * colScale_cpp (documentation)
//' @title Apply / by column
//' @description Fast computation of sweep(X, MARGIN = 1, FUN = "/", STATS = scale)
//' @name colScale_cpp
//'
//' @param X A matrix.
//' @param scale a numeric vector of length equal to the number of rows of \code{x}
//' 
//' @return A matrix of same size as X.
//' @author Brice Ozenne <broz@@sund.ku.dk>
//' @examples
//' x <- matrix(1,6,5)
//' sweep(x, MARGIN = 1, FUN = "/", STATS = 1:6)
//' colScale_cpp(x, 1:6 )


// * colScale_cpp (code)
//' @rdname colScale_cpp
//' @export
// [[Rcpp::export]]
arma::mat colScale_cpp(arma::mat X, const arma::colvec& scale){
  X.each_col() /= scale;
  return(X);
}

// * rowScale_cpp (documentation)
//' @title Apply / by row
//' @description Fast computation of sweep(X, MARGIN = 2, FUN = "/", STATS = scale)
//' @name rowScale_cpp
//' 
//' @param X A matrix.
//' @param scale a numeric vector of length equal to the number of rows of \code{x}
//' 
//' @return A matrix of same size as X.
//' @author Brice Ozenne <broz@@sund.ku.dk>
//' @examples
//' x <- matrix(1,6,5)
//' sweep(x, MARGIN = 2, FUN = "/", STATS = 1:5)
//' rowScale_cpp(x, 1:5 )
//' 
//' rowScale_cpp(x, colMeans(x) )

// * rowScale_cpp (code)
//' @rdname rowScale_cpp
//' @export
// [[Rcpp::export]]
arma::mat rowScale_cpp(arma::mat X, const arma::rowvec& scale){
  X.each_row() /= scale;
  return(X);
}

// * colMultiply_cpp (documentation)
//' @title Apply * by column
//' @description Fast computation of sweep(X, MARGIN = 1, FUN = "*", STATS = scale)
//' @name colMultiply_cpp
//' 
//' @param X A matrix.
//' @param scale a numeric vector of length equal to the number of rows of \code{x}
//' 
//' @return A matrix of same size as X.
//' @author Brice Ozenne <broz@@sund.ku.dk>
//' @examples
//' x <- matrix(1,6,5)
//' sweep(x, MARGIN = 1, FUN = "*", STATS = 1:6)
//' colMultiply_cpp(x, 1:6 )
 
// * colMultiply_cpp (code)
//' @name colMultiply_cpp
//' @export
// [[Rcpp::export]]
arma::mat colMultiply_cpp(arma::mat X, const arma::colvec& scale){
  X.each_col() %= scale;
  return(X);
}

// * rowMultiply_cpp (documentation)
//' @title Apply * by row
//' @description Fast computation of sweep(X, MARGIN = 2, FUN = "*", STATS = scale)
//' @name rowMultiply_cpp
//' 
//' @param X A matrix.
//' @param scale a numeric vector of length equal to the number of rows of \code{x}
//' 
//' @return A matrix of same size as X.
//' @author Brice Ozenne <broz@@sund.ku.dk>
//' @examples
//' x <- matrix(1,6,5)
//' sweep(x, MARGIN = 2, FUN = "*", STATS = 1:5)
//' rowMultiply_cpp(x, 1:5 )
//' 
//' rowMultiply_cpp(x, 1/colMeans(x) )
//' 
// * rowMultiply_cpp (code)
//' @name rowMultiply_cpp
//' @export
// [[Rcpp::export]]
arma::mat rowMultiply_cpp(arma::mat X, const arma::rowvec& scale){
  X.each_row() %= scale;
  return(X);
}

