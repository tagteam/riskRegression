// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//' @title Apply - by column
//'
//' @description Fast computation of sweep(X, MARGIN = 1, FUN = "-", STATS = center)
//' 
//' @param X A matrix.
//' @param center a numeric vector of length equal to the number of rows of \code{x}
//' 
//' @return A matrix of same size as X.
//' @author Brice Ozenne <broz@@sund.ku.dk>
//' @examples
//' x <- matrix(1,6,5)
//' sweep(x, MARGIN = 1, FUN = "-", STATS = 1:6)
//' colCenter_cpp(x, 1:6 )
//' 
//' @export
// [[Rcpp::export]]
arma::mat colCenter_cpp(arma::mat X, arma::colvec& center){
  X.each_col() -= center;
  return(X);
}


//' @title Apply - by row
//'
//' @description Fast computation of sweep(X, MARGIN = 2, FUN = "-", STATS = center)
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
//' 
//' @export
// [[Rcpp::export]]
arma::mat rowCenter_cpp(arma::mat X, arma::rowvec& center){
  X.each_row() -= center;
  return(X);
}

//' @title Apply / by column
//'
//' @description Fast computation of sweep(X, MARGIN = 1, FUN = "/", STATS = scale)
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
//' 
//' @export
// [[Rcpp::export]]
arma::mat colScale_cpp(arma::mat X, arma::colvec& scale){
  X.each_col() /= scale;
  return(X);
}


//' @title Apply / by row
//'
//' @description Fast computation of sweep(X, MARGIN = 2, FUN = "/", STATS = scale)
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
//' 
//' @export
// [[Rcpp::export]]
arma::mat rowScale_cpp(arma::mat X, arma::rowvec& scale){
  X.each_row() /= scale;
  return(X);
}

//' @title Apply * by column
//'
//' @description Fast computation of sweep(X, MARGIN = 1, FUN = "*", STATS = scale)
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
//' 
//' @export
// [[Rcpp::export]]
arma::mat colMultiply_cpp(arma::mat X, arma::colvec& scale){
  X.each_col() %= scale;
  return(X);
}


//' @title Apply * by row
//'
//' @description Fast computation of sweep(X, MARGIN = 2, FUN = "*", STATS = scale)
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
//' @export
// [[Rcpp::export]]
arma::mat rowMultiply_cpp(arma::mat X, arma::rowvec& scale){
  X.each_row() %= scale;
  return(X);
}