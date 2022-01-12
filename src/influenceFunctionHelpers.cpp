// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Same as colMeans(A*b) for a matrrix A and vector b
// [[Rcpp::export]]
NumericVector columnMeanWeight(NumericMatrix A, NumericVector x){
  int nrows = A.nrow();
  int ncols = A.ncol();
  NumericVector ans(ncols);
  for (int j = 0; j < ncols; j++){
    double sum = 0.0;
    for (int i = 0; i < nrows; i++){
      sum += A(i,j)*x[i];
    }
    ans[j] = sum/nrows;
  }
  return ans;
}


// Same as colSums(b*(A+1)) for a matrix A and vector b
// [[Rcpp::export]]
NumericVector T3CalculationHelper(NumericVector x, NumericMatrix A){
  int nrows = A.nrow();
  int ncols = A.ncol();
  NumericVector ans(ncols);
  for (int j = 0; j < ncols; j++){
    double sum = 0.0;
    for (int i = 0; i < nrows; i++){
      sum += (A(i,j)+1.0)*x[i];
    }
    ans[j] = sum;
  }
  return ans;
}
// 
// // [[Rcpp::export]]
// NumericMatrix htijCalculationHelper(NumericMatrix mcase, NumericMatrix mcontrol,NumericMatrix wcase,NumericMatrix wcontrol,NumericVector n){
//   int nrows = mcase.nrow();
//   int ncols = mcase.ncol();
//   NumericMatrix ans(nrows,ncols);
//   for (int j = 0; j < ncols; j++){
//     for (int i = 0; i < nrows; i++){
//       if (mcase(i,j) > mcontrol(i,j)) {
//         ans(i,j) = wcase(i,j)*wcontrol(i,j)*n[0]*n[0];
//         // if (i < 50 && j < 50) {
//         //   Rcout << "1 s to " << i << "," << j << " with value" << ans(i,j) << "\n";
//         // }
//       }
//       else if (mcase(i,j) == mcontrol(i,j)) {
//         ans(i,j) = 0.5* wcase(i,j)*wcontrol(i,j)*n[0]*n[0];
//       }
//       else {
//         ans(i,j) = 0.0;
//       }
//     }
//   }
//   return ans;
// }
// 

// [[Rcpp::export]]
NumericMatrix htijCalculationHelper(NumericVector mcase, NumericVector mcontrol,NumericVector wcase,NumericVector wcontrol,int n, int nrows, int ncols){
  NumericMatrix ans(nrows,ncols);
  for (int j = 0; j < ncols; j++){
    for (int i = 0; i < nrows; i++){
      if (mcase(i) > mcontrol(j)) {
        ans(i,j) = wcase(i)*wcontrol(j)*n*n;
        // if (i < 50 && j < 50) {
        //   Rcout << "1 s to " << i << "," << j << " with value" << ans(i,j) << "\n";
        // }
      }
      else if (mcase(i) == mcontrol(j)) {
        ans(i,j) = 0.5* wcase(i)*wcontrol(j)*n*n;
      }
      else {
        ans(i,j) = 0.0;
      }
    }
  }
  return ans;
}

// [[Rcpp::export]]
NumericMatrix rowSumsCrossprodSpec(arma::mat &X, arma::mat &Y){
  return(wrap(arma::sum(X,1).t()*(Y+1)));
}

