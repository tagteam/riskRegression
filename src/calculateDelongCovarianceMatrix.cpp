// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// C++ functions for calculating the asymptotic covariance matrix for the delongtest function.

// calculates midrank
vec calculateMidrank(vec& z){
  int m = z.size();
  vec wtemp = sort(z);
  uvec index = sort_index(z);
  vec w(m+1);
  for (int i = 0; i <m;i++ ){
    w[i]=wtemp[i];
  }
  w[m]=wtemp[m-1]+1;
  vec t(m);
  int i = 0;
  int a, b, j;
  while (i < m){
    a = i;
    j = a;
      while (w[j]==w[a]){
        j+=1;
      }
      b = j - 1;
    for (int k = a; k <= b; k++){
      t[k] = (double)(a+b)/2;
    }
    i = b+1;
  }
  int k;
  vec tk(m);
  for (i = 0; i < m; i++){
    k = index[i];
    tk[k] = t[i]+1;
  }
  return tk;
}

// Fast implementation of the calculation of the covariance matrix.
// Number of rows is the number of observations, for X and Y respectively
// Number of columns is the number of experiments
// we note here that pointers and references should be used in order
// to make efficient use of memory
// [[Rcpp::export]]
NumericMatrix calculateDelongCovarianceFast(NumericMatrix& Xs, NumericMatrix& Ys){
  int m = Xs.nrow();
  int n = Ys.nrow();
  if (Xs.ncol()!=Ys.ncol()){
    stop("Incompatible matrix dimensions.");
  }
  int k = Xs.ncol();
  mat V10(k,m);
  mat V01(k,n);
  //can be used for more efficient computation of cov. matrix, but that code isn't working
  for (int r = 0; r < k; r++){
    // Make them into armadillo vectors; might be an inefficient and superfluous operation
    // strangely enough we cannot write Xr = as<vec>(Xs(_,r))
    // for now it seems to work well enough
    NumericVector Xx = Xs(_,r);
    NumericVector Yy = Ys(_,r);
    vec Xr = as<vec>(Xx);
    vec Yr = as<vec>(Yy);

    // concatenate
    vec Zr = join_cols(Xr,Yr);
    // calculate midranks
    vec TZr = calculateMidrank(Zr);
    vec TXr = calculateMidrank(Xr);
    vec TYr = calculateMidrank(Yr);
    for (int i = 0; i < m; i++){
      V10(r,i)=(TZr[i]-TXr[i])/((double) n);
    }
    for (int j = 0; j < n; j++){
      V01(r,j)=1.0-(TZr[j+m]-TYr[j])/((double) m);
    }
  }
  mat S(k,k);
  mat s10 = arma::cov(V10.t());
  mat s01 = arma::cov(V01.t());
  S = s01/((double) n)+s10/((double) m);
  return wrap(S);
}
