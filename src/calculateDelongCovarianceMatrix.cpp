// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Calculates the asymptotic covariance matrix for the delongtest function.
// First function is useful in the case where dolist has length 0 and we don't want the entire covariance matrix
// [[Rcpp::export]]
NumericVector delongtestHelper(int nauc, int nCases, int nControls, NumericMatrix tmn, NumericMatrix tmp){
  //std::cout << nauc << "\n\n";;
  NumericVector SEOfS(nauc);
  for (int j = 0; j < nauc; j++){
    vec temp1(nCases,fill::zeros);
    for (int r = 0; r < nCases; r++){
      for (int k = 0; k < nControls; k++){
        if (tmn(j,k) < tmp(j,r)) {
          temp1[r]+=1.0;
        }
        else if (tmn(j,k) == tmp(j,r)) {
          temp1[r]+=0.5;
        }
      }
    }
    //std::cout << "Temp1: " << temp1.n_rows << "\n";
    //temp1.print();
    //std::cout << "\n\n\n\n\n\n\n";

    vec temp2(nControls,fill::zeros);
    for (int r = 0;  r < nControls; r++){
      for (int k = 0; k < nCases; k++){
        if (tmp(j,k) > tmn(j,r)) {
          temp2[r]+=1.0;
        }
        else if (tmp(j,k) == tmn(j,r)) {
          temp2[r]+=0.5;
        }
      }
    }
    //std::cout << arma::var(temp1);

    //temp2.print();
    /*std::cout << "Temp2: ";
    temp1.print();
    std::cout << "\n";
    std::cout << "length of it: " << temp1.n_rows;*/

    SEOfS[j] = sqrt(arma::var(temp1)/(nCases*nControls*nControls) + arma::var(temp2)/(nCases*nCases*nControls));
    //std::cout << "SE: " << SEOfS[j] << "\n";
  }
  return wrap(SEOfS);
}


// [[Rcpp::export]]
NumericMatrix rowSumsAlt1(NumericMatrix V, NumericMatrix tmn, NumericMatrix tmp) {
  for (int r = 0; r < V.nrow(); r++) {
    uvec ans(V.ncol(),fill::zeros);
    for (int j = 0; j < tmn.ncol(); j++) {
      for (int i = 0; i < tmn.nrow(); i++) {
        if (tmn(i,j)<tmp(i,r)) {
          ans[i] += 1.0;
        }
        else if (tmn(i,j) == tmp(i,r)) {
          ans[i] += 0.5;
        }
      }
    }
    for (int j = 0; j < V.ncol();j++){
      V(r,j) = ans(j);
    }
  }
  return wrap(V);
}

// [[Rcpp::export]]
NumericMatrix rowSumsAlt2(NumericMatrix V, NumericMatrix tmn, NumericMatrix tmp) {
  for (int r = 0; r < V.nrow(); r++) {
    uvec ans(V.ncol(),fill::zeros);
    for (int j = 0; j < tmp.ncol(); j++) {
      for (int i = 0; i < tmp.nrow(); i++) {
        if (tmp(i,j)>tmn(i,r)) {
          ans[i] += 1.0;
        }
        else if (tmp(i,j) == tmn(i,r)) {
          ans[i] += 0.5;
        }
      }
    }
    for (int j = 0; j < V.ncol();j++){
      V(r,j) = ans(j);
    }
  }
  return wrap(V);
}

