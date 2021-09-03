// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
/*
 *     V10 <- matrix(0, nrow = nCases, ncol = nauc)
 V01 <- matrix(0, nrow = nControls, ncol = nauc)
 tmn <- t(riskcontrols)
 tmp <- t(riskcases)
 for (i in 1:nCases) {
 V10[i, ] <- rowSums(tmn < tmp[, i]) + 0.5 * rowSums(tmn == tmp[, i])
 }
 for (i in 1:nControls) {
 V01[i, ] <- rowSums(tmp > tmn[, i]) + 0.5 * rowSums(tmp == tmn[, i])
 }
 browser()
 */

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
