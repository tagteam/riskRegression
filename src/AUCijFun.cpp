
// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix AUCijFun(NumericVector riskCase, NumericVector riskControl){

  int nCase = riskCase.size();
  int nControl = riskControl.size();
  DoubleVector Great(nControl), Equal(nControl);
  NumericMatrix out(nCase,nControl);
  
  for (int x = 0; x<nCase; x++){
    Great = riskCase[x] > riskControl;
    Equal = riskCase[x] == riskControl;
    out(x,_) = Great + Equal*0.5;
  }

  return out;
}

