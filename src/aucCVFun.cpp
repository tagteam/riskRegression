#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
List aucLoobFun(IntegerVector IDCase, IntegerVector IDControl, NumericMatrix riskMat, LogicalMatrix splitMat, NumericVector weights){
  int nCases = IDCase.length();
  int nControls = IDControl.length();
  NumericVector ic0Case(nCases);
  NumericVector ic0Control(nControls);
  bool warn = false;
  for (int i = 0; i < nCases; i++){
    for (int j = 0; j < nControls; j++){
      int idCase=IDCase[i]-1; // R indexing to C++ indexing. 
      int idControl=IDControl[j]-1; // R indexing to C++ indexing.
      int B = splitMat.ncol();
      int ibij = 0;
      double aucij = 0;
      for (int u = 0; u<B;u++){
        if (splitMat(idCase,u) && splitMat(idControl,u)){
          ibij+=1;
          if (riskMat(idCase,u) > riskMat(idControl,u)){
            aucij += 1.0;
          }
          else if (riskMat(idCase,u) == riskMat(idControl,u)){
            aucij+=0.5;
          }
        }
      }
      if (ibij == 0){ // the pair is not oob
        warn = true;  
      }
      else {
        ic0Case[i] += weights[idCase]*weights[idControl]*aucij / ((double) ibij);
        ic0Control[j] += weights[idCase]*weights[idControl]*aucij / ((double) ibij);
      }
    }
  }
  return(List::create(Named("warn") = warn,
                      Named("ic0Case") = ic0Case,
                      Named("ic0Control") = ic0Control));
}
