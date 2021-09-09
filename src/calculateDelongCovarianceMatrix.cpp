// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Calculates the asymptotic covariance matrix for the delongtest function.

// [[Rcpp::export]]
List calculateDelongSE(NumericMatrix risk, NumericVector nauc, NumericVector cause, NumericVector response, NumericVector doList, bool keepVCOV, NumericVector& alpha, bool fitSE) {
  Rcpp::Environment myEnv = Rcpp::Environment::global_env();
  if (fitSE || doList.length() > 0) {

    // calculating nCases and nControls

    //bool cases [cause.length()];
    //bool controls [cause.length()];
    unsigned int nCases = 0;
    unsigned int nControls = 0;
    for (int i = 0; i < cause.length(); i++){
        if (cause[i] == response[i]){
          //cases[i] = true;
          nCases++;
        }
        else {
          //controls[i] = true;
          nControls++;
        }
    }

    // calculate tmn, tmp
    NumericMatrix tmn(nauc[1], nControls);
    NumericMatrix tmp(nauc[1], nCases);
    int j = 0;
    int k = 0;
    for (int i = 0; i < cause.length(); i++){
      if (cause[i] == response[i]){
        tmn(_,j) = risk(i,_);
        j++;
      }
      else {
        tmp(_,k) = risk(i,_);
        k++;
      }
    }

    // easier method if we only need the variances and not the covariances
    if (!keepVCOV && doList.length() == 0){
      vec varOfS(nauc[1],fill::zeros);
      for (int j = 0; j < nauc[1]; j++){
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

        // TBD
        vec temp2(nCases,fill::zeros);
        for (int r = 0; r < nCases; r++){
          for (int k = 0; k < nControls; k++){
            if (tmn(j,k) < tmp(j,r)) {
              temp2[j]+=1.0;
            }
            else if (tmn(j,k) == tmp(j,r)) {
              temp2[j]+=0.5;
            }
          }
        }

        varOfS[j] = arma::var(temp1)/(nCases*nControls*nControls) + arma::var(temp2)/(nCases*nCases*nControls);
      }
      List L = List::create(varOfS);
      return wrap(L);
    }


    else if



  }

  uvec Cases = response == cause

  if (needEntireCovariance){

  }
  else {
    // needs to be implemented
    return wrap(cause);
  }
}

NumericMatrix rowSumsAlt1(NumericMatrix& V, NumericMatrix& tmn, NumericMatrix& tmp) {
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



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

