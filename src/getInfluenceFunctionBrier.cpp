// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"
#include "IC-Nelson-Aalen-cens-time.h"

using namespace Rcpp;
using namespace arma;

// Calculate influence function for competing risk case/survival case with Nelson-Aalen censoring.
// see https://github.com/eestet75/riskRegressionStudy/blob/master/PicsForImplementation/BrierTrainTest.png
// Should be used with loob estimates and generally 
// [[Rcpp::export(rng=false)]]
NumericVector getInfluenceFunctionBrierKMCensoringUseSquared(double tau,
                                                             NumericVector time,
                                                             NumericVector residuals,
                                                             NumericVector status) {
  checkNAs(tau, GET_VARIABLE_NAME(tau));
  checkNAs(time, GET_VARIABLE_NAME(time));
  checkNAs(residuals, GET_VARIABLE_NAME(residuals));
  checkNAs(status, GET_VARIABLE_NAME(status));

  int n = time.size();
  NumericVector ic(n);
  arma::uvec sindex(n,fill::zeros);
  arma::vec utime=unique(time);
  int nu=utime.size();
  arma::vec atrisk(nu);
  arma::vec MC_term2(nu,fill::zeros);
  getInfluenceFunctionKM(time,status,atrisk,MC_term2,sindex,utime);
  
  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  auto lower = std::upper_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower) -1;
  double icpart1 = 0;
  double icpart2 = 0;
  double icpart = 0;
  for (int i = 0; i < n; i++){
    if ((status[i] == 1 || status[i] == 2) && time[i] <= tau){
      icpart2 += residuals[i];
    }
    else if (time[i] > tau){
      icpart += residuals[i];
    }
  }
  icpart = icpart / ( (double) n);
  double brier = mean(residuals);
  int tieIter = 0;
  // can do while loops together
  while ((tieIter < n) && (time[tieIter] == time[0])) {
    if ((time[tieIter] <= tau) && (status[tieIter]==1 || status[tieIter]==2)){
      icpart2 -= residuals[tieIter];
    }
    tieIter++;
  }
  int upperTie = tieIter-1;
  double icterm{}, icterm2{}, fihattau{};
  int j;
  for (int i = 0; i<n; i++){
    if (i > firsthit){
      j = firsthit;
    }
    else {
      j = i;
    }
    if (j==-1){
      fihattau = 0.0;
    }
    else if (utime[sindex[j]] < time[i]){
      fihattau = - MC_term2[sindex[j]];
    }
    else {
      fihattau =  (1-(status[i] != 0))*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
    }
    icterm = 1.0 / n * (icpart1+fihattau*icpart2);
    icterm2 = icpart * fihattau;
    if (upperTie == i){
      int tieIter = i+1;
      while ((tieIter < n) && (time[tieIter] == time[i+1])) {
        if ((time[tieIter] <= tau) && (status[tieIter]==1 || status[tieIter]==2)){
          icpart1 -= residuals[tieIter]*MC_term2[sindex[i]];
          icpart2 -= residuals[tieIter];
        }
        tieIter++;
      }
      upperTie = tieIter-1;
    }
    ic[i] = residuals[i] - brier + icterm+icterm2;
  }
  return ic;
}