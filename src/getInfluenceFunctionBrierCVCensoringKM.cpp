// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export(rng=false)]]
NumericVector getInfluenceFunctionBrierCVCensoringKM(double tau,
                                                     NumericVector time,
                                                     NumericVector residuals,
                                                     NumericVector status) {
  // Thomas' code from IC of Nelson-Aalen estimator
  //initialize first time point t=0 with data of subject i=0
  int n = time.size();
  arma::uvec sindex(n,fill::zeros);
  arma::vec utime=unique(time);
  int nu=utime.size();
  arma::vec atrisk(nu);
  arma::vec Cens(nu,fill::zeros);
  arma::vec hazardC(nu,fill::zeros);
  arma::vec MC_term2(nu,fill::zeros);
  int t=0;
  double Y = (double) n;
  atrisk[0]=Y;
  Cens[0]=(1-(status[0] != 0));
  hazardC[0]=Cens[0]/Y;
  MC_term2[0]+=hazardC[0];
  //loop through time points until last subject i=(n-1)
  for (int i=1;i<=n;i++) {
    if (i<n && time[i]==time[i-1]){// these are tied values
      Cens[t] +=(1-(status[i] != 0));
      Y-=1;
      sindex[i]=t;    // index pointer from subject i to unique time point t
    }else{
      utime[t]=time[i-1];
      hazardC[t]=Cens[t]/atrisk[t];
      MC_term2[t]=hazardC[t]*n/atrisk[t];
      //initialize next time point with data of current subject i
      if (i<n){
        t++;
        sindex[i]=t;    // index pointer from subject i to unique time point t
        Y-=1;
        atrisk[t]=Y;
        Cens[t]=(1-(status[i] != 0));
      }
    }
  }
  MC_term2 = arma::cumsum(MC_term2);
  NumericVector icG(n);
  
  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  auto lower = std::upper_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower) -1;
  double icpart1 = 0;
  double icpart2 = 0;
  double icpart = 0;
  for (int i = 0; i < n; i++){
    if (status[i] == 1 && time[i] <= tau){
      icpart2 += residuals[i];
    }
    else if (time[i] > tau){
      icpart += residuals[i];
    }
  }
  icpart = icpart / ( (double) n);
  int tieIter = 0;
  // can do while loops together
  while ((tieIter < n) && (time[tieIter] == time[0])) {
    if ((time[tieIter] <= tau) && (status[tieIter]==1)){
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
        if ((time[tieIter] <= tau) && (status[tieIter]==1)){
          icpart1 -= residuals[tieIter]*MC_term2[sindex[i]];
          icpart2 -= residuals[tieIter];
        }
        tieIter++;
      }
      upperTie = tieIter-1;
    }
    icG[i] = icterm+icterm2;
  }
  return icG;
}