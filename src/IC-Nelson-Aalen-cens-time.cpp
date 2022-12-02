// [[Rcpp::depends(RcppArmadillo)]]
#include "IC-Nelson-Aalen-cens-time.h"

using namespace Rcpp;
using namespace arma;

// The above for usage in C++ functions, as the above function is not memory efficient.
void getInfluenceFunctionKM(Rcpp::NumericVector& time, Rcpp::NumericVector& status,arma::vec& atrisk,arma::vec& MC_term2,arma::uvec& sindex,arma::vec& utime){
  int nu = atrisk.size();
  int n = time.size();
  arma::vec Cens(nu,fill::zeros);
  arma::vec hazardC(nu,fill::zeros);
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
}
