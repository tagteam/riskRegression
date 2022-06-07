// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"
using namespace Rcpp;
using namespace arma;

//' @title Influence function for Nelson-Aalen estimator.
//' 
//' @description Fast computation of influence function for Nelson-Aalen estimator of the censoring times
//' @param time sorted vector of event times. Sorted according to time and -status so that events come first a tied times.
//' @param status sorted vector of 0 = censored or 1 = event (any cause). Sorted according to time and -status so that events come first a tied times.
//' @return A square matrix where each column corresponds to a subject and each row to a time point. 
//' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
//' @examples
//' time = c(1,3,3,4)
//' status = c(1,0,1,1)
//' IC_Nelson_Aalen_cens_time(time,status)
//' @export
// [[Rcpp::export]]
NumericMatrix IC_Nelson_Aalen_cens_time(NumericVector time,
					NumericVector status){
  uword n = time.size();
  arma::uvec sindex(n,fill::zeros);
  arma::vec utime=unique(time);
  uword nu=utime.size();
  arma::vec atrisk(nu);
  arma::vec Cens(nu,fill::zeros);
  arma::vec hazardC(nu,fill::zeros);
  arma::vec MC_term2(nu,fill::zeros);
  arma::mat MC(nu,n);
  MC.zeros();
  //initialize first time point t=0 with data of subject i=0
  arma::uword t=0;
  double Y = (double) n;
  atrisk[0]=Y;
  Cens[0]=(1-status[0]);
  hazardC[0]=Cens[0]/Y;
  MC_term2[0]+=hazardC[0];
  //loop through time points until last subject i=(n-1)
  for (uword i=1;i<=n;i++) {
    if (i<n && time[i]==time[i-1]){// these are tied values    
	Cens[t] +=(1-status[i]);
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
	Cens[t]=(1-status[i]);
      }
    }
  }
  // for (uword t=0;t<nu;t++)Rprintf("t=%d\tutime[t]=%1.2f\tCens[t]=%1.2f\tatrisk[t]=%1.2f\thazardC[t]=%1.2f\tMC_term2[t]=%1.2f\t\n",t,utime[t],Cens[t],atrisk[t],hazardC[t],MC_term2[t]);
  // integrate compensator over time
  MC_term2 = arma::cumsum(MC_term2);
  // for (uword t=0;t<nu;t++) Rprintf("t=%d\tMC_term2[t]=%1.2f\t\n",t,MC_term2[t]);
  for (uword i=0;i<n;i++){
    for (uword t=0;t<nu;t++){
      // Rprintf("i=%d\tt=%d\tsindex[i]=%d\tCens[sindex[i]]=%1.2f\tatrisk[sindex[i]]=%1.2f\tMC_term2[sindex[i]]=%1.2f\tMC_term2[t]=%1.2f\t\n",i,t,sindex[i],Cens[sindex[i]],atrisk[sindex[i]],MC_term2[sindex[i]],MC_term2[t]);
      if (utime[t]<time[i]){
	MC(t,i) = - MC_term2[t];
      } else{
	MC(t,i) = (1-status[i])*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
      }
    }
  }
  return wrap(MC);
}


