// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Apply cumsum in each row 
//'
//' @description Fast computation of influence function for Nelson-Aalen estimator of the censoring times
//' @param time event times
//' @param status binary 
//' @return A square matrix where each column corresponds to a subject and each row to a time point. 
    //' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
    //' @examples
    //' time = c(1,3,4)
    //' status = c(1,0,1)
    //' IC_Nelson_Aalen_cens_time(time,status)
    //' @export
    // [[Rcpp::export]]
NumericMatrix IC_Nelson_Aalen_cens_time(NumericVector time,
					NumericVector status){
  double n = time.size();
  arma::vec utime=unique(time);
  double nu=utime.size();
  // intialize
  arma::vec atrisk(n);
  arma::vec hazardC(n,fill::zeros);
  // double MC2;
  arma::vec MC_term2(n,fill::zeros);
  arma::mat MC(nu,n);
  MC.zeros();
  for (uword t=0;t<n;t++) {atrisk[i]=n-t;}
  for (uword t=0;t<n;t++) {hazardC[t]=(1-status[t])/atrisk[t];}
  for (uword t=0;t<n;t++) {
    MC_term2[t]+=hazardC[t]*n/atrisk[t];
  }
  MC_term2 = arma::cumsum(MC_term2);
  for (uword i=0;i<n;i++){
    for (uword t=0;t<n;t++){
      if (t<i){
 	MC(t,i) = - MC_term2[t];
      } else{
 	MC(t,i) = (1-status[i])*n/atrisk[i]- MC_term2[i];
      }
    }
  }
  return wrap(MC);
}
