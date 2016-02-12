// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream> 
#include <Rmath.h> 

using namespace Rcpp;
// using namespace std;
// using namespace arma;

// const double CST_EPSILON = std::numeric_limits < double>::epsilon();
// NumericVector caclEventtimes_cpp(const IntegerVector& status, const NumericVector& alltimes, int n_patients, int cause);
// NumericVector sortNV_hpp(NumericVector x) ;
//' @export
// [[Rcpp::export]]
List BaselineHazardFast_cpp(const IntegerVector& status, const NumericVector& pattimes, const NumericVector& Xb, int n_patients,
                            const NumericVector& eventtimes, int cause){
  
  // NOTE : status, alltimes, Xb must refer to the patients in the same order
  //        eventtimes must be sorted in ascending order
  
  //// intialisation
  // sort status, pattimes, Xb in ascending order relative to pattimes
  // sort eventtimes in ascending order
  int n_eventtimes = eventtimes.size();
  
  double W, W_next = 0;
  for(int iter_pat = 0 ; iter_pat < n_patients ; iter_pat++){
    W_next += exp(Xb[iter_pat]);
  }
  NumericVector Lambda0(n_eventtimes+1);
  
  double death;
  
  //// main: loop over eventtimes
  Lambda0[0] = 0;
  
  int iter_pat = 0 ;
  
  for( int iter_time = 0 ; iter_time < n_eventtimes ; iter_time ++){
   W = W_next;
    
    // patients with event at time t
    death = 0;
    do{
      
      if(status[iter_pat] == cause){
        death ++;
      }
      W_next -= exp(Xb[iter_pat]); // no more at risk at for next time (already censored or dead)
      iter_pat++;
        
    }while(iter_pat < n_patients && pattimes[iter_pat] == eventtimes[iter_time]);
   
    // udpate Lambda0
    Lambda0[iter_time+1] = Lambda0[iter_time] + death / W;
    
  }
  
  //// export
  Lambda0.erase(Lambda0.begin());
  
  return(List::create(Named("hazard")  = Lambda0,
                      Named("time")  = eventtimes)
  );
  
}
