// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <iostream> 
#include <Rmath.h> 

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

const double CST_EPSILON = std::numeric_limits < double>::epsilon();
NumericVector caclEventtimes_cpp(const IntegerVector& status, const NumericVector& alltimes, int n_patients, int cause);
NumericVector sortNV_hpp(NumericVector x) ;

//  1 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List BaselineHazard_cpp(const IntegerVector& status, const NumericVector& alltimes, int cause,
                        const NumericVector& Xb, NumericVector times){
  
  // NOTE : status, alltimes, cause, Xb must refer to the patients in the same order
  
  //// intialisation
  int n_patients = Xb.size() ; 
  int n_times = times.size() ; 
  times = sortNV_hpp(times);
  
  // define eventtimes
  NumericVector eventtimes = caclEventtimes_cpp(status, alltimes, n_patients, cause);
  eventtimes = sortNV_hpp(eventtimes);
  int n_eventtimes = eventtimes.size() ; 
  
  // prepare Lambda
  double W = 0, Lambda_tempo = 0;
  for(int iter_pat = 0 ; iter_pat < n_patients ; iter_pat++){
    W += exp(Xb[iter_pat]);
  }
  NumericVector Lambda0(n_times);
  
  double d;
  IntegerVector atRisk(n_patients,1);
  double time_tempo;
  
  int iter_time = 0;
  int iter_eventtime = 0;
  bool test_eqtime;
  bool test_inftime;
  
  //// main: loop over eventtimes
  
  // before the first event
  while(iter_time < n_times && times[iter_time] < eventtimes[0]){
    Lambda0[iter_time] = 0;
    iter_time++;
  }
  
  // first to last event
  while((iter_time < n_times) && (iter_eventtime < n_eventtimes)){
    time_tempo = eventtimes[iter_eventtime];
    
    // update W and d
    d = 0;
    for(int iter_pat = 0 ; iter_pat < n_patients ; iter_pat++){
      
      if(status[iter_pat] == cause && alltimes[iter_pat] == time_tempo){ // dead at time_tempo
        d ++;
      } else if(alltimes[iter_pat] < time_tempo && atRisk[iter_pat] == 1){ // no more at risk at time_tempo (already censored or dead)
        W -= exp(Xb[iter_pat]);
        atRisk[iter_pat] = 0;
      }
      
    }
    
    // udpate Lambda0
    Lambda_tempo += d / W;
    
    // update Lambda0 if requested time in [ eventtimes[iter_eventtime] ; eventtimes[iter_eventtime + 1] [
    test_eqtime = times[iter_time] == time_tempo ;
    test_inftime = (iter_eventtime + 1) < n_eventtimes &&  times[iter_time] < eventtimes[iter_eventtime + 1];
    
    while((iter_time < n_times) && (test_eqtime || test_inftime) ){ 
      Lambda0[iter_time] = Lambda_tempo ;
      
      iter_time++;
      if(iter_time < n_times){
        test_eqtime = times[iter_time] == time_tempo ;
        test_inftime = (iter_eventtime + 1) < n_eventtimes &&  times[iter_time] < eventtimes[iter_eventtime + 1];
      }
      
    }
    
    // move to next time
    iter_eventtime++;
  }
  
  // after the first event
  while(iter_time < n_times){
    Lambda0[iter_time] = Lambda_tempo;
    iter_time++;
  }
  
  
  //// export
  return(List::create(Named("hazard")  = Lambda0,
                      Named("time")  = times)
  );
  
  
  //   return(List::create(Named("hazard")  = Lambda0,
  //                       Named("d")  = d,
  //                       Named("W")  = W,
  //                       Named("atRisk")  = atRisk,
  //                       Named("Lambda_tempo")  = Lambda_tempo,
  //                       Named("time")  = times,
  //                       Named("alltimes")  = alltimes, 
  //                       Named("eventtimes")  = eventtimes));
}


//  2 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector caclEventtimes_cpp(const IntegerVector& status, const NumericVector& alltimes, int n_patients, int cause){
  
  NumericVector eventtimes(0);
  bool test_tie;
  
  for(int iter_pat=0 ; iter_pat < n_patients ; iter_pat++){
    
    if(status[iter_pat] == cause){
      
      // check ties
      test_tie = false;
      for(int iter_time = 0 ; iter_time < eventtimes.size(); iter_time++){
        if(eventtimes[iter_time] == alltimes[iter_pat]){
          test_tie = true;
          break;
        }
      }
      if(test_tie == false){
      eventtimes.push_back (alltimes[iter_pat]);
      }
      
    }
    
  }
  
  return(eventtimes);
}

//  3 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector sortNV_hpp(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}
