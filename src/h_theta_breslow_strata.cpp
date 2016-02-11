// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <iostream> 
#include <Rmath.h> 

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

const double CST_EPSILON = std::numeric_limits < double>::epsilon();
NumericVector sortNV_hpp(NumericVector x) ;

//  1 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List BaselineHazard_cpp(const IntegerVector& status, const NumericVector& alltimes, int cause,
                        const NumericVector& Xb, NumericVector times){
  
  //// for the moment do not handle ties
  
  // intialisation
  size_t n_patients = Xb.size() ; 
  
  times = sortNV_hpp(times);
  size_t n_times = times.size() ; 
  
  vector<double> eventtimes(0);
  for(size_t iter_pat=0 ; iter_pat < n_patients ; iter_pat++){
    if(status[iter_pat] == cause){
      eventtimes.push_back (alltimes[iter_pat]);
    }
  }
  std::sort (eventtimes.begin(), eventtimes.end()); 
  size_t n_eventtimes = eventtimes.size() ; 

  vector < double > Lambda0(n_times);
  
  double W = 0, Lambda_tempo = 0;
  for(size_t iter_pat = 0 ; iter_pat < n_patients ; iter_pat++){
    W += exp(Xb[iter_pat]);
  }
  
  double d;
  vector < bool > atRisk(n_patients,true);
  double time_tempo;
  
  size_t iter_time = 0;
  size_t iter_eventtime = 0;
  bool test_eqtime;
  bool test_inftime;
  
  // specific case before the first event
  if(iter_time < n_times && times[iter_time] < eventtimes[0]){
    Lambda0[iter_time] = 0;
    iter_time++;
  }
  
  // main: loop over eventtimes
  while((iter_time < n_times) && (iter_eventtime < n_eventtimes)){
    time_tempo = eventtimes[iter_eventtime];

    // update W and d
      d = 0;
      for(size_t iter_pat = 0 ; iter_pat < n_patients ; iter_pat++){
        
        if(status[iter_pat] == cause && alltimes[iter_pat] == time_tempo){ // dead at time_tempo
          d ++;
        } else if(alltimes[iter_pat] < time_tempo && atRisk[iter_pat] == true){ // no more at risk at time_tempo (already censored or dead)
          W -= exp(Xb[iter_pat]);
          atRisk[iter_pat] = false;
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
  
  
  // export
  return(List::create(Named("hazard")  = Lambda0,
                      Named("time")  = times)
  );
  
  
//   return(List::create(Named("Lambda0")  = Lambda0,
//                       Named("d")  = d,
//                       Named("W")  = W,
//                       Named("atRisk")  = atRisk,
//                       Named("Lambda_tempo")  = Lambda_tempo,
//                       Named("times")  = times,
//                       Named("alltimes")  = alltimes, 
//                       Named("eventtimes")  = eventtimes));
}


//  2 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
inline NumericVector sortNV_hpp(NumericVector x) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

//  2 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
vector<double> mainX () {
  int myints[] = {32,71,12,45,26,80,53,33};
  std::vector<double> myvector (myints, myints+8);               // 32 71 12 45 26 80 53 33
  
  // using default comparison (operator <):
  std::sort (myvector.begin(), myvector.end());           //(12 32 45 71)26 80 53 33
  
  return myvector;
}
