// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat predictCIF_cpp(const std::vector<arma::mat>& hazard, 
                         const std::vector<arma::mat>& cumhazard, 
                         const arma::mat& eXb_h, 
                         const arma::mat& eXb_cumH, 
                         const arma::mat& strata, 
                         const std::vector<double>& newtimes, 
                         const std::vector<double>& etimes, 
                         const std::vector<double>& etimeMax, 
                         double t0,
                         int nEventTimes, 
                         int nNewTimes, 
                         int nData, 
                         int cause, 
                         int nCause){
  
  arma::mat pred_CIF(nData, nNewTimes);
  pred_CIF.fill(0);
  
  double hazard_it;
  double Allcumhazard_it;
  double survival_t0=1;
  int iterP; // index of the prediction time
  rowvec strataI(nCause);
  
  for(int iterI=0 ; iterI<nData; iterI++){ // index of the patient
    R_CheckUserInterrupt();
    
    iterP = 0;
    strataI = strata.row(iterI);
    
    for(int iterT=0 ; iterT<nEventTimes; iterT++){ // index of the time in the integral (event time number)
      // update position 
      while(iterP < nNewTimes && newtimes[iterP]<etimes[iterT]){
        iterP++;
        pred_CIF(iterI,iterP) = pred_CIF(iterI,iterP-1);
      }
      if(iterP >= nNewTimes){break;}
      // get hazard for the cause of interest
      hazard_it = hazard[cause](iterT,strataI[cause])*eXb_h(iterI,cause);
      
      // sum all cumhazard for all causes and exp the result
      Allcumhazard_it = 0; 
      for(int iterC=0 ; iterC<nCause; iterC++){
        Allcumhazard_it += cumhazard[iterC](iterT,strataI[iterC])*eXb_cumH(iterI,iterC);
      }
      
      // update the integral
      if(R_IsNA(t0)){
        pred_CIF(iterI,iterP) += exp(-Allcumhazard_it) * hazard_it;  
      }else{// [only for conditional CIF]
        
        if(etimes[iterT]<t0 && iterT<(nEventTimes+1) && etimes[iterT+1]>=t0){
          survival_t0 = exp(-Allcumhazard_it); // get the survival up to t0
        }
        
        if(etimes[iterT] >= t0){ // not needed  newtimes[iterP]>=t0  because newtimes >= etimes see update position above 
          // Rcout << hazard_it << " ("<< iterP << ","<< etimes[iterT] << ","<< newtimes[iterP] << ")";
          pred_CIF(iterI,iterP) += exp(-Allcumhazard_it) * hazard_it / survival_t0;
        }
        
      }
    }
    // Rcout << endl;
    
    if(iterP < nNewTimes){ // deal with prediction times after the last event
      
      if(newtimes[iterP]>etimeMax[iterI]){ // was the computation complete for this event
        pred_CIF(iterI,iterP) = NA_REAL;  
      }
      if(iterP < nNewTimes-1){
        for(int iterPP = iterP+1; iterPP<nNewTimes ; iterPP++){
          if(newtimes[iterPP]>etimeMax[iterI]){
            pred_CIF(iterI,iterPP) = NA_REAL;  
          }else{
            pred_CIF(iterI,iterPP) = pred_CIF(iterI,iterPP-1);  
          }
        }
      }
      
    }
    
    if(R_IsNA(t0) == false){ // before t0 fill with NA
      iterP = 0;
      while(iterP < nNewTimes && newtimes[iterP]<t0){
        pred_CIF(iterI,iterP) = NA_REAL;
        iterP++;
      }
      
    }
    
  }
  
  return(pred_CIF);
}