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
  pred_CIF.fill(NA_REAL);
  
  double hazard_it;
  double Allcumhazard_it;
  double CIF_it;
  double survival_t0=1;
  int iterP; // index of the prediction time
  rowvec strataI(nCause);
   
  for(int iterI=0 ; iterI<nData; iterI++){ // index of the patient
    R_CheckUserInterrupt();
    
    CIF_it = 0;
    iterP = 0;
    strataI = strata.row(iterI);

    for(int iterT=0 ; iterT<nEventTimes; iterT++){ // index of the time in the integral (event time number)
      // update position 
      while(iterP < nNewTimes && newtimes[iterP]<etimes[iterT]){
        if(newtimes[iterP] <= etimeMax[iterI]){
          pred_CIF(iterI,iterP) = CIF_it;
        }
        iterP++;
      }
      
      // if CIF has been calculated for all patients no need to continue to loop
      // if the next prediction time is after the last event no need to continue (all NA)
      if(iterP >= nNewTimes || newtimes[iterP] > etimeMax[iterI]){
        break;
       } 
      
      // get hazard for the cause of interest
      hazard_it = hazard[cause](iterT,strataI[cause])*eXb_h(iterI,cause);
      
      // sum all cumhazard for all causes and exp the result
      Allcumhazard_it = 0; 
      for(int iterC=0 ; iterC<nCause; iterC++){
        Allcumhazard_it += cumhazard[iterC](iterT,strataI[iterC])*eXb_cumH(iterI,iterC);
      }
      
      // update the integral
      if(R_IsNA(t0)){
        CIF_it += exp(-Allcumhazard_it) * hazard_it;  
      }else{// [only for conditional CIF]
        
        // get the survival up to t0 i.e. the survival at etimes just before t0
        if(etimes[iterT]<t0 && iterT<(nEventTimes-1) && etimes[iterT+1]>=t0){
          // NOTE: if iterT = nEventTimes-1 and etimes[iterT]<t0 then the landmark (t0) is after the last event so the CIF will be set to NA (since always etimes[iterT] < t0)
          survival_t0 = exp(-Allcumhazard_it); 
        }
        
        if(etimes[iterT] >= t0){ // not needed  newtimes[iterP]>=t0  because newtimes >= etimes see update position above 
          // Rcout << hazard_it << " ("<< iterP << ","<< etimes[iterT] << ","<< newtimes[iterP] << ")";
          CIF_it += exp(-Allcumhazard_it) * hazard_it / survival_t0;
        }
        
      }
    }
    
    
    if(iterP < nNewTimes){ // deal with prediction times before or equal to the last event 
     //> censored event are not in etimes thus prediction time after the last death and before the last censored event should be CIF_it and not NA
     //> prediction time exactly equal to the last event will not be assigned any value in the previous loop (because newtimes[iterP]<etimes[iterT]). It will be updated here.
     
      for(int iterPP = iterP; iterPP<nNewTimes ; iterPP++){
        if(newtimes[iterPP] <= etimeMax[iterI]){
          pred_CIF(iterI,iterPP) = CIF_it;  
        }else{
          break;
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