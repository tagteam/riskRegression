// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat predictCIF_cpp(const std::vector<arma::mat>& hazard, const std::vector<arma::mat>& cumHazard, const arma::mat& eXb_h, const arma::mat& eXb_cumH, const arma::mat& strata, 
                         const std::vector<double>& newtimes, const std::vector<double>& etimes, const std::vector<double>& etimeMax,
                         int nTimes, int nNewTimes, int nData, int cause, int nCause){
  
  arma::mat pred_CIF(nData, nNewTimes);
  pred_CIF.fill(0);
  
  double hazard_it;
  double AllcumHazard_it;
  size_t iterP; // index of the prediction time
  rowvec strataI(nCause);
  
  for(size_t iterI=0 ; iterI<nData; iterI++){ // index of the patient
    R_CheckUserInterrupt();
    
    iterP = 0;
    strataI = strata.row(iterI);
    
    for(size_t iterT=0 ; iterT<nTimes; iterT++){ // index of the time in the integral (event time number)
      // update position 
      while(iterP < nNewTimes && newtimes[iterP]<etimes[iterT]){
        iterP++;
        pred_CIF(iterI,iterP) = pred_CIF(iterI,iterP-1);
      }
      if(iterP >= nNewTimes){break;}
      // get hazard for the cause of interest
      hazard_it = hazard[cause](iterT,strataI[cause])*eXb_h(iterI,cause);
      
      // sum all cumHazard for all causes and exp the result
      AllcumHazard_it = 0; 
      for(size_t iterC=0 ; iterC<nCause; iterC++){
        AllcumHazard_it += cumHazard[iterC](iterT,strataI[iterC])*eXb_cumH(iterI,iterC);
      }
      // update the integral
      pred_CIF(iterI,iterP) += exp(-AllcumHazard_it) * hazard_it;  
      // Rcout << iterI << "*" << iterP << " | " << eXb(iterI,0) << " " << AllcumHazard_it << " " <<  exp(-AllcumHazard_it) << " " << hazard_it << endl;
    }
    
    if(iterP < nNewTimes){ // deal with prediction times after the last event
      
      if(newtimes[iterP]>etimeMax[iterI]){ // was the computation complete for this event
        pred_CIF(iterI,iterP) = NA_REAL;  
      }
      if(iterP < nNewTimes-1){
        for(size_t iterPP = iterP+1; iterPP<nNewTimes ; iterPP++){
          if(newtimes[iterPP]>etimeMax[iterI]){
            pred_CIF(iterI,iterPP) = NA_REAL;  
          }else{
            pred_CIF(iterI,iterPP) = pred_CIF(iterI,iterPP-1);  
          }
        }
      }
      
    }
    
    
  }
  
  return(pred_CIF);
}
