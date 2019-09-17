// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat predictCIF_cpp(const std::vector<arma::mat>& hazard, 
                         const std::vector<arma::mat>& cumhazard, 
                         const arma::mat& eXb,
                         const arma::mat& strata, 
                         const std::vector<double>& newtimes, 
                         const std::vector<double>& etimes, 
                         const std::vector<double>& etimeMax, 
                         double t0,
                         int nEventTimes, 
                         int nNewTimes, 
                         int nData, 
                         int cause, 
                         int nCause,
						 bool survtype,
                         bool productLimit,
						 bool diag){
  
  arma::mat pred_CIF;
  if (diag) {
	pred_CIF.resize(nData, 1);  
  }else{
	pred_CIF.resize(nData, nNewTimes);
  }
  
  pred_CIF.fill(NA_REAL);
  
  double hazard_it; // hazard for the cause of interest at time t for individual i
  double hazard_tempo; // hazard for the cause of interest at time t for individual i
  double survival_it; // overall survival at time t for individual i
  double CIF_it;// cumulative incidence at time t for individual i
  double survival_it0=1;// survival at time t0 for individual i
  int iterP; // index of the prediction time
  rowvec strataI(nCause);

  int iNNewTimes;
  vector<double> iNewTimes;
  if(diag){
	iNewTimes.resize(1);
  }
  
  for(int iterI=0 ; iterI<nData; iterI++){ // index of the patient
    R_CheckUserInterrupt();
    
    CIF_it = 0;
    iterP = 0;
    survival_it = 1;
    strataI = strata.row(iterI);

	if(diag){
	  iNNewTimes = 1;
	  iNewTimes[0] = newtimes[iterI];
	}else{
	  iNNewTimes = nNewTimes;
	  iNewTimes = newtimes;
	}
	
    for(int iterT=0 ; iterT<nEventTimes; iterT++){ // index of the time in the integral (event time number)
      // update position 
      while(iterP < iNNewTimes && iNewTimes[iterP]<etimes[iterT]){
        if(iNewTimes[iterP] <= etimeMax[iterI]){
		  pred_CIF(iterI,iterP) = CIF_it;
        }
        iterP++;
      }
      
      // if CIF has been calculated for all patients no need to continue to loop
      // if the next prediction time is after the last event no need to continue (all NA)
      if(iterP >= iNNewTimes || iNewTimes[iterP] > etimeMax[iterI]){
        break;
	  } 
      
      // get hazard for the cause of interest
      hazard_it = hazard[cause](iterT,strataI[cause])*eXb(iterI,cause);
      
      // sum all cumhazard for all causes times the linear predictor and then take the exponential
      if(iterT>0){ // it is survival at t- which is computed i.e. the survival at the previous eventtime (censoring does not affect survival)
        
        if(productLimit){ // product limit - equivalent to mstate

		  if(survtype){ // product limit 
            survival_it *= (1-eXb(iterI,1)*hazard[1](iterT-1,strataI[1]));
		  }else{	    
			hazard_tempo = 0;
			for(int iterC=0 ; iterC<nCause; iterC++){
			  hazard_tempo += eXb(iterI,iterC)*hazard[iterC](iterT-1,strataI[iterC]);
			}
			survival_it *= (1-hazard_tempo);
		  }
	  
		}else{
	  
		  if(survtype){
			survival_it = exp(-cumhazard[1](iterT-1,strataI[1])*eXb(iterI,1));
		  }else{	     
			survival_it = 0; 
			for(int iterC=0 ; iterC<nCause; iterC++){
			  survival_it += cumhazard[iterC](iterT-1,strataI[iterC])*eXb(iterI,iterC);
			}
			survival_it = exp(-survival_it);
		  }
	  
        }
	
      } // otherwise the survival stays at 1
      
      // update the integral
      if(R_IsNA(t0)){	
        CIF_it += survival_it * hazard_it;
	
      }else{// [only for conditional CIF]

		// get the survival up to t0 i.e. the survival at etimes just before t0
        if(etimes[iterT]>=t0 && ((iterT>1 && etimes[iterT-1]<t0) || iterT == 0)){
          // NOTE: if iterT = nEventTimes-1 and etimes[iterT]<t0 then the landmark (t0) is after the last event so the CIF will be set to NA (since always etimes[iterT] < t0)
          survival_it0 = survival_it;
		}
        
        if(etimes[iterT] >= t0){ // not needed  iNewTimes[iterP]>=t0  because iNewTimes >= etimes see update position above 
          // Rcout << hazard_it << " ("<< iterP << ","<< etimes[iterT] << ","<< iNewTimes[iterP] << ")";
	  CIF_it += survival_it * hazard_it / survival_it0;
        }
        
      }
    }
    
    
    if(iterP < iNNewTimes){ // deal with prediction times before or equal to the last event 
     //> censored event are not in etimes thus prediction time after the last death and before the last censored event should be CIF_it and not NA
     //> prediction time exactly equal to the last event will not be assigned any value in the previous loop (because iNewTimes[iterP]<etimes[iterT]). It will be updated here.
     
      for(int iterPP = iterP; iterPP<iNNewTimes ; iterPP++){
        if(iNewTimes[iterPP] <= etimeMax[iterI]){
          pred_CIF(iterI,iterPP) = CIF_it;  
        }else{
          break;
        }
      }
      
    }
    
    if(R_IsNA(t0) == false){ // before t0 fill with NA
      iterP = 0;
      while(iterP < nNewTimes && iNewTimes[iterP]<t0){
        pred_CIF(iterI,iterP) = NA_REAL;
        iterP++;
      }
      
    }
    
  }
  
  return(pred_CIF);
}
