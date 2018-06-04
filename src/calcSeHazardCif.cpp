// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


inline double calcIFhazard(double time,
			   double sampleTime,
			   const rowvec& IFbeta,
			   const colvec& Ehazard0,
			   const rowvec& X,
			   double hazard_iS0,
			   double newEXb,
			   double sampleEXb,
			   double hazard0,
			   double iS0,
			   bool sameStrata,
			   bool jumpTime,
			   int p,
			   bool hazard);
  
// {{{ calcSeHazard_cpp : on the log scale
// [[Rcpp::export]]
List calcSeHazard_cpp(const NumericVector& seqTau, // horizon time for the predictions
		      const IntegerVector& indexTau, // position of the prediction times relatively to the event times in the strata
		      const IntegerVector& indexJump, // position of the event time of the sample relatively to the event times in the strata
		      const arma::mat& IFbeta, // pre-computed iidCox
		      const arma::mat& cumEhazard0, // pre-computed iidCox
		      const NumericVector& iS0, // pre-computed iidCox
		      const NumericVector& cumhazard_iS0, // pre-computed iidCox
		      const NumericVector& newEXb, // linear predictor fit
		      const NumericVector& sampleEXb, // linear predictor new observation
		      const arma::mat& X, // covariates to condition on to make prediction
		      const LogicalVector& sameStrata, // strata of observations used to fit the Cox model
		      const NumericVector& sampleTime, // event times relative to the observations used to fit the Cox model
		      const NumericVector& cumhazard0, // baseline cumulative hazard
		      const arma::mat& newSurvival, //
		      double firstJumpTime, double lastSampleTime,
		      int nTau, int nNewObs, int nSample, int p,
		      bool exportSE, bool exportIF, bool exportIFsum_cumhazard, bool exportIFsum_survival){

  // define objects
  double IF_cumhazard;
  int iTau0=0,iTauMax=nTau-1;
  int indexTempo;
  
  arma::mat SEcumhazard;
  if(exportSE){
    SEcumhazard.resize(nNewObs, nTau);
    SEcumhazard.fill(0);
  }
  arma::cube IFcumhazard;
  if(exportIF){
    IFcumhazard.resize(nNewObs, nTau, nSample);
    IFcumhazard.fill(0);
  }
  arma::mat IFsum_cumhazard;
  if(exportIFsum_cumhazard){
    IFsum_cumhazard.resize(nSample, nTau);
    IFsum_cumhazard.fill(0);
  }
  arma::mat IFsum_survival;
  if(exportIFsum_survival){
    IFsum_survival.resize(nSample, nTau);
    IFsum_survival.fill(0);
  }
  
  // narrow prediction times
  while((iTau0 < nTau) && firstJumpTime>seqTau[iTau0]){ // start at the first event or after 
    iTau0++;    
  }
  while((iTauMax >= 0) && seqTau[iTauMax]>lastSampleTime){ // end at the last event or before
    iTauMax--;
  }

  // interesting times
  if(iTau0 < nTau && iTauMax >= 0){
    for(int iNewObs=0; iNewObs<nNewObs ; iNewObs++){
     
      R_CheckUserInterrupt();
      
      for(int iTime=iTau0; iTime<=iTauMax; iTime++){

	for(int iSample=0; iSample<nSample ; iSample++){
   	  indexTempo = std::min(indexJump[iSample],indexTau[iTime]);
	  IF_cumhazard = calcIFhazard(seqTau[iTime],
				      sampleTime[iSample],
				      IFbeta.row(iSample),
				      cumEhazard0.col(iTime),
				      X.row(iNewObs),
				      cumhazard_iS0[indexTempo],
				      newEXb[iNewObs],
				      sampleEXb[iSample],  
				      cumhazard0[iTime],
				      iS0[iSample],
				      sameStrata[iSample],
				      false, // wrong value but not used when computing the cumulative hazard
				      p, false);
	  if(exportSE){
	    SEcumhazard(iNewObs,iTime) += pow(IF_cumhazard,2);
	  }
	  if(exportIF){
	    IFcumhazard(iNewObs,iTime,iSample) = IF_cumhazard;	    
	  }
	  if(exportIFsum_cumhazard){
	    IFsum_cumhazard(iSample,iTime) += IF_cumhazard;	    
	  }
	  if(exportIFsum_survival){
	    IFsum_survival(iSample,iTime) += -IF_cumhazard*newSurvival(iNewObs,iTime);
	  }
	}
      }
    }

    // return standard error instead of variance
    if(exportSE){
      for(int iTime=iTau0; iTime<=iTauMax; iTime++){
	SEcumhazard.col(iTime) = sqrt(SEcumhazard.col(iTime));
      }      
    }
  }

  // set IF/variance to NA after the last event time
  if(iTauMax+1 < nTau){
    for(int iTime=iTauMax+1; iTime<nTau; iTime++){
      if(exportIF){
	for(int iSample=0; iSample<nSample ; iSample++){
	  IFcumhazard.slice(iSample).col(iTime).fill(NA_REAL);
	}
      }
      if(exportSE){
	SEcumhazard.col(iTime).fill(NA_REAL);
      }
      if(exportIFsum_cumhazard){
	IFsum_cumhazard.col(iTime).fill(NA_REAL);
      }
      if(exportIFsum_survival){
	IFsum_survival.col(iTime).fill(NA_REAL);
      }
    }
  }
 
  // SEcumhazard
  return(List::create(Named("iid") = IFcumhazard,
		      Named("iidsum_cumhazard") = IFsum_cumhazard,
		      Named("iidsum_survival") = IFsum_survival,
                      Named("se") = SEcumhazard));
}
// }}}

// {{{ calcSeCif_cpp
// [[Rcpp::export]]
List calcSeCif_cpp(const NumericVector& seqTau, // horizon time for the predictions
		   const NumericVector& jumpTime, 
		   const LogicalVector& jumpTheCause, 
		   const arma::mat& indexJump,
		   const arma::mat& indexSample,
		   const std::vector< arma::mat >& IFbeta,
                   const arma::mat& cif,
		   const std::vector< arma::mat >& Ehazard0,
		   const std::vector< arma::mat >& cumEhazard0,
		   const arma::mat& iS0,
		   const std::vector< NumericVector >& cumhazard_iS0,
   		   const std::vector< NumericVector >& hazard_iS0,
		   const arma::mat& newEXb,
                   const arma::mat& sampleEXb,
		   const std::vector< arma::mat >& X,
		   const arma::mat sameStrata,
		   const NumericVector& sampleTime,
  		   const std::vector< NumericVector>& hazard0,
		   const std::vector< NumericVector>& cumhazard0,
		   int theCause, double firstJumpTime, double lastSampleTime,
 		   int nTau, int nJump, int nNewObs, int nSample, int nCause, const IntegerVector& p,
		   bool survtype,
 		   bool exportSE, bool exportIF, bool exportIFsum){

  // // define objects
  NumericVector veci_IF_risk(nSample);
  double i_IF_cumhazard; // cumul for all causes
  double i_IF_hazard; // for the cause of interest
  double i_survival; // for all causes
  double i_hazard; // for the cause of interest
  
  int iTau0=0,iTauMax=nTau-1;
  int indexTempo;
  double IF_risk_tempo; // for the export
  
  arma::mat SErisk;
  if(exportSE){
    SErisk.resize(nNewObs, nTau);
    SErisk.fill(0);
  }  
  arma::cube IFrisk;
  if(exportIF){
    IFrisk.resize(nNewObs, nTau, nSample);
    IFrisk.fill(0);
  }
  arma::mat IFsumrisk;
  if(exportIFsum){
    IFsumrisk.resize(nSample, nTau);
    IFsumrisk.fill(0);
  }

  // narrow prediction times
  while((iTau0 < nTau) && firstJumpTime>seqTau[iTau0]){ // start at the first event or after 
    iTau0++;    
  }
  while((iTauMax >= 0) && seqTau[iTauMax]>lastSampleTime){ // end at the last event or before
    iTauMax--;
  }
  int iTau;

   // interesting times
  if(iTau0 < nTau && iTauMax >= 0){
	
    for(int iNewObs=0; iNewObs<nNewObs ; iNewObs++){

      R_CheckUserInterrupt();    

      i_survival = 1;
      iTau = iTau0; 
      std::fill(veci_IF_risk.begin(), veci_IF_risk.end(), 0);
     	
      for(int iJump=0; iJump<nJump; iJump++){
	if(jumpTime[iJump]>lastSampleTime){break;}

	//// compute the hazard for the cause of interest ////
	i_hazard = hazard0[theCause][iJump]*newEXb(iNewObs,theCause);

	//// compute the survival ////
	// it is survival at t- which is computed i.e. the survival at the previous eventtime (censoring does not affect survival)
	if(iJump>0){
	  if(survtype){ // product limit
	    i_survival *= exp(-cumhazard0[1](iJump-1)*newEXb(iNewObs,1));
	  }else{	    
	    i_survival = 0; 
	    for(int iterC=0 ; iterC<nCause; iterC++){
	      i_survival += cumhazard0[iterC](iJump-1)*newEXb(iNewObs,iterC);
	    }
	    i_survival = exp(-i_survival);
	  }
	}// otherwise the survival stays at 1

	for(int iSample=0; iSample<nSample ; iSample++){

	  indexTempo = std::min(indexSample(iSample,theCause),indexJump(iJump,theCause));
	  // compute the influence function for the hazard
  	  i_IF_hazard = calcIFhazard(jumpTime[iJump],
				     sampleTime[iSample],
  				     IFbeta[theCause].row(iSample),
  				     Ehazard0[theCause].col(iJump),
				     X[theCause].row(iNewObs),
				     hazard_iS0[theCause][indexTempo],
				     newEXb(iNewObs,theCause),
				     sampleEXb(iSample,theCause),  
				     hazard0[theCause][iJump],
				     iS0(iSample,theCause),
				     sameStrata(iSample,theCause),
				     jumpTheCause[iJump], 
				     p[theCause], true);

	  // cumulate the influence function for the cumulative hazard
	  // it is at t-
	  i_IF_cumhazard = 0;
	  if(iJump>0){
	  for(int iCause=0; iCause<nCause; iCause++){
	    indexTempo = std::min(indexSample(iSample,iCause),indexJump(iJump-1,iCause));

	    i_IF_cumhazard += calcIFhazard(jumpTime[iJump-1],
					   sampleTime[iSample],
					   IFbeta[iCause].row(iSample),
					   cumEhazard0[iCause].col(iJump-1),
					   X[iCause].row(iNewObs),
					   cumhazard_iS0[iCause][indexTempo],
					   newEXb(iNewObs,iCause),
					   sampleEXb(iSample,iCause),  
					   cumhazard0[iCause][iJump-1],
					   iS0(iSample,iCause),
					   sameStrata(iSample,iCause),
					   true, // this argument is ignored when computing the influence function for the cumulative hazard
					   p[iCause], false);
	    
	  }
	  }

	  // update the influence function for the absolute risk
	  veci_IF_risk[iSample] += i_survival *(i_IF_hazard - i_hazard * i_IF_cumhazard);


	} // end loop over sample

	//// export results ////
	// while there are remaining times to export
	// AND the prediction time is before the next event
	//     OR its the last jump, i.e. the prediction time coincide with the last jump (otherwise it would have been removed at the begining)

	while((iTau <= iTauMax) && ( (((iJump+1)<nJump) && (seqTau[iTau] < jumpTime[iJump+1])) || (iJump+1==nJump))){

	  for(int iSample=0; iSample<nSample ; iSample++){
            IF_risk_tempo = veci_IF_risk[iSample];
	    
	    if(exportSE){
	      SErisk(iNewObs,iTau) += pow(IF_risk_tempo,2);
	    }
	    if(exportIF){
	      IFrisk(iNewObs,iTau,iSample) = IF_risk_tempo;
	    }
	    if(exportIFsum){
	      IFsumrisk(iSample,iTau) += IF_risk_tempo;
	    }
	  }
	  iTau++;	  
	}      
    } // end loop over jump times
	  
      // set IF/variance to NA after the last event time
      if(iTau < nTau){
        for(int iTime=iTau; iTime<nTau; iTime++){
	  if(exportIF){
	    for(int iSample=0; iSample<nSample ; iSample++){
	      IFrisk(iNewObs,iTime,iSample) = NA_REAL;
	    }
	  }
	  if(exportSE){
	    SErisk(iNewObs,iTime) = NA_REAL;
	  }	  
	}
      }
    } // end loop over new observations

    // return standard error instead of variance
    if(exportSE){
      for(int iTime=iTau0; iTime<nTau; iTime++){
	SErisk.col(iTime) = sqrt(SErisk.col(iTime));
      }
    }    
  }

   // set IF/variance to NA after the last event time
  if(iTauMax+1 < nTau){
    for(int iTime=iTauMax+1; iTime<nTau; iTime++){
      if(exportIF){
	for(int iSample=0; iSample<nSample ; iSample++){
	  IFrisk.slice(iSample).col(iTime).fill(NA_REAL);
	}
      }
      if(exportSE){
	SErisk.col(iTime).fill(NA_REAL);
      }
      if(exportIFsum){
	IFsumrisk.col(iTime).fill(NA_REAL);
      }
    }
  }
  
  // export
  return(List::create(Named("iid") = IFrisk,
		      Named("iidsum") = IFsumrisk,
		      Named("se") = SErisk));
}
// }}}

// {{{ calcIFhazard
inline double calcIFhazard(double time,
			   double sampleTime,			   
			   const rowvec& IFbeta,
			   const colvec& Ehazard0,
			   const rowvec& X,
			   double hazard_iS0,
			   double newEXb,
			   double sampleEXb,
			   double hazard0,
			   double iS0,
			   bool sameStrata,
			   bool jumpTime,
			   int p, bool hazard){

  double IFbetaE=0;
  double XIFbeta=0;
  double IF_hazard0=0;
  double IF_hazard;

  if(p>0){
    for(int iX = 0; iX < p; iX++){     
      IFbetaE += IFbeta[iX] * Ehazard0[iX];
      XIFbeta += X[iX] * IFbeta[iX];
    }
  }

  if(sameStrata){
    if(hazard==true){
      IF_hazard0 = - IFbetaE;      
      if(sampleTime==time){
        IF_hazard0 += iS0;
      }
      if(jumpTime && time <= sampleTime){
        IF_hazard0 -= sampleEXb * hazard_iS0;
      }
    }else{ // cumulative hazard
      if(sampleTime<=time){
       IF_hazard0 = - IFbetaE - sampleEXb * hazard_iS0 + iS0;
      }else{
       IF_hazard0 = - IFbetaE - sampleEXb * hazard_iS0;
      }
    }
  }else{
    IF_hazard0 = - IFbetaE;
  }
   
  if(p>0){	    
    IF_hazard = newEXb*(IF_hazard0 + hazard0 * XIFbeta);	  
  }else{
    IF_hazard = IF_hazard0;
  }

  return(IF_hazard);    
}
// }}}
