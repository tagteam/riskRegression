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
  
// * calcSeHazard_cpp: compute IF for survival (method 1)
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

// * calcAIFsurv_cpp: compute average IF for survival (method 3)
// [[Rcpp::export]]
std::vector< std::vector<arma::mat> > calcAIFsurv_cpp(const std::vector<arma::mat>& ls_IFcumhazard,
				       const arma::mat& IFbeta,
						      const std::vector<arma::rowvec>& cumhazard0,
				       const arma::mat& survival,
				       const arma::colvec& eXb,
				       const arma::mat& X,
				       const NumericVector& prevStrata,
						      const std::vector<arma::uvec>& ls_indexStrata,
 				       const std::vector<arma::mat>& factor,
				       int nTimes,
				       int nObs,
				       int nStrata,
				       int nVar,
				       bool exportCumHazard,
				       bool exportSurvival){

  // ** prepare output
  int nFactor = factor.size();
  std::vector< std::vector<arma::mat> > out(2);
  if(exportCumHazard){
    out[0].resize(nFactor);
    for(int iFactor=0; iFactor<nFactor; iFactor++){
    out[0][iFactor].resize(nObs, nTimes);
    out[0][iFactor].fill(0.0);
    }
  }
  if(exportSurvival){
        out[1].resize(nFactor);
    for(int iFactor=0; iFactor<nFactor; iFactor++){
    out[1][iFactor].resize(nObs, nTimes);
    out[1][iFactor].fill(0.0);
    }
  }

  // ** initialize
  arma::rowvec rowvec_tempo;
  arma::mat IFtempo1;
  arma::mat IFtempo2;
  int nObs_strata;
  arma::mat factor_tempo;
  
  if(nVar == 0){
    // ** non-parameteric
    // <IF>(cumhazard) = IF_cumhazard0
    // <IF>(Surv) = E[Surv] IF_cumhazard0

    for(int iStrata=0; iStrata<nStrata; iStrata++){

      nObs_strata = ls_indexStrata[iStrata].size();
  
       for(int iFactor=0; iFactor<nFactor; iFactor++){

	 factor_tempo = factor[iFactor].rows(ls_indexStrata[iStrata]);
	 
         if(exportCumHazard){
   	   IFtempo1 = ls_IFcumhazard[iStrata];
	   // compute average factor (=1 by default)
	   rowvec_tempo = sum(factor_tempo, 0)/nObs_strata; // colSums
           // multiply by prevalence
	   IFtempo1.each_row() %= (prevStrata[iStrata] * rowvec_tempo);
	   // update cumhazard
	   out[0][iFactor] += IFtempo1;
      }
      if(exportSurvival){
 	   IFtempo1 = ls_IFcumhazard[iStrata];	      
           // compute average survival times factor
           rowvec_tempo = sum(survival.rows(ls_indexStrata[iStrata]) % factor_tempo,0)/nObs_strata; // colSums
          // multiply by prevalence
 	  IFtempo1.each_row() %= (prevStrata[iStrata] * rowvec_tempo);
	  // update survival
	  out[1][iFactor] -= IFtempo1;
      }
      }
      
    }

    
  }else{
    // ** semi-parametric
    // <IF>(cumhazard) = E[eXb] IF_cumhazard0 + E[eXb * cumhazard0 * X] IF_beta
    // <IF>(Surv) = E[Surv * eXb] IF_cumhazard0 + E[Surv * eXb * cumhazard0 * X] IF_beta
    int iObsStrata;
    arma::rowvec E_eXb;
    arma::rowvec E_S_eXb;
    arma::mat E_eXb_cumhazard0_X;
    arma::mat E_S_eXb_cumhazard0_X;

    arma::rowvec X_eXb;
    arma::rowvec S_eXb;
    
    if(exportCumHazard){
      E_eXb.resize(nTimes);     
      E_eXb_cumhazard0_X.resize(nVar, nTimes);      
    }
	
    if(exportSurvival){
      E_S_eXb.resize(nTimes);      
      E_S_eXb_cumhazard0_X.resize(nVar, nTimes);
    }

    // for(int iStrata=0; iStrata<nStrata; iStrata++){
    for(int iStrata=0; iStrata<nStrata; iStrata++){
	  
      nObs_strata = ls_indexStrata[iStrata].size();

     for(int iFactor=0; iFactor<nFactor; iFactor++){

         factor_tempo = factor[iFactor].rows(ls_indexStrata[iStrata]);
         factor_tempo.each_col() %= eXb.rows(ls_indexStrata[iStrata]);
  
      if(exportCumHazard){	
	E_eXb.fill(0.0);
	E_eXb_cumhazard0_X.fill(0.0);

	// compute n*expectation
	for(int iObs=0; iObs<nObs_strata; iObs++){
	  iObsStrata = ls_indexStrata[iStrata][iObs];
		    
	  E_eXb += factor_tempo.row(iObs); // note: it is iObs and not iObs_strata

	  E_eXb_cumhazard0_X += (X.row(iObsStrata)).t() * (cumhazard0[iStrata] % factor_tempo.row(iObs));
	}
	// divide by n (for the average) and multiply first term by the prevalence
	E_eXb *= prevStrata[iStrata] / nObs_strata;
	E_eXb_cumhazard0_X *= prevStrata[iStrata] / nObs_strata;

	// update first term
	IFtempo1 = ls_IFcumhazard[iStrata];
	IFtempo1.each_row() %= E_eXb;
	
	// export
	out[0][iFactor] += IFtempo1 + IFbeta * E_eXb_cumhazard0_X;		  
      }

      if(exportSurvival){
	E_S_eXb.fill(0.0);
	E_S_eXb_cumhazard0_X.fill(0.0);
      
	// compute n*expectation
	for(int iObs=0; iObs<nObs_strata; iObs++){
	  iObsStrata = ls_indexStrata[iStrata][iObs];

	  E_S_eXb += survival.row(iObsStrata) % factor_tempo.row(iObs);

	  E_S_eXb_cumhazard0_X += (X.row(iObsStrata)).t() * (cumhazard0[iStrata] % survival.row(iObsStrata) % factor_tempo.row(iObs));
	}
	// divide by n (for the average) and multiply first term by the prevalence
	E_S_eXb *= prevStrata[iStrata] / nObs_strata;
	E_S_eXb_cumhazard0_X *= prevStrata[iStrata] / nObs_strata;

	// update first term
	IFtempo1 = ls_IFcumhazard[iStrata];
	IFtempo1.each_row() %= E_S_eXb;

	// export		 
	out[1][iFactor] -= IFtempo1 + IFbeta * E_S_eXb_cumhazard0_X;
        }
      }
    }
  }

  return(out);
  
}


// * calcSeCif_cpp
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

// * calcIFhazard
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

// * calcSeCif2_cpp
// [[Rcpp::export]]
List calcSeCif2_cpp(std::vector<arma::mat> ls_IFbeta, std::vector<arma::mat> ls_X,
		    std::vector<arma::mat> ls_cumhazard, std::vector<arma::mat> ls_hazard,
		    std::vector< std::vector<arma::mat> > ls_IFcumhazard, std::vector< std::vector<arma::mat> > ls_IFhazard,
		    NumericMatrix eXb,
		    arma::vec timeIndex, int nIndex, arma::vec time,
		    int nObs, int nJumpTime, arma::vec maxJumpTime,
		    int theCause, int nCause, bool hazardType, arma::vec nVar,
		    int nNewObs, NumericMatrix strata,
		    bool exportSE, bool exportIF, bool exportIFsum){

   arma::mat X_IFbeta;

   arma::mat IFhazard;
   arma::mat IFcumhazard;
   arma::colvec hazard;
   arma::colvec cumhazard;

   double ieXb;
   int iStrataCause;

   arma::colvec IF_tempo;
   arma::colvec cumIF_tempo;

   // ** initialize for export
   arma::mat outSE;
   if(exportSE){
     outSE.resize(nNewObs,nIndex);
     outSE.fill(NA_REAL);
   }
   arma::cube outIF;
   if(exportIF){
     outIF.resize(nNewObs,nObs,nIndex);
     outIF.fill(NA_REAL);
   }
   arma::mat outIFsum;
   if(exportIFsum){
     outIFsum.resize(nObs,nIndex);
     outIFsum.fill(NA_REAL);
   }

   // ** prepare arguments
   std::vector<arma::mat> ls_tcumhazard(nCause);
   std::vector<arma::mat> ls_thazard(nCause);
   for(int iCause=0; iCause<nCause; iCause ++){
           if(nVar[iCause] > 0){
	     ls_thazard[iCause] = ls_hazard[iCause].t();
             ls_tcumhazard[iCause] = ls_cumhazard[iCause].t();
	   }
   }
   
   // ** take care of the prediction times before the first jump
   int iTau = 0;
   // Rcout << "1";
   while( (timeIndex[iTau] < 0) && (iTau < nIndex) ){
     if(exportSE){
       outSE.col(iTau) = zeros<colvec>(nNewObs);
     }
     if(exportIF){
       outIF.slice(iTau) = zeros<mat>(nNewObs,nObs);
     }
     if(exportIFsum){
       outIFsum.col(iTau) = zeros<colvec>(nNewObs);
     }
     iTau++;
   }
   // Rcout << "-end ";
   
   // ** prepare the influence function
   int iiTau;
   
   for(int iNewObs=0; iNewObs<nNewObs; iNewObs++){
     iiTau = iTau;
     IFcumhazard = zeros<mat>(nObs,nJumpTime);
     cumhazard = zeros<colvec>(nJumpTime);
     ieXb = NA_REAL;
   
   for(int iCause=0; iCause<nCause; iCause ++){

     iStrataCause = strata(iNewObs,iCause);

        // Rcout << "2 ";
       if(nVar[iCause]>0){
	 X_IFbeta = ls_IFbeta[iCause] * (ls_X[iCause].row(iNewObs)).t();
         ieXb = eXb(iNewObs,iCause);
       }

               // Rcout << "3 ";
       if(hazardType || (iCause != theCause)){
	   if(nVar[iCause]>0){
	        cumhazard += ieXb * ls_cumhazard[iCause].col(iStrataCause);
	        // IFcumhazard += ieXb * (ls_IFcumhazard[iCause][iStrataCause] + X_IFbeta * ls_cumhazard[iCause].col(iStrataCause).t());
	        IFcumhazard += ieXb * (ls_IFcumhazard[iCause][iStrataCause] + X_IFbeta * ls_tcumhazard[iCause].row(iStrataCause));
	   }else{
		  cumhazard += ls_cumhazard[iCause].col(iStrataCause);
		  IFcumhazard += ls_IFcumhazard[iCause][iStrataCause];
	   }
       }

               // Rcout << "4 ";
       if(iCause == theCause){
           if(nVar[iCause] > 0){
               hazard = ieXb * ls_hazard[iCause].col(iStrataCause);
	       // IFhazard = ieXb * (ls_IFhazard[iCause][iStrataCause] + X_IFbeta * ls_hazard[iCause].col(iStrataCause).t());
	       IFhazard = ieXb * (ls_IFhazard[iCause][iStrataCause] + X_IFbeta * ls_thazard[iCause].row(iStrataCause));
           } else{
	       hazard = ls_hazard[iCause].col(iStrataCause);
	       IFhazard = ls_IFhazard[iCause][iStrataCause];
           }
       }
   }
     
   // ** loop over time
   cumIF_tempo = zeros<colvec>(nObs);   

   for(int iTime=0; iTime<nJumpTime; iTime++){
               // Rcout << "5 " ;
     // prepare IF
     if(hazard[iTime]>0){
     if(iTime==0){
       IF_tempo = IFhazard.col(iTime);
     }else{
       // cumhazard is evaluated just before the jump
       IF_tempo = (IFhazard.col(iTime) - IFcumhazard.col(iTime-1) * hazard[iTime]) * exp(-cumhazard[iTime-1]);
     }
      cumIF_tempo = cumIF_tempo + IF_tempo;
     }

                     // Rcout << "6";
      // store
   while((iTime == timeIndex[iiTau]) && (time[iiTau] <= maxJumpTime[iNewObs])){
                          // Rcout << "*";
     if(exportSE){
       // Rcout << "a";
       outSE.row(iNewObs).col(iiTau) = sqrt(accu(pow(cumIF_tempo,2)));
     }
     if(exportIF){
       // Rcout << "b";
       outIF.slice(iiTau).row(iNewObs) = cumIF_tempo.t();
     }
     if(exportIFsum){
       // Rcout << "c";
       if(iNewObs==0){
         outIFsum.col(iiTau) = cumIF_tempo/nNewObs; // necessary because outIFsum is initialized as NA
       }else{
	 outIFsum.col(iiTau) += cumIF_tempo/nNewObs;
       }
     }
     iiTau++;
     if(iiTau == nIndex){break;}
   }
   // Rcout << "-end " << endl;
   if(iiTau == nIndex){break;} // not necessary because nJumpTime should ensure that it is never triggered - so just for safety

   }
   // Rcout << "endend" << endl;	
   }
   
   // ** export
  return(List::create(Named("se") = outSE,
		      Named("iid") = outIF,
		      Named("average.iid") = outIFsum));

}

