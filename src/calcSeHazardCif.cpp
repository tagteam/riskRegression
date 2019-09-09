// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
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

// * calcSeCox
// ** calcSeHazard_cpp: compute IF for the cumlative hazard / survival (method 1)
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

// ** calcAIFsurv_cpp: compute average IF for the cumlative hazard / survival (method 3)
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
													  bool diag,
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
			
		  if(diag){
			// apply individual weights (a weight for each individual, i.e. each time)
			IFtempo1.each_row() %= trans(factor_tempo);
			// average and multiply by prevalence
 			IFtempo1 = sum(IFtempo1,1) * prevStrata[iStrata]/nObs_strata;
		  }else{		
			// compute average factor (=1 by default)
			rowvec_tempo = sum(factor_tempo, 0)/nObs_strata; // colSums
			// multiply by prevalence
			IFtempo1.each_row() %= (prevStrata[iStrata] * rowvec_tempo);
		  }
		  // update cumhazard
		  out[0][iFactor] += IFtempo1;
		}
		if(exportSurvival){
		  IFtempo1 = ls_IFcumhazard[iStrata];	      
		  if(diag){
			// apply individual weights (a weight for each individual, i.e. each time)
			IFtempo1.each_row() %= trans(survival.rows(ls_indexStrata[iStrata]) % factor_tempo);
			// average and multiply by prevalence
 			IFtempo1 = sum(IFtempo1,1) * prevStrata[iStrata]/nObs_strata;
		  }else{
			// compute average survival times factor
			rowvec_tempo = sum(survival.rows(ls_indexStrata[iStrata]) % factor_tempo,0)/nObs_strata; // colSums
			// multiply by prevalence
			IFtempo1.each_row() %= (prevStrata[iStrata] * rowvec_tempo);
		  }
		  // update survival
		  out[1][iFactor] -= IFtempo1;
		}
      }
      
    }

    
  }else{
    // ** semi-parametric

    int iObsStrata;
	arma::rowvec E_eXb; // when diag=FALSE
	arma::colvec E2_eXb; // when diag=TRUE
    arma::rowvec E_S_eXb; // when diag=FALSE
    arma::colvec E2_S_eXb; // when diag=TRUE
    arma::mat E_eXb_cumhazard0_X;
    arma::mat E_S_eXb_cumhazard0_X;

    arma::rowvec X_eXb;
    arma::rowvec S_eXb;
    
    if(exportCumHazard){
	  if(diag){
		// <IF>(cumhazard) = E[eXb IF_cumhazard0] + E[eXb * cumhazard0 * X] IF_beta
		E2_eXb.resize(nObs);     
	  }else{
		// <IF>(cumhazard) = E[eXb] IF_cumhazard0 + E[eXb * X] * cumhazard0 * IF_beta
		E_eXb.resize(nTimes);     
	  }
		E_eXb_cumhazard0_X.resize(nVar, nTimes);
    }
	
    if(exportSurvival){
	  if(diag){
		E2_S_eXb.resize(nObs);
		// <IF>(Surv) = E[Surv * eXb IF_cumhazard0] + E[Surv * eXb * cumhazard0 * X] IF_beta
	  }else{
		E_S_eXb.resize(nTimes);
	  // <IF>(Surv) = E[Surv * eXb] IF_cumhazard0 + E[Surv * eXb * X] * cumhazard0 IF_beta
	  }
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

		  if(diag){
			
			// compute n*expectation
			for(int iObs=0; iObs<nObs_strata; iObs++){
			  iObsStrata = ls_indexStrata[iStrata][iObs];
			  E2_eXb += ls_IFcumhazard[iStrata].col(iObs) * factor_tempo(iObs,0); // note: it is iObs and not iObs_strata
			  E_eXb_cumhazard0_X += trans(X.row(iObsStrata)) * (cumhazard0[iStrata](iObs) * factor_tempo(iObs,0));
			}
			// divide by n (for the average) and multiply first term by the prevalence
			E2_eXb *= prevStrata[iStrata] / nObs_strata;
			E_eXb_cumhazard0_X *= prevStrata[iStrata] / nObs_strata;
			
			// export							 
			out[0][iFactor] += E2_eXb + IFbeta * E_eXb_cumhazard0_X;		  

		  }else{

			// compute n*expectation
			for(int iObs=0; iObs<nObs_strata; iObs++){
			  iObsStrata = ls_indexStrata[iStrata][iObs];
			  E_eXb += factor_tempo.row(iObs); // note: it is iObs and not iObs_strata
			  E_eXb_cumhazard0_X += trans(X.row(iObsStrata)) * factor_tempo.row(iObs);
			}
			// divide by n (for the average) and multiply first term by the prevalence
			E_eXb *= prevStrata[iStrata] / nObs_strata;
			E_eXb_cumhazard0_X.each_row() %= cumhazard0[iStrata] * (prevStrata[iStrata] / nObs_strata);


			// update first term
			IFtempo1 = ls_IFcumhazard[iStrata];
			IFtempo1.each_row() %= E_eXb;
	
			// export
			out[0][iFactor] += IFtempo1 + IFbeta * E_eXb_cumhazard0_X;		  
		  }
		}

		if(exportSurvival){
		  E_S_eXb.fill(0.0);
		  E_S_eXb_cumhazard0_X.fill(0.0);
      
		  if(diag){
			// compute n*expectation
			for(int iObs=0; iObs<nObs_strata; iObs++){
			  iObsStrata = ls_indexStrata[iStrata][iObs];
			  E2_S_eXb += ls_IFcumhazard[iStrata].col(iObs) * (survival(iObsStrata,0) * factor_tempo(iObs,0));
			  E_S_eXb_cumhazard0_X += trans(X.row(iObsStrata)) * (cumhazard0[iStrata](iObs) * survival(iObsStrata,0) * factor_tempo(iObs,0));
			}
			// divide by n (for the average) and multiply first term by the prevalence
			E2_S_eXb *= prevStrata[iStrata] / nObs_strata;
			E_S_eXb_cumhazard0_X *= (prevStrata[iStrata] / nObs_strata);

			// export							 
			out[1][iFactor] -= E2_S_eXb + IFbeta * E_S_eXb_cumhazard0_X;		  
		  }else{
			// compute n*expectation
			for(int iObs=0; iObs<nObs_strata; iObs++){
			  iObsStrata = ls_indexStrata[iStrata][iObs];

			  E_S_eXb += survival.row(iObsStrata) % factor_tempo.row(iObs);

			  E_S_eXb_cumhazard0_X += trans(X.row(iObsStrata)) * (survival.row(iObsStrata) % factor_tempo.row(iObs));
			}
			// divide by n (for the average) and multiply first term by the prevalence
			E_S_eXb *= prevStrata[iStrata] / nObs_strata;
			E_S_eXb_cumhazard0_X.each_row() %= cumhazard0[iStrata] * (prevStrata[iStrata] / nObs_strata);

			// update first term
			IFtempo1 = ls_IFcumhazard[iStrata];
			IFtempo1.each_row() %= E_S_eXb;

			// export		 
			out[1][iFactor] -= IFtempo1 + IFbeta * E_S_eXb_cumhazard0_X;
		  }
        }
      }
    }
  }

  return(out);
  
}


// * calcSeCSC
// ** calcSeCif_cpp: compute IF for the absolute risk (method 1)
// [[Rcpp::export]]
List calcSeCif_cpp(const NumericVector& seqTau, // horizon time for the predictions
				   const NumericVector& jumpTime, 
				   const LogicalVector& jumpTheCause, 
				   const arma::mat& indexJump,
				   const arma::mat& indexSample,
				   const std::vector< arma::mat >& IFbeta,
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


// ** calcSeCif2_cpp: compute IF for the absolute risk (method 2)
// [[Rcpp::export]]
List calcSeCif2_cpp(std::vector<arma::mat> ls_IFbeta, std::vector<arma::mat> ls_X,
					std::vector<arma::mat> ls_cumhazard, arma::mat ls_hazard,
					std::vector< std::vector<arma::mat> > ls_IFcumhazard, std::vector<arma::mat> ls_IFhazard,
					NumericMatrix eXb,
					int nJumpTime, NumericVector JumpMax,
					NumericVector tau, arma::vec tauIndex, int nTau,
					int nObs,  
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
    outSE.resize(nNewObs,nTau);
    outSE.fill(0);
  }
  arma::cube outIF;
  if(exportIF){
    outIF.resize(nNewObs,nObs,nTau);
    outIF.fill(0);
  }
  arma::mat outIFsum;
  if(exportIFsum){
    outIFsum.resize(nObs,nTau);
    outIFsum.fill(0);
  }

  // ** prepare arguments
  std::vector<arma::mat> ls_tcumhazard(nCause);
  arma::mat ls_thazard;
  if(nVar[theCause] > 0){
    ls_thazard = ls_hazard.t();
  }
  for(int iCause=0; iCause<nCause; iCause ++){
    if(nVar[iCause] > 0){
      ls_tcumhazard[iCause] = ls_cumhazard[iCause].t();
    }
  }
   
  // ** skip before first event (out is initialized at 0)
  int iTau = 0;
  while(iTau < nTau && tauIndex[iTau] < 0){
    iTau++;
  }
  int iiTau;

  // ** prepare the influence function
  for(int iNewObs=0; iNewObs<nNewObs; iNewObs++){

    R_CheckUserInterrupt();
     
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
		  hazard = ieXb * ls_hazard.col(iStrataCause);
		  // IFhazard = ieXb * (ls_IFhazard[iCause][iStrataCause] + X_IFbeta * ls_hazard[iCause].col(iStrataCause).t());
		  IFhazard = ieXb * (ls_IFhazard[iStrataCause] + X_IFbeta * ls_thazard.row(iStrataCause));
		} else{
		  hazard = ls_hazard.col(iStrataCause);
		  IFhazard = ls_IFhazard[iStrataCause];
		}
      }
    }
     
    // ** loop over time
    cumIF_tempo = zeros<colvec>(nObs);   

    for(int iJump=0; iJump<nJumpTime; iJump++){

      // Rcout << "5 " ;
      // prepare IF
      if(hazard[iJump]>0){
		if(iJump==0){
		  IF_tempo = IFhazard.col(iJump);
		}else{
		  // cumhazard is evaluated just before the jump
		  IF_tempo = (IFhazard.col(iJump) - IFcumhazard.col(iJump-1) * hazard[iJump]) * exp(-cumhazard[iJump-1]);
		}
		cumIF_tempo = cumIF_tempo + IF_tempo;
      }

      // Rcout << "6";
      // store
      while(iiTau < nTau && tauIndex[iiTau] == iJump && tau[iiTau] <= JumpMax[iNewObs]){
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
		  outIFsum.col(iiTau) += cumIF_tempo/nNewObs;
		}
		iiTau++;
      }
      // Rcout << "-end " << endl;
      if(iiTau == nTau){break;} // not necessary because nJumpTime should ensure that it is never useful - so just for safety

    }

    // ** fill remaining columns with NA
    while(iiTau < nTau){
      if(exportSE){
		outSE.row(iNewObs).col(iiTau).fill(NA_REAL);
      }
      if(exportIF){
		// Rcout << "b";
		outIF.slice(iiTau).row(iNewObs).fill(NA_REAL);
      }
      if(exportIFsum){
		outIFsum.col(iiTau).fill(NA_REAL);
      }
      iiTau++;
    }
	 
    // Rcout << "endend" << endl;	
  }
   
  // ** export
  return(List::create(Named("se") = outSE,
					  Named("iid") = outIF,
					  Named("average.iid") = outIFsum));

}

// ** calcAIFcif_cpp: compute average IF for the absolute risk (method 3)
// [[Rcpp::export]]
std::vector<arma::mat> calcAIFcif_cpp(const arma::mat& hazard1,
									  const std::vector<arma::mat>& ls_cumhazard,
									  const std::vector<arma::mat>& ls_tX,
									  const arma::mat& eXb,
									  const std::vector<arma::mat>& ls_IFbeta,
									  const std::vector<arma::mat>& ls_IFhazard,
									  const std::vector< std::vector<arma::mat> >& ls_IFcumhazard,
									  int nCause, int theCause, bool hazardType,
									  const NumericVector& tau, int nTau, const IntegerVector& tauIndex, 
									  int nJumpTime, const NumericVector& JumpMax,  
									  int nObs,
									  int nNewObs,
									  const arma::mat& levelStrata, int nStrata, const std::vector<IntegerVector>& ls_indexStrata,
									  const IntegerVector& nVar,
									  const arma::mat& factor){
  // note: tauIndex must be -1 before the first event and greater than nTau after the last event
  //       ls_tX must be a list of design matrix where the observations are in columns
  // Rcout << "Start!" << endl;
	
  // ** prepare output
  int nFactor = factor.n_cols;

  std::vector<arma::mat> out(nFactor);
  for(int iFactor=0; iFactor<nFactor; iFactor++){
    out[iFactor].resize(nObs, nTau);
    out[iFactor].fill(0.0);
  }

  // ** intialization
  arma::rowvec Esurv_eXb1(nFactor);
  arma::mat Esurv_hazard1_eXb(nCause,nFactor);
  std::vector<arma::mat> Esurv_hazard1_X_cumhazard(nCause);

  std::vector<arma::mat> int_bmissingIF(nCause);
  arma::mat int_hazard1(nObs,nFactor);
  arma::mat int_cumhazard(nObs,nFactor);
  arma::mat int_b(nObs,nFactor);
  arma::mat int_total(nObs,nFactor);
  
  for(int iCause=0; iCause<nCause; iCause++){
    if(nVar[iCause]>0){
      Esurv_hazard1_X_cumhazard[iCause].resize(nVar[iCause],nFactor);
      int_bmissingIF[iCause].resize(nVar[iCause],nFactor);
    }
  }

  double survivalTempo;
  NumericVector cumhazardTempo(nCause);
  double hazard1Tempo;
  double factorTempo;
  int iObs;
  int nObs_strata;
  															 
  // ** skip before first event (out is initialized at 0)
  int iTau = 0;
  while(iTau < nTau && tauIndex[iTau] < 0){
    iTau ++;
  }
  int iiTau;
  
  // ** loop
  for(int iStrata=0; iStrata<nStrata; iStrata++){
    // Rcout << "Strata:" << iStrata << endl << endl;
    nObs_strata = ls_indexStrata[iStrata].size();

    // start from iTau
    iiTau = iTau;
     
    // initialize integral terms to 0
    int_hazard1.fill(0.0);
    int_cumhazard.fill(0.0);
    for(int iCause=0; iCause<nCause; iCause++){
      if(nVar[iCause]>0){
		int_bmissingIF[iCause].fill(0.0);
      }
    }
	 
    for(int iJump=0; iJump<nJumpTime; iJump++){
      // Rcout << "Jump " << iJump << ": ";

      R_CheckUserInterrupt();
       
      // after last jump or last prediciton time stop computing the integral
      if(iiTau == nTau || tau[iiTau] > JumpMax[iStrata]){break;}
	   
      // initialize expectation
      Esurv_eXb1.fill(0.0);
      Esurv_hazard1_eXb.fill(0.0);
      for(int iCause=0; iCause<nCause; iCause++){
		if(nVar[iCause]>0){
		  Esurv_hazard1_X_cumhazard[iCause].fill(0.0);
		}
      }
	   
      // *** compute expectation
      // Rcout << "(expectation ";	   
      for(int iObsStrata=0; iObsStrata<nObs_strata; iObsStrata++){
		iObs = ls_indexStrata[iStrata][iObsStrata];

		if(hazard1(iJump,levelStrata(iStrata,theCause))>0){
		  // Rcout << iJump << " " << levelStrata(iStrata,theCause) << " " << iObs << " " << theCause;	   

		  hazard1Tempo = hazard1(iJump,levelStrata(iStrata,theCause)) * eXb(iObs,theCause);
		  survivalTempo = 0.0;
		   
		  if(iJump>0){
			if(hazardType){			 
			  for(int iCause=0; iCause<nCause; iCause++){
				cumhazardTempo[iCause] = ls_cumhazard[iCause](iJump-1,levelStrata(iStrata,iCause)) * eXb(iObs,iCause);
				survivalTempo += cumhazardTempo[iCause];
			  }
			}else{
			  cumhazardTempo[1] = ls_cumhazard[1](iJump-1,levelStrata(iStrata,1)) * eXb(iObs,1);
			  survivalTempo += cumhazardTempo[1];
			}
		  }else{
			cumhazardTempo.fill(0.0);
			// for(int iCause=0; iCause<nCause; iCause++){
			// cumhazardTempo[iCause] = 0.0; // used instead of .fill(0.0) because it was not doing the job
			// }
		  }
		 
		  survivalTempo = exp(-survivalTempo);
		  // Rcout << " | survival: " << survivalTempo;	
		  for(int iFactor=0; iFactor<nFactor; iFactor++){
			factorTempo = factor(iObs,iFactor);

            Esurv_eXb1[iFactor] += survivalTempo * eXb(iObs,theCause) * factorTempo;
			if(hazardType){
			  for(int iCause=0; iCause<nCause; iCause++){
				// Rcout << ".";	   
				Esurv_hazard1_eXb(iCause,iFactor) += survivalTempo * hazard1Tempo * eXb(iObs,iCause) * factorTempo;
				if(nVar[theCause]>0){
				  if(iCause==theCause){
					Esurv_hazard1_X_cumhazard[iCause].col(iFactor) += survivalTempo * hazard1Tempo * ls_tX[iCause].col(iObs) * (cumhazardTempo[iCause] - 1) * factorTempo;
				  }else{
					Esurv_hazard1_X_cumhazard[iCause].col(iFactor) += survivalTempo * hazard1Tempo * ls_tX[iCause].col(iObs) * cumhazardTempo[iCause] * factorTempo;
				  }
				}
			  }						   
			}else{
			  Esurv_hazard1_eXb(1,iFactor) += survivalTempo * hazard1Tempo * eXb(iObs,1) * factorTempo;
			  if(nVar[0]>0){
				Esurv_hazard1_X_cumhazard[0].col(iFactor) -= survivalTempo * hazard1Tempo * ls_tX[0].col(iObs) * factorTempo;
			  }
			  if(nVar[1]>0){
				Esurv_hazard1_X_cumhazard[1].col(iFactor) += survivalTempo * hazard1Tempo * ls_tX[1].col(iObs) * cumhazardTempo[1] * factorTempo;
			  }
			}
		  }
		}
      }
      // Rcout << "done) ";
		   
      // *** update integral
      // Rcout << "(integral ";
      int_hazard1 += ls_IFhazard[levelStrata(iStrata,theCause)].col(iJump) * Esurv_eXb1;
      for(int iCause=0; iCause<nCause; iCause++){
		if(iJump>0){
		  int_cumhazard += ls_IFcumhazard[iCause][levelStrata(iStrata,iCause)].col(iJump-1) * Esurv_hazard1_eXb.row(iCause);
		}
		if(nVar[iCause]>0){
		  int_bmissingIF[iCause] += Esurv_hazard1_X_cumhazard[iCause];
		}
      }	 
      // Rcout << " done) ";
	  	 
      // *** export
      // Rcout << "(update ";
      if(iiTau < nTau && tauIndex[iiTau]==iJump && tau[iiTau] <= JumpMax[iStrata]){
		int_b.fill(0.0);
		for(int iCause=0; iCause<nCause; iCause++){
		  if(nVar[iCause]>0){
			// Rcout << int_bmissingIF[iCause](0,0) << " ";
			int_b += ls_IFbeta[iCause] * int_bmissingIF[iCause];
		  }
		}

		int_total = (int_hazard1 - int_cumhazard - int_b);
		// Rcout << int_hazard1(0,0) << "/" <<  int_cumhazard(0,0) << "/" << int_b(0,0) << ":" << int_total(0,0) << endl;

		// int should be multiplied by the prevalence (n.strata/n.total) and divided by n.strata to change the sum in average for the Eterms
		// so in total we just divide by n.total
		int_total /= nNewObs;

		while(iiTau < nTau && tauIndex[iiTau]==iJump && tau[iiTau] <= JumpMax[iStrata]){
	    
		  // Rcout << "*" ;
		  for(int iFactor=0; iFactor<nFactor; iFactor++){
			out[iFactor].col(iiTau) += int_total.col(iFactor);
		  }
		  iiTau ++;
		}
      }
	   
      // Rcout << "done) " << endl;		   
    }

    // ** fill remaining columns with NA
    if(iiTau < nTau){
      for(int iFactor=0; iFactor<nFactor; iFactor++){
		out[iFactor].cols(iiTau, nTau-1).fill(NA_REAL);
      }
    } 

  }


  return(out);
}


// * additional functions


// ** calcIFhazard
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
