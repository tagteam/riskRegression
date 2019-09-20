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

arma::mat calcSurvBeforeJump_cpp(const std::vector<arma::mat>& ls_cumhazard,
						   const arma::mat& eXb,
						   unsigned int nCause, unsigned int theCause, bool hazardType,
						   int nNewObs, int nJumpTime,
						   const arma::umat& Ustrata, unsigned int nStrata, const std::vector<arma::uvec>& ls_indexStrata);

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
    SEcumhazard.fill(0.0);
  }
  arma::cube IFcumhazard;
  if(exportIF){
    IFcumhazard.resize(nNewObs, nTau, nSample);
    IFcumhazard.fill(0.0);
  }
  arma::mat IFsum_cumhazard;
  if(exportIFsum_cumhazard){
    IFsum_cumhazard.resize(nSample, nTau);
    IFsum_cumhazard.fill(0.0);
  }
  arma::mat IFsum_survival;
  if(exportIFsum_survival){
    IFsum_survival.resize(nSample, nTau);
    IFsum_survival.fill(0.0);
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
													  int diag,
													  bool exportCumHazard,
													  bool exportSurvival){

  // ** prepare output
  int nFactor = factor.size();
  std::vector< std::vector<arma::mat> > out(2);
  if(exportCumHazard){
    out[0].resize(nFactor);
    for(int iFactor=0; iFactor<nFactor; iFactor++){
	  if(diag){
		out[0][iFactor].resize(nObs, 1);
	  }else{
		out[0][iFactor].resize(nObs, nTimes);
	  }
      out[0][iFactor].fill(0.0);
    }
  }
  if(exportSurvival){
    out[1].resize(nFactor);
    for(int iFactor=0; iFactor<nFactor; iFactor++){
	  if(diag){
		out[1][iFactor].resize(nObs, 1);
	  }else{
		out[1][iFactor].resize(nObs, nTimes);
	  }
      out[1][iFactor].fill(0.0);
    }
  }

  // ** initialization
  int nObs_strata;

  arma::mat iW_eXb;
  arma::mat iW_eXb_S;
  
  arma::mat iAIF_H;
  arma::mat iAIF_S;
  
  arma::mat iE_W_eXb_cumhazard0_X;
  arma::mat iE_W_eXb_cumhazard0_X_S;
  
  for(int iStrata=0; iStrata<nStrata; iStrata++){

	nObs_strata = ls_indexStrata[iStrata].size();
	
	for(int iFactor=0; iFactor<nFactor; iFactor++){

	  iW_eXb = factor[iFactor].rows(ls_indexStrata[iStrata]);
	  if(nVar>0){
		iW_eXb.each_col() %= eXb.rows(ls_indexStrata[iStrata]);
	  }
	  if(exportSurvival){
		iW_eXb_S = survival.rows(ls_indexStrata[iStrata]) % iW_eXb;
	  }
	  
	  if(diag==1){
		// <IF>(cumhazard) = E[w * eXb * IF_cumhazard0] + E[w * eXb * cumhazard0 * X] * IF_beta
		// <IF>(survival) = -(E[w * Surv * eXb * IF_cumhazard0] + E[w * Surv * eXb * cumhazard0 * X] * IF_beta)

		// first term
		if(exportCumHazard){
		  // Rcout << "b: ";
		  iAIF_H = ls_IFcumhazard[iStrata].cols(ls_indexStrata[iStrata]);
		  iAIF_H.each_row() %= trans(iW_eXb);
		  iAIF_H = sum(iAIF_H,1)/nObs_strata;
		  // Rcout << endl;
		}
		if(exportSurvival){
		  iAIF_S = ls_IFcumhazard[iStrata].cols(ls_indexStrata[iStrata]);
		  iAIF_S.each_row() %= trans(iW_eXb_S);
		  iAIF_S = sum(iAIF_S,1)/nObs_strata;
		}

		// second term
		if(nVar>0){
		  if(exportCumHazard){
			// Rcout << "c: ";
			iE_W_eXb_cumhazard0_X = X.rows(ls_indexStrata[iStrata]);
			iE_W_eXb_cumhazard0_X.each_col() %= trans(cumhazard0[iStrata].cols(ls_indexStrata[iStrata])) % iW_eXb;
			iE_W_eXb_cumhazard0_X = sum(iE_W_eXb_cumhazard0_X,0)/nObs_strata;
			iAIF_H += IFbeta * trans(iE_W_eXb_cumhazard0_X);
			// Rcout << endl;
		  }
		  if(exportSurvival){
			iE_W_eXb_cumhazard0_X_S = X.rows(ls_indexStrata[iStrata]);
			iE_W_eXb_cumhazard0_X_S.each_col() %= trans(cumhazard0[iStrata].cols(ls_indexStrata[iStrata])) % iW_eXb_S;
			iE_W_eXb_cumhazard0_X_S = sum(iE_W_eXb_cumhazard0_X_S,0)/nObs_strata;
			iAIF_S += IFbeta * trans(iE_W_eXb_cumhazard0_X_S);
		  }		  
		}
		  
	  }else{
		// <IF>(cumhazard) = E[w * eXb] IF_cumhazard0 + E[w * eXb * X] * cumhazard0 * IF_beta
		// <IF>(survival) = -(E[w * Surv * eXb] IF_cumhazard0 + E[w * Surv * eXb * X] * cumhazard0 * IF_beta)

		// first term
		if(exportCumHazard){
		  iAIF_H = ls_IFcumhazard[iStrata];
		  iAIF_H.each_row() %= sum(iW_eXb,0)/nObs_strata;
		}
		if(exportSurvival){
		  iAIF_S = ls_IFcumhazard[iStrata];		  
		  iAIF_S.each_row() %= sum(iW_eXb_S,0)/nObs_strata;
		}
		  
		// second term
		if(nVar>0){
		  if(exportCumHazard){
			iE_W_eXb_cumhazard0_X = (trans(X.rows(ls_indexStrata[iStrata])) * iW_eXb) / nObs_strata;			
			iE_W_eXb_cumhazard0_X.each_row() %= cumhazard0[iStrata];
			iAIF_H += IFbeta * iE_W_eXb_cumhazard0_X;
		  }
		  if(exportSurvival){
			iE_W_eXb_cumhazard0_X_S = (trans(X.rows(ls_indexStrata[iStrata])) * iW_eXb_S) / nObs_strata;			
			iE_W_eXb_cumhazard0_X_S.each_row() %= cumhazard0[iStrata];
			iAIF_S += IFbeta * iE_W_eXb_cumhazard0_X_S;
		  }
		}
	  }
			
	  // update - multiplying by the prevalence of the strata
	  if(exportCumHazard){
		out[0][iFactor] += iAIF_H * prevStrata[iStrata];
	  }
	  if(exportSurvival){
		out[1][iFactor] -= iAIF_S * prevStrata[iStrata];
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
    SErisk.fill(0.0);
  }  
  arma::cube IFrisk;
  if(exportIF){
    IFrisk.resize(nNewObs, nTau, nSample);
    IFrisk.fill(0.0);
  }
  arma::mat IFsumrisk;
  if(exportIFsum){
    IFsumrisk.resize(nSample, nTau);
    IFsumrisk.fill(0.0);
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
					bool exportSE, bool exportIF, bool exportIFsum, bool diag){

  arma::mat X_IFbeta;

  arma::mat IFhazard;
  arma::mat IFcumhazard;
  arma::colvec hazard;
  arma::colvec cumhazard;

  double ieXb;
  unsigned int iStrataCause;

  arma::colvec IF_tempo;
  arma::colvec cumIF_tempo;

  // ** initialize for export
  arma::mat outSE;
  if(exportSE){
	outSE.resize(nNewObs,nTau);
    outSE.fill(0.0);
  }
  arma::cube outIF;
  if(exportIF){
	if(diag){
	  outIF.resize(nNewObs,nObs,1);
	}else{
	  outIF.resize(nNewObs,nObs,nTau);
	}
    outIF.fill(0.0);
  }
  arma::mat outIFsum;
  if(exportIFsum){
	if(diag){
	  outIFsum.resize(nObs,1);
	}else{
	  outIFsum.resize(nObs,nTau);
	}
    outIFsum.fill(0.0);
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
  int iNTau;
  int iNJumpTime;
  NumericVector iVTau(1);
  arma::vec iTauIndex(1);
  arma::uvec iUvec_linspace;
  arma::uvec iUvec_strata;
  
  // ** prepare the influence function
  for(int iNewObs=0; iNewObs<nNewObs; iNewObs++){

    R_CheckUserInterrupt();
     
    if(diag){
	  iiTau = iTau;
	  iNTau = 1;
	  iVTau[0] = tau[iNewObs];
	  iTauIndex[0] = tauIndex[iNewObs];
	  iNJumpTime = iTauIndex[0]+1;
	}else{
	  iiTau = iTau;
	  iNTau = nTau;
	  iVTau = tau;
	  iTauIndex = tauIndex;
	  iNJumpTime = nJumpTime;
	}
	
	if(iTau>=iNTau){continue;}
	
	if(diag){
	  iUvec_linspace = linspace<uvec>(0, iNJumpTime-1, iNJumpTime);
	  IFcumhazard = zeros<mat>(nObs,iNJumpTime);
	  cumhazard = zeros<colvec>(iNJumpTime);
	}else{
	  IFcumhazard = zeros<mat>(nObs,nJumpTime);
	  cumhazard = zeros<colvec>(nJumpTime);
	}
    ieXb = NA_REAL;
   
    for(int iCause=0; iCause<nCause; iCause ++){

      iStrataCause = strata(iNewObs,iCause);
	  if(diag){
		iUvec_strata = {iStrataCause};
	  }
	  
      // Rcout << "2 ";
      if(nVar[iCause]>0){
		X_IFbeta = ls_IFbeta[iCause] * (ls_X[iCause].row(iNewObs)).t();
		ieXb = eXb(iNewObs,iCause);
      }

      // Rcout << "3 ";
      if(hazardType || (iCause != theCause)){
		if(nVar[iCause]>0){
		  if(diag){
			cumhazard += ieXb * ls_cumhazard[iCause].submat(iUvec_linspace,iUvec_strata);
			IFcumhazard += ieXb * (ls_IFcumhazard[iCause][iStrataCause].cols(iUvec_linspace) + X_IFbeta * ls_tcumhazard[iCause].submat(iUvec_strata,iUvec_linspace));
		  }else{
			cumhazard += ieXb * ls_cumhazard[iCause].col(iStrataCause);
			IFcumhazard += ieXb * (ls_IFcumhazard[iCause][iStrataCause] + X_IFbeta * ls_tcumhazard[iCause].row(iStrataCause));
		  }
		}else{
		  if(diag){
			cumhazard += ls_cumhazard[iCause].submat(iUvec_linspace,iUvec_strata);
			IFcumhazard += ls_IFcumhazard[iCause][iStrataCause].cols(iUvec_linspace);
		  }else{
			cumhazard += ls_cumhazard[iCause].col(iStrataCause);
			IFcumhazard += ls_IFcumhazard[iCause][iStrataCause];
		  }
		}
      }

      // Rcout << "4 ";
      if(iCause == theCause){
		if(nVar[iCause] > 0){
		  if(diag){
			hazard = ieXb * ls_hazard.submat(iUvec_linspace,iUvec_strata);
			IFhazard = ieXb * (ls_IFhazard[iStrataCause].cols(iUvec_linspace) + X_IFbeta * ls_thazard.submat(iUvec_strata,iUvec_linspace));
		  }else{
			hazard = ieXb * ls_hazard.col(iStrataCause);
			IFhazard = ieXb * (ls_IFhazard[iStrataCause] + X_IFbeta * ls_thazard.row(iStrataCause));
		  }
		}else{
		  if(diag){
			hazard = ls_hazard.submat(iUvec_linspace,iUvec_strata);
			IFhazard = ls_IFhazard[iStrataCause].cols(iUvec_linspace);
		  }else{
			hazard = ls_hazard.col(iStrataCause);
			IFhazard = ls_IFhazard[iStrataCause];		
		  }
		}
	  }
	}
     
    // ** loop over time
    cumIF_tempo = zeros<colvec>(nObs);   

    for(int iJump=0; iJump<iNJumpTime; iJump++){

      // Rcout << "5" ;
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
      while(iiTau < iNTau && iTauIndex[iiTau] == iJump && iVTau[iiTau] <= JumpMax[iNewObs]){

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
		  outIFsum.col(iiTau) += cumIF_tempo;
		}
		iiTau++;
      }
      // Rcout << "-end " << endl;
      if(iiTau == iNTau){break;} 

    }

    // ** fill remaining columns with NA
    while(iiTau < iNTau){
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

  if(exportIFsum){
	outIFsum /= nNewObs;
  }

   
  // ** export
  return(List::create(Named("se") = outSE,
					  Named("iid") = outIF,
					  Named("average.iid") = outIFsum));

}

// ** calcAIFcif_cpp: compute average IF for the absolute risk (method 3)
// [[Rcpp::export]]
std::vector<arma::mat> calcAIFcif_cpp(const arma::rowvec& hazard1,
									  const arma::mat& IFhazard1,
									  const std::vector<arma::rowvec>& cumhazard,
									  const std::vector<arma::mat>& IFcumhazard,
									  const std::vector<arma::mat>& IFbeta,
									  const arma::mat& eXb1_S,
									  const std::vector<arma::mat>& eXb1_S_eXbj,
									  const std::vector<arma::mat>& eXb1_S_X1,
									  const std::vector<std::vector<arma::mat>>& eXb1_S_Xj_eXbj,
									  int nObs, int nJump, int nTimes,
									  arma::uvec validTimes, arma::uvec NATimes, arma::uvec sindexTimes, 
									  int nCause, std::vector<bool> test_allCause, std::vector<bool> test_theCause,
									  std::vector<int> nVar,
									  const std::vector<arma::mat>& factor,
									  bool diag, bool is_factor2){
	
  // ** prepare
  int nFactor = factor.size();
  int nFactor2;
  std::vector<arma::mat> out(nFactor);

  arma::mat Mtempo;
  
  arma::rowvec iW;
  arma::uvec iJump;

  arma::mat iTerm1;
  arma::mat E_W_eXb1_S;
    
  arma::mat iTerm3;
  arma::mat E_W_eXb1_S_eXbj;
  
  arma::mat iTerm24;
  arma::mat E_W_eXb1_S_X1;
  arma::mat E_W_eXb1_S_Xj_eXbj;
  
  arma::mat iAIF;
  
  // ** loop over factors (level 1)
  for(int iFactor=0; iFactor<nFactor; iFactor++){
	
	// *** initialization
	if(diag){
	  out[iFactor].resize(nObs, 1);

	  // fill
	  out[iFactor].fill(0.0);
	  if(NATimes.size()>0){
		Mtempo.resize(NATimes.size(),1);
		Mtempo.fill(NA_REAL);
		out[iFactor].rows(NATimes) = Mtempo;
	  }
	}else{
      out[iFactor].resize(nObs, nTimes);

	  // fill
	  out[iFactor].fill(0.0);
	  if(NATimes.size()>0){
		Mtempo.resize(nObs,NATimes.size());
		Mtempo.fill(NA_REAL);
		out[iFactor].cols(NATimes) = Mtempo;
	  }
	}

	nFactor2 = factor[iFactor].n_cols;
	iAIF.resize(nObs,nJump);
	iAIF.fill(0.0);
	
	// ** loop over factors (level 2)
	for(int iFactor2=0; iFactor2<nFactor2; iFactor2++){
	  iW = factor[iFactor].col(iFactor2);

	  // *** term 1
	  E_W_eXb1_S = eXb1_S;
	  E_W_eXb1_S.each_col() %= iW;
	  E_W_eXb1_S = sum(E_W_eXb1_S,1)/nObs;

	  iTerm1 = IFhazard1;
	  iTerm1 %= E_W_eXb1_S;
	  iAIF += iTerm1;
	  
	  for(int iCause=0; iCause<nCause; iCause++){

		// *** term 3
		if(test_allCause[iCause]){
		  E_W_eXb1_S_eXbj = eXb1_S_eXbj[iCause];
		  E_W_eXb1_S_eXbj.each_col() %= iW;
		  E_W_eXb1_S_eXbj = sum(E_W_eXb1_S_eXbj,1) % hazard1 /nObs;
			
		  iTerm3 = IFcumhazard[iCause];
		  iTerm3 %= E_W_eXb1_S_eXbj;
		  iAIF += iTerm3;
		}

		if(nVar[iCause]>0){
		  iTerm24.resize(nVar[iCause],nJump);
		  iTerm24.fill(0.0);

		  // *** term 2
		  if(test_theCause[iCause]){
			for(int iP=0; iP<nVar[iCause]; iP++){
			  E_W_eXb1_S_X1 = eXb1_S_X1[iP];
			  E_W_eXb1_S_X1.each_col() %= iW;
			  E_W_eXb1_S_X1 = sum(E_W_eXb1_S_X1,1)/nObs;

			  iTerm24.row(iP) += E_W_eXb1_S_X1 * hazard1;
			}
		  }

		  // *** term 4
		  if(test_allCause[iCause]){
			for(int iP=0; iP<nVar[iCause]; iP++){
			  E_W_eXb1_S_Xj_eXbj = eXb1_S_Xj_eXbj[iCause][iP];
			  E_W_eXb1_S_Xj_eXbj.each_col() %= iW;
			  E_W_eXb1_S_Xj_eXbj = sum(E_W_eXb1_S_Xj_eXbj,1)/nObs;

			  iTerm24.row(iP) += E_W_eXb1_S_Xj_eXbj * cumhazard[iCause] * hazard1;
			}
		  }

		  iAIF += IFbeta[iCause] * iTerm24;
		}
	  }

	  iAIF = cumsum(iAIF,1);
	  
	  if(diag){
		
	  }else{
		out[iFactor].cols(validTimes) += iAIF.cols(sindexTimes);
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

