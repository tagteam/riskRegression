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
				   const arma::mat& survival,
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

      iTau = iTau0; 
      std::fill(veci_IF_risk.begin(), veci_IF_risk.end(), 0);
     	
      for(int iJump=0; iJump<nJump; iJump++){
		if(jumpTime[iJump]>lastSampleTime){break;}

		//// compute the hazard for the cause of interest ////
		i_hazard = hazard0[theCause][iJump]*newEXb(iNewObs,theCause);

		//// compute the survival ////
		// it is survival at t- which is stored i.e. the survival at the previous eventtime (censoring does not affect survival)
		i_survival = survival(iNewObs,iJump);

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
List calcSeCif2_cpp(const std::vector<arma::mat>& ls_IFbeta, const std::vector<arma::mat>& ls_X,
					const std::vector<arma::mat>& ls_cumhazard, const arma::mat& ls_hazard, const arma::mat& survival,
					const std::vector< std::vector<arma::mat> >& ls_IFcumhazard, const std::vector<arma::mat>& ls_IFhazard,
					const NumericMatrix& eXb,
					int nJumpTime, const NumericVector& JumpMax,
					const NumericVector& tau, const arma::vec& tauIndex, int nTau,
					int nObs,  
					int theCause, int nCause, bool hazardType, arma::vec nVar,
					int nNewObs, NumericMatrix strata,
					bool exportSE, bool exportIF, bool exportIFsum, bool diag){

  arma::mat X_IFbeta;

  arma::mat IFhazard;
  arma::mat IFcumhazard;
  arma::colvec hazard;

  double ieXb;
  unsigned int iStrataCause;

  arma::colvec IF_tempo;
  arma::colvec cumIF_tempo;

  // ** initialize for export
  arma::mat outSE;
  if(exportSE){
	if(diag){
      outSE.resize(nNewObs,1);
    }else{
      outSE.resize(nNewObs,nTau);
    }
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

  int startObs;
  int iiTau;
  int iNTau = 0; // false initialization to avoid warning
  int iNJumpTime = 0; // false initialization to avoid warning
  if(diag){
    startObs = iTau;
  }else{
    startObs = 0;
    iNTau = nTau;
    iNJumpTime = nJumpTime;
  }
  
  // ** prepare the influence function
  arma::uvec iUvec_linspace(1);
  arma::uvec iUvec_strata(1);

  for(int iNewObs=startObs; iNewObs<nNewObs; iNewObs++){

    R_CheckUserInterrupt();

    if(iTau>=nTau){continue;}
    
    if(diag){
      iiTau = iNewObs;
      iNTau = iNewObs + 1;
      iNJumpTime = tauIndex[iNewObs]+1;
    }else{
      iiTau = iTau;
    }
    // Rcout << "start: " << iiTau << " " << iNTau << endl;
	
    if(diag){
      iUvec_linspace = linspace<uvec>(0, iNJumpTime-1, iNJumpTime);
      IFcumhazard = zeros<mat>(nObs,iNJumpTime);
    }else{
      IFcumhazard = zeros<mat>(nObs,nJumpTime);
    }
    ieXb = NA_REAL;
   
    for(int iCause=0; iCause<nCause; iCause ++){

      iStrataCause = strata(iNewObs,iCause);
      if(diag){
		iUvec_strata(0) = iStrataCause;
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
	    IFcumhazard += ieXb * (ls_IFcumhazard[iCause][iStrataCause].cols(iUvec_linspace) + X_IFbeta * ls_tcumhazard[iCause].submat(iUvec_strata,iUvec_linspace));
	  }else{
	    IFcumhazard += ieXb * (ls_IFcumhazard[iCause][iStrataCause] + X_IFbeta * ls_tcumhazard[iCause].row(iStrataCause));
	  }
	}else{
	  if(diag){
	    IFcumhazard += ls_IFcumhazard[iCause][iStrataCause].cols(iUvec_linspace);
	  }else{
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

      // Rcout << "5 " ;
      // prepare IF
      if(hazard[iJump]>0){
	if(iJump==0){
	  IF_tempo = IFhazard.col(iJump);
	}else{
	  // survival is evaluated just before the jump
	  IF_tempo = (IFhazard.col(iJump) - IFcumhazard.col(iJump-1) * hazard[iJump]) * survival(iNewObs,iJump);
	}
	cumIF_tempo = cumIF_tempo + IF_tempo;
      }

      // Rcout << "6";
      // store
      // Rcout << "test: " << tauIndex[iiTau] << " " << iJump << "/ " << nJumpTime << " | " << tau[iiTau] << " " << JumpMax[iNewObs] << endl;
      while((iiTau < iNTau) && (tauIndex[iiTau] == iJump) && (tau[iiTau] <= JumpMax[iNewObs])){

	if(exportSE){
	  // Rcout << "a";
	  if(diag){
	    outSE.row(iNewObs).col(0) = sqrt(accu(pow(cumIF_tempo,2)));
	  }else{
	    outSE.row(iNewObs).col(iiTau) = sqrt(accu(pow(cumIF_tempo,2)));
	  }
	}
	if(exportIF){
	  // Rcout << "b";
	  if(diag){
	    outIF.slice(0).row(iNewObs) = cumIF_tempo.t();
	  }else{
	    outIF.slice(iiTau).row(iNewObs) = cumIF_tempo.t();
	  }
	}
	if(exportIFsum){
	  // Rcout << "c";
	  if(diag){
	    outIFsum.col(0) += cumIF_tempo;
	  }else{
	    outIFsum.col(iiTau) += cumIF_tempo;
	  }
	}
	// Rcout << "increment: " << iiTau << " " << iNTau << endl;
	iiTau++;
      }
      // Rcout << "-end " << endl;
      if(iiTau == iNTau){break;} 

    }

    // Rcout << "end: " << iiTau << " " << iNTau << endl;

    // ** fill remaining columns with NA
    while(iiTau < iNTau){
      if(exportSE){
	if(diag){
	  outSE.row(iNewObs).col(0).fill(NA_REAL);
	}else{
	  outSE.row(iNewObs).col(iiTau).fill(NA_REAL);
	}
      }
      if(exportIF){
	// Rcout << "b";
	if(diag){
	  outIF.slice(0).row(iNewObs).fill(NA_REAL);
	}else{
	  outIF.slice(iiTau).row(iNewObs).fill(NA_REAL);
	}
      }
      if(exportIFsum){
	if(diag){
	  outIFsum.col(0).fill(NA_REAL);
	}else{
	  outIFsum.col(iiTau).fill(NA_REAL);
	}
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
// (not finished, no valid, slower than R version)
// arma::mat calcAIFcif_cpp(const arma::rowvec& hazard0_cause,
// 						 const std::vector<arma::rowvec>& cumhazard0,
// 						 const arma::mat& IFhazard0_cause,
// 						 const std::vector<arma::mat>& IFcumhazard0,
// 						 const std::vector<arma::mat>& IFbeta,
// 						 const arma::mat& eXb1_S,
// 						 const std::vector<arma::mat>& eXb1_S_eXbj,
// 						 const std::vector<arma::mat>& eXb1_S_X1,
// 						 const std::vector<std::vector<arma::mat>>& eXb1_S_Xj_eXbj,
// 						 const arma::vec weight, const arma::mat factor,
// 						 int nJump, const arma::uvec subsetJump,
// 						 int nCause, std::vector<bool> test_allCause, std::vector<bool> test_theCause,
// 						 std::vector<int> nVar){
	

//   arma::mat AIF;
//   arma::mat AIFtempo;
//   arma::mat Ep;
//   arma::mat Etempo;
  
//   // *** term 1
//   // Rcout << "term1: " ;
//   AIF = IFhazard0_cause.cols(subsetJump);
//   AIF.each_row() %= sum(eXb1_S.cols(subsetJump) % factor,0) / weight;
//   // Rcout << "end " << endl;
  
//   for(int iCause=0; iCause<nCause; iCause++){

// 	// *** term 3
// 	if(test_allCause[iCause]){
// 	  // Rcout << "term3: " ;
//   	  AIFtempo = IFcumhazard0[iCause].cols(subsetJump);
// 	  AIFtempo.each_row() %= hazard0_cause.cols(subsetJump) % sum(eXb1_S_eXbj[iCause].cols(subsetJump) % factor,0) / weight;
// 	  AIF += AIF - AIFtempo;
// 	  // Rcout << "end " << endl;
//   	}
             
// 	if(nVar[iCause]>0){ 
// 	  Ep.resize(nVar[iCause],nJump);
// 	  Ep.fill(0.0);

// 	  // *** term 2
// 	  if(test_theCause[iCause]){
// 		// Rcout << "term2: " ;
// 		for(int iP=0; iP<nVar[iCause]; iP++){
// 		  Ep.row(iP) += sum(eXb1_S_X1[iP].cols(subsetJump) % factor, 0) / weight;
// 		}
// 		// Rcout << "end " << endl ;
// 	  }
                 
// 	  // *** term 4
// 	  if(test_allCause[iCause]){
// 		// Rcout << "term4: " ;
  
// 		for(int iP=0; iP<nVar[iCause]; iP++){
// 		  Etempo = sum(eXb1_S_Xj_eXbj[iCause][iP].cols(subsetJump) % factor, 0) / weight;
// 		  Etempo.each_row() %= cumhazard0[iCause].cols(subsetJump);
// 		  Ep.row(iP) -= Etempo;
// 		}
// 		// Rcout << "end " << endl ;
// 	  }

// 	  // Rcout << "store: " ;
// 	  Ep.each_row() %= hazard0_cause.cols(subsetJump);
// 	  AIF += IFbeta[iCause] * Ep;
// 	  // Rcout << "end " << endl ;
// 	}
//   }

//   AIF = sum(AIF,1);
//   return(AIF);
// }


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

