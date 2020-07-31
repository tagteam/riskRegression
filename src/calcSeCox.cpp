// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// * calcSeMinimalCox_cpp: compute IF/sumIF/se for the hazard / cumlative hazard / survival (method 1)
// [[Rcpp::export]]
// J: number of jump times
// n: number of observations in the training set
// N: number of observations in the prediction set
// p: number of regressors
// S: number of strata
// T: number of prediction times
List calcSeMinimalCox_cpp(const arma::vec& seqTau, // horizon time for the predictions (T)
						  const arma::mat& newSurvival, // predicted survival for all observations at each horizon time (NxT)
						  const std::vector< arma::vec >& hazard0, // baseline hazard for each strata  S:T
						  const std::vector< arma::vec > & cumhazard0, // baseline cumulative hazard for each strata S:T
						  const arma::mat& newX, // design matrix (Nxp)
						  const arma::vec& neweXb, // exponential of the linear predictor (N)
						  const arma::mat& IFbeta, //influence function for the regression coefficients (nxp)
						  const std::vector< arma::mat >& Ehazard0, // to compute the influence function of the baseline hazard S:(pxJ)
						  const std::vector< arma::mat >& cumEhazard0, // to compute the influence function of the baseline hazard S(pxJ)
						  const std::vector< arma::vec >& hazard_iS0, // to compute the influence function of the baseline hazard S:J
						  const std::vector< arma::vec >& cumhazard_iS0, // to compute the influence function of the baseline hazard S:J
						  const arma::mat& delta_iS0, // to compute the influence function of the baseline hazard (nxS, S because multiplied by the strata indicator)
						  const arma::mat& sample_eXb, // to compute the influence function of the baseline hazard (nxS, S because multiplied by the strata indicator)
						  const arma::vec& sample_time, // event times of the training set (n)
						  const std::vector< arma::uvec>& indexJumpSample_time, // index of the jump time corresponding to the sample times S:J
						  const std::vector< arma::vec>& jump_time, // jump times for each strata S:J
						  const std::vector< arma::uvec >& indexJumpTau, // index of the jump time corresponding to the horizon times S:J
						  const arma::vec& lastSampleTime, // time after which the prediction are NA (S)
						  const std::vector< arma::uvec>& newdata_index, // index of the observations within strata
						  const std::vector<arma::mat>& factor,
						  int nTau, int nNewObs, int nSample, int nStrata, int p, 
						  bool diag, bool exportSE, bool exportIF, bool exportIFmean,
						  bool exportHazard, bool exportCumhazard, bool exportSurvival,
						  int debug){

  // ** prepare
  if(debug>0){Rcpp::Rcout << "Prepare" << std::endl;}
  int iStrata_tauMin,iStrata_tauMax; 
  arma::vec iStrata_seqTau;
  arma::uvec iStrata_indexJumpTau;
  int iStrata_nNewObs; 
  int iJump; 
  int iTau,iTauStore;
  
  int iNewObs2; 
  arma::colvec iStrata_IFhazard0, iStrata_IFcumhazard0;
  arma::colvec iStrata_IFhazard, iStrata_IFcumhazard, iStrata_IFsurvival;
  
  arma::uvec index_timestop(nSample);
  arma::uvec tempo_uvec(1);
	
  int nFactor = factor.size();
  arma::vec iStrata_weXb,iStrata_weXbS;
  arma::mat iStrata_weXbX,iStrata_weXbXS;
  arma::vec iSurvival;

  // ** initialize
  if(debug>0){Rcpp::Rcout << "Initialize" << std::endl;}
  arma::cube IF_hazard;
  std::vector< arma::mat > IFmean_hazard(nFactor);
  if(exportHazard){
	if(exportIF){
	  IF_hazard.resize(nSample, nNewObs, nTau);
	  IF_hazard.fill(0.0);
	}
	if(exportIFmean){
	  for(int iFactor=0; iFactor<nFactor; iFactor++){
		IFmean_hazard[iFactor].resize(nSample, nTau);
		IFmean_hazard[iFactor].fill(0.0);
	  }
	}
  }
  
  arma::cube IF_cumhazard;
  arma::mat SE_cumhazard;
  std::vector< arma::mat > IFmean_cumhazard(nFactor);
  if(exportCumhazard){
	if(exportIF){
	  IF_cumhazard.resize(nSample, nNewObs, nTau);
	  IF_cumhazard.fill(0.0);
	}
	if(exportSE){
	  SE_cumhazard.resize(nNewObs, nTau);
	  SE_cumhazard.fill(0.0);
	}
	if(exportIFmean){
	  for(int iFactor=0; iFactor<nFactor; iFactor++){
		IFmean_cumhazard[iFactor].resize(nSample, nTau);
		IFmean_cumhazard[iFactor].fill(0.0);
	  }
	}
  }

  arma::cube IF_survival;
  arma::mat SE_survival;
  std::vector< arma::mat > IFmean_survival(nFactor);
  if(exportSurvival){
	if(exportIF){
	  IF_survival.resize(nSample, nNewObs, nTau);
	  IF_survival.fill(0.0);
	}
	if(exportSE){
	  SE_survival.resize(nNewObs, nTau);
	  SE_survival.fill(0.0);
	}
	if(exportIFmean){
	  for(int iFactor=0; iFactor<nFactor; iFactor++){
		IFmean_survival[iFactor].resize(nSample, nTau);
		IFmean_survival[iFactor].fill(0.0);
	  }
	}
  }
  
  
  // ** compute IF within each strata
  if(debug>0){Rcpp::Rcout << "Compute IF" << std::endl;}
  for(int iStrata=0; iStrata<nStrata ; iStrata++){
	if(debug>0){Rcpp::Rcout << " > strata " << iStrata << "/" << (nStrata-1) << " ";}

	// *** check if any observation in strata
	iStrata_nNewObs = newdata_index[iStrata].size();
	if(iStrata_nNewObs==0){continue;}

	// *** jump times
	if(jump_time[iStrata].size()==0){continue;}

	// *** narrow down prediction times
    iStrata_tauMin=0;
    if(diag){
	  iStrata_tauMax = iStrata_nNewObs-1;
	  iStrata_seqTau = seqTau(newdata_index[iStrata]);
	  iStrata_indexJumpTau = indexJumpTau[iStrata](newdata_index[iStrata]);
	  iStrata_nNewObs = 1;
	}else{
	  iStrata_tauMax = nTau-1;
	  iStrata_seqTau = seqTau;
	  iStrata_indexJumpTau = indexJumpTau[iStrata];
	}

	// WARNING: this part needs to be before iStrata_tauMax is modified (i.e. the next while loop)
	while((iStrata_tauMin <= iStrata_tauMax) && jump_time[iStrata](0)>iStrata_seqTau(iStrata_tauMin)){ // start at the first event or after 
	  iStrata_tauMin++;    
	}
	if(iStrata_tauMin > iStrata_tauMax){continue;}

	while((iStrata_tauMax >= 0) && seqTau(iStrata_tauMax)>lastSampleTime(iStrata)){ // end at the last event or before
	  for(int iNewObs=0; iNewObs<iStrata_nNewObs ; iNewObs++){		
		if(diag){
		  iTauStore = 0;
		  iNewObs2 = newdata_index[iStrata](iStrata_tauMax);
		}else{
		  iTauStore = iStrata_tauMax;
		  iNewObs2 = newdata_index[iStrata](iNewObs);
		}
		
		if(exportHazard){
		  if(exportIF){IF_hazard.slice(iTauStore).col(iNewObs2).fill(NA_REAL);}
		  if(exportIFmean){
			for(int iFactor=0; iFactor<nFactor; iFactor++){
			  IFmean_hazard[iFactor].col(iTauStore).fill(NA_REAL);
			}
		  }
		}

		if(exportCumhazard){
		  if(exportIF){IF_cumhazard.slice(iTauStore).col(iNewObs2).fill(NA_REAL);}
		  if(exportSE){SE_cumhazard(iNewObs2,iTauStore) = NA_REAL;}
		  if(exportIFmean){
			for(int iFactor=0; iFactor<nFactor; iFactor++){
			  IFmean_cumhazard[iFactor].col(iTauStore).fill(NA_REAL);
			}
		  }
		}

		if(exportSurvival){
		  if(exportIF){IF_survival.slice(iTauStore).col(iNewObs2).fill(NA_REAL);}
		  if(exportSE){SE_survival(iNewObs2,iTauStore) = NA_REAL;}
		  if(exportIFmean){
			for(int iFactor=0; iFactor<nFactor; iFactor++){
			  IFmean_survival[iFactor].col(iTauStore).fill(NA_REAL);
			}
		  }
		}
	  }
	  iStrata_tauMax--;
	}
	if(iStrata_tauMax < 0){continue;}
	
    R_CheckUserInterrupt();

	if(debug>1){Rcpp::Rcout << " (tau=" << iStrata_tauMin << "-" << iStrata_tauMax << ") ";}
	
	// *** compute IF/SE/IFmean at each time point
	for(int iTime=iStrata_tauMin; iTime<=iStrata_tauMax; iTime++){
	  if(diag){
		iTau = newdata_index[iStrata](iTime);
		iTauStore = 0;
	  }else{
		iTau = iTime;
		iTauStore = iTime;
	  }

	  // jump
	  iJump = iStrata_indexJumpTau(iTime);

	  if(debug>1){Rcpp::Rcout << " IF0 " ;}

	  // **** IF baseline hazard 
	  if(exportHazard && (iStrata_seqTau(iTime) == jump_time[iStrata](iJump))){
	  	iStrata_IFhazard0 = delta_iS0.col(iStrata) % (sample_time == iStrata_seqTau(iTime)) - sample_eXb.col(iStrata) % (iStrata_seqTau(iTime) <= sample_time) * hazard_iS0[iStrata](iJump);
	  	if(p>0){
		  iStrata_IFhazard0 -= IFbeta * Ehazard0[iStrata].col(iTau);
	  	}
	  }else{
		iStrata_IFhazard0.resize(nSample);
		iStrata_IFhazard0.fill(0.0);
	  }

	  // **** IF baseline cumulative hazard hazard 
	  if(exportCumhazard || exportSurvival){
		index_timestop = indexJumpSample_time[iStrata];
		index_timestop.elem(find(index_timestop > iJump)).fill(iJump);

		iStrata_IFcumhazard0 = delta_iS0.col(iStrata) % (sample_time <= iStrata_seqTau(iTime)) - sample_eXb.col(iStrata) % cumhazard_iS0[iStrata](index_timestop);
	  	if(p>0){
		  iStrata_IFcumhazard0 -= IFbeta * cumEhazard0[iStrata].col(iTau);
	  	}
	  }
	  
	  // **** IF/SE hazard/cumhazard/survival
	  if(exportIF || exportSE || (exportIFmean && diag)){
		if(debug>1){Rcpp::Rcout << " IF " ;}
		for(int iNewObs=0; iNewObs<iStrata_nNewObs; iNewObs++){
		  if(diag){
			iNewObs2 = newdata_index[iStrata](iTime);
		  }else{
			iNewObs2 = newdata_index[iStrata](iNewObs);
		  }

		  if(p>0){
			if(exportHazard){
			  iStrata_IFhazard = neweXb(iNewObs2)*(iStrata_IFhazard0 + hazard0[iStrata](iTau) * IFbeta * trans(newX.row(iNewObs2)));
			}
			if(exportCumhazard || exportSurvival){
			  iStrata_IFcumhazard = neweXb(iNewObs2)*(iStrata_IFcumhazard0 + cumhazard0[iStrata](iTau) * IFbeta * trans(newX.row(iNewObs2)));
			}
		  }else{
			if(exportHazard){
			  iStrata_IFhazard = iStrata_IFhazard0;
			}
			if(exportCumhazard || exportSurvival){
			  iStrata_IFcumhazard = iStrata_IFcumhazard0;
			}
		  }
		  if(exportSurvival){
			iStrata_IFsurvival = -iStrata_IFcumhazard*newSurvival(iNewObs2,iTauStore);
		  }
		  
		  // store
		  if(exportIF){
			if(exportHazard){IF_hazard.slice(iTauStore).col(iNewObs2)= iStrata_IFhazard;}
			if(exportCumhazard){IF_cumhazard.slice(iTauStore).col(iNewObs2) = iStrata_IFcumhazard;}
			if(exportSurvival){IF_survival.slice(iTauStore).col(iNewObs2) = iStrata_IFsurvival;}
		  }
		  
		  if(exportSE){
			if(exportCumhazard){SE_cumhazard(iNewObs2,iTauStore) = arma::sum(iStrata_IFcumhazard % iStrata_IFcumhazard);}
			if(exportSurvival){SE_survival(iNewObs2,iTauStore) = arma::sum(iStrata_IFsurvival % iStrata_IFsurvival);}
		  }
			
		  if(exportIFmean && diag){
			for(int iFactor=0; iFactor<nFactor; iFactor++){
			  
			  if(exportHazard){IFmean_hazard[iFactor].col(iTauStore) += iStrata_IFhazard * factor[iFactor](iNewObs2,iTauStore);}
			  if(exportCumhazard){IFmean_cumhazard[iFactor].col(iTauStore) += iStrata_IFcumhazard * factor[iFactor](iNewObs2,iTauStore);}
			  if(exportSurvival){IFmean_survival[iFactor].col(iTauStore) += iStrata_IFsurvival * factor[iFactor](iNewObs2,iTauStore);}
			}
		  }
		}
	  }

	  // **** IF mean hazard/cumhazard/survival
	  if(exportIFmean && diag == false){
		// <IF>(hazard) = E[w * eXb] IF_hazard0 + E[w * eXb * X] * hazard0 * IF_beta
		// <IF>(cumhazard) = E[w * eXb] IF_cumhazard0 + E[w * eXb * X] * cumhazard0 * IF_beta
		// <IF>(survival) = -(E[w * Surv * eXb] IF_cumhazard0 + E[w * Surv * eXb * X] * cumhazard0 * IF_beta)

		if(debug>1){Rcpp::Rcout << " IF mean ";}
		tempo_uvec(0) = iTau;
		
		for(int iFactor=0; iFactor<nFactor; iFactor++){

		  iStrata_weXb = factor[iFactor].submat(newdata_index[iStrata],tempo_uvec) % neweXb(newdata_index[iStrata]);
		  if(p>0){
			iStrata_weXbX = newX.rows(newdata_index[iStrata]);
			iStrata_weXbX.each_col() %= iStrata_weXb;
		  }

		  if(exportHazard){
			IFmean_hazard[iFactor].col(iTauStore) += iStrata_IFhazard0 * arma::sum(iStrata_weXb);
			if(p>0){
			  IFmean_hazard[iFactor].col(iTauStore) += IFbeta % arma::trans(arma::sum(iStrata_weXbX,0)) * hazard0[iStrata](iTau);
			}
		  }
		
		  if(exportCumhazard){
			IFmean_cumhazard[iFactor].col(iTauStore) += iStrata_IFcumhazard0 * arma::sum(iStrata_weXb);
			if(p>0){
			  IFmean_cumhazard[iFactor].col(iTauStore) += IFbeta * arma::trans(arma::sum(iStrata_weXbX,0)) * cumhazard0[iStrata](iTau);
			}
		  }

		  if(exportSurvival){
			iStrata_weXbS = iStrata_weXb % newSurvival.submat(newdata_index[iStrata],tempo_uvec);
			IFmean_survival[iFactor].col(iTauStore) -= iStrata_IFcumhazard0 * arma::sum(iStrata_weXbS,0);

			if(p>0){
			  iStrata_weXbXS = iStrata_weXbX;
			  iStrata_weXbXS.each_col() %= newSurvival.submat(newdata_index[iStrata],tempo_uvec);
			  IFmean_survival[iFactor].col(iTauStore) -= IFbeta * arma::trans(arma::sum(iStrata_weXbXS,0)) * cumhazard0[iStrata](iTau);
			}
		  }
			
		} // end IFactor
	  } // end if
	} // end iTime
	if(debug>1){Rcpp::Rcout << std::endl;}
  } // end iStrata

  // ** Post process
  if(debug>0){Rcpp::Rcout << "Post process" << std::endl;}
  if(exportSE){
	if(exportCumhazard){SE_cumhazard = sqrt(SE_cumhazard);}
	if(exportSurvival){SE_survival = sqrt(SE_survival);}
  }

  if(exportIFmean == true){
	for(int iFactor=0; iFactor<nFactor; iFactor++){
	  if(exportHazard){IFmean_cumhazard[iFactor] /= nNewObs;}
	  if(exportCumhazard){IFmean_cumhazard[iFactor] /= nNewObs;}
	  if(exportSurvival){IFmean_survival[iFactor] /= nNewObs;}
	}
  }

  // ** Export
  if(debug>0){Rcpp::Rcout << "Export" << std::endl;}
  return(List::create(Named("IF_hazard") = IF_hazard,
					  Named("IFmean_hazard") = IFmean_hazard,
					  Named("IF_cumhazard") = IF_cumhazard,
					  Named("SE_cumhazard") = SE_cumhazard,
					  Named("IFmean_cumhazard") = IFmean_cumhazard,
					  Named("IF_survival") = IF_survival,
					  Named("SE_survival") = SE_survival,
					  Named("IFmean_survival") = IFmean_survival));
}

// * calcAIFsurv_cpp: compute average IF for the cumlative hazard / survival (method 3)
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
