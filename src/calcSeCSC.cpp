// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


// * calcSeMinimalCSC_cpp: compute IF/sumIF/se for the cif (method 1)
// [[Rcpp::export]]
// J: number of jump times
// n: number of observations in the training set
// N: number of observations in the prediction set
// p: number of regressors
// S: number of strata
// T: number of prediction times
List calcSeMinimalCSC_cpp(const arma::vec& seqTau, // horizon time for the predictions (T)
						  const arma::mat& newSurvival, // predicted survival for all observations at each jump time (NxJ)
						  const arma::mat& hazard0, // baseline hazard of the event of interest for each strata JxS
						  const std::vector< arma::mat >& cumhazard0, // baseline cumulative hazard for each strata and cause S:(JxC)
						  const std::vector< arma::mat >& newX, // design matrix for each cause C:(Nxp)
						  const arma::mat& neweXb, // exponential of the linear predictor for each cause (NxC)
						  const std::vector< arma::mat >& IFbeta, //influence function for the regression coefficients for each cause C:(nxp)
						  const std::vector< arma::mat >& Ehazard0, // to compute the influence function of the baseline hazard S:(pxJ)
						  const std::vector< std::vector< arma::mat > >& cumEhazard0, // to compute the influence function of the baseline hazard C:S:(pxJ)
						  const std::vector< arma::vec >& hazard_iS0, // to compute the influence function of the baseline hazard S:J
						  const std::vector< std::vector< arma::vec > >& cumhazard_iS0, // to compute the influence function of the baseline hazard C:S:J
						  const std::vector< arma::mat> & delta_iS0, // to compute the influence function of the baseline hazard C:(nxS, S because multiplied by the strata indicator)
						  const std::vector< arma::mat> & sample_eXb, // to compute the influence function of the baseline hazard C:(nxS, S because multiplied by the strata indicator)
						  const arma::vec& sample_time, // event times of the training set (n)
						  const std::vector< std::vector< arma::uvec > > & indexJumpSample_time, // index of the jump time corresponding to the sample times n C:S:J
						  const arma::vec& jump_time, // jump times for each strata J
						  const arma::mat& isJump_time1, // jump times for each strata J
						  const std::vector< std::vector< arma::vec > > & jump2jump, // index of the jump time in each strata/cause corresponding to the jump time over all strata C:S:J
						  const arma::vec& firstTime1theCause, // time before which the prediction are 0 (S)
						  const arma::vec& lastSampleTime, // time after which the prediction are NA (S)
						  const std::vector< arma::uvec >& newdata_index, // index of the observations within strata
						  const std::vector< arma::mat >& factor,
						  const arma::mat& grid_strata,
						  int nTau, int nNewObs, int nSample, int nStrata, int nCause, const arma::vec& p, 
						  int theCause, bool diag, bool survtype,
						  bool exportSE, bool exportIF, bool exportIFmean,
						  int debug){

  // ** prepare
  if(debug>0){Rcpp::Rcout << "Prepare" << std::endl;}
  int iStrata_tau,iStrata_tau2,iStrata_tauMax; 
  arma::vec iStrata_seqTau;
  arma::uvec iStrata_indexJumpTau;
  int iStrata_nNewObs, iiJump;
  int nJump = jump_time.size();
  int iTauStore = 0; // factice initialization to avoid warning
  int nFactor = factor.size();
  
  int iNewObs2;
  arma::colvec iStrata_IFhazard0;
  std::vector< arma::vec > iStrata_IFcumhazard0(nCause);
  arma::colvec iStrata_IFint;
  arma::mat iStrata_AIFint;
  double iSlambda1;
  
  arma::uvec index_timestop(nSample);  
  arma::uvec tempo_uvec(1), tempo_uvec2(1);

  arma::vec iStrata_wSeXb1,iStrata_wSeXb1eXbj;
  arma::mat iStrata_wSeXb1X1,iStrata_wSeXb1eXbjXj;
  
  // ** initialize
  if(debug>0){Rcpp::Rcout << "Initialize" << std::endl;}
  arma::cube IF_cif;
  arma::mat SE_cif;
  std::vector< arma::mat > IFmean_cif(nFactor);
  if(exportIF || exportSE){
	IF_cif.resize(nSample, nNewObs, nTau);
	IF_cif.fill(0.0);
  }
  if(exportSE){
	SE_cif.resize(nNewObs, nTau);
	SE_cif.fill(0.0);
  }
  if(exportIFmean){
	for(int iFactor=0; iFactor<nFactor; iFactor++){
	  IFmean_cif[iFactor].resize(nSample, nTau);
	  IFmean_cif[iFactor].fill(0.0);
	}
  }
  
  // ** compute IF within each strata
  if(debug>0){Rcpp::Rcout << "Compute IF" << std::endl;}
  for(int iStrata=0; iStrata<nStrata ; iStrata++){
	if(debug>0){Rcpp::Rcout << " > strata " << iStrata << "/" << (nStrata-1) << " ";}

	// *** check if any observation in strata
	iStrata_nNewObs = newdata_index[iStrata].size();
	if(iStrata_nNewObs==0){continue;}

	// *** narrow down prediction times
    iStrata_tau=0;
    if(diag){
	  iStrata_tauMax = iStrata_nNewObs-1;
	  iStrata_seqTau = seqTau(newdata_index[iStrata]);
	  iStrata_nNewObs = 1;
	}else{
	  iStrata_tauMax = nTau-1;
	  iStrata_seqTau = seqTau;
	}

	// WARNING: this part needs to be before iStrata_tauMax is modified (i.e. the next while loop)
	while((iStrata_tau <= iStrata_tauMax) && firstTime1theCause(iStrata)>iStrata_seqTau(iStrata_tau)){ // start at the first event or after 
	  iStrata_tau++;    
	}
	if(iStrata_tau > iStrata_tauMax){continue;}

	while((iStrata_tauMax >= 0) && seqTau(iStrata_tauMax)>lastSampleTime(iStrata)){ // end at the last event or before
	  for(int iNewObs=0; iNewObs<iStrata_nNewObs ; iNewObs++){		
		if(diag){
		  iTauStore = 0;
		  iNewObs2 = newdata_index[iStrata](iStrata_tauMax);
		}else{
		  iTauStore = iStrata_tauMax;
		  iNewObs2 = newdata_index[iStrata](iNewObs);
		}
		
		if(exportIF || exportSE){IF_cif.slice(iTauStore).col(iNewObs2).fill(NA_REAL);}
		if(exportIFmean){
		  for(int iFactor=0; iFactor<nFactor; iFactor++){
			IFmean_cif[iFactor].col(iTauStore).fill(NA_REAL);
		  }
		}
	  }
	  iStrata_tauMax--;
	}
	if(iStrata_tauMax < 0){continue;}
	
    R_CheckUserInterrupt();
	iStrata_tau2=iStrata_tau;
	if(exportIFmean && diag == false){
	  iStrata_AIFint.resize(nSample,nFactor);
	  iStrata_AIFint.fill(0.0);
	}
	if(debug>1){Rcpp::Rcout << " (tau=" << iStrata_tau << "-" << iStrata_tauMax << ") " << endl;}
	
	// *** compute IF/SE/IFmean at each time point
	for(int iJump=0; iJump<nJump; iJump++){
	  if(debug>1){Rcpp::Rcout << " IF0 ";}

	  // **** IF baseline hazard
	  if(isJump_time1(iJump,iStrata)){ // only update for jumps corresponding to the event of interest in the strata
		iiJump = jump2jump[theCause][grid_strata(iStrata,theCause)](iJump);
	  
		iStrata_IFhazard0 = delta_iS0[theCause].col(grid_strata(iStrata,theCause)) % (sample_time == jump_time(iJump)) - sample_eXb[theCause].col(grid_strata(iStrata,theCause)) % (jump_time(iJump) <= sample_time) * hazard_iS0[iStrata](iiJump);

		if(p(theCause)>0){
		  iStrata_IFhazard0 -= IFbeta[theCause] * Ehazard0[iStrata].col(iJump);
		}
	  }

	  // **** IF baseline cumulative hazard hazard
	  if(isJump_time1(iJump,iStrata) && iJump>0){
		for(int iCause=0; iCause<nCause; iCause++){
		  iiJump = jump2jump[iCause][grid_strata(iStrata,theCause)](iJump-1);

		  index_timestop = indexJumpSample_time[iCause][grid_strata(iStrata,theCause)];
		  index_timestop.elem(find(index_timestop > iiJump)).fill(iiJump);

		  iStrata_IFcumhazard0[iCause] = delta_iS0[iCause].col(grid_strata(iStrata,iCause)) % (sample_time <= jump_time(iJump-1)) - sample_eXb[iCause].col(grid_strata(iStrata,iCause)) % cumhazard_iS0[iCause][grid_strata(iStrata,theCause)](index_timestop);
		  if(p(iCause)>0){
			iStrata_IFcumhazard0[iCause] -= IFbeta[iCause] * cumEhazard0[iCause][grid_strata(iStrata,iCause)].col(iJump-1);
		  }
		}
	  }
	  
	  // **** IF/SE cif
	  if(exportIF || exportSE || (exportIFmean && diag)){ 
		if(debug>1){Rcpp::Rcout << " IF " ;}

		if(isJump_time1(iJump,iStrata)){ // only update for jumps corresponding to the event of interest in the strata
		  for(int iNewObs=0; iNewObs<iStrata_nNewObs; iNewObs++){
			iNewObs2 = newdata_index[iStrata](iNewObs);
			if(diag){
			  iTauStore = 0;
			  if(jump_time[iJump]>seqTau[iNewObs2]){continue;}
			}else{
			  iTauStore = iStrata_tau;
			}

		 	iSlambda1 = newSurvival(iNewObs2,iJump) * hazard0(iJump,grid_strata(iStrata,theCause)) * neweXb(iNewObs2,theCause);

			iStrata_IFint = newSurvival(iNewObs2,iJump) * iStrata_IFhazard0 * neweXb(iNewObs2,theCause);
			if(p(theCause)>0){
			  iStrata_IFint += iSlambda1 * IFbeta[theCause] * arma::trans(newX[theCause].row(iNewObs2));
			}
			if(iJump>0){
			  for(int iCause=0; iCause<nCause; iCause++){
				iStrata_IFint -= iSlambda1 * iStrata_IFcumhazard0[iCause] * neweXb(iNewObs2,iCause);
				if(p(iCause)>0){
				  iiJump = jump2jump[iCause][grid_strata(iStrata,theCause)](iJump-1);
				  iStrata_IFint -= iSlambda1 * cumhazard0[iCause](iJump-1,grid_strata(iStrata,iCause)) * neweXb(iNewObs2,iCause) * IFbeta[iCause] * arma::trans(newX[iCause].row(iNewObs2));
				} 
			  }
			}
			
			// store		  
			if(exportIF || exportSE){
			  IF_cif.slice(iTauStore).col(iNewObs2) += iStrata_IFint;
			}
			if(exportIFmean && diag){
			  for(int iFactor=0; iFactor<nFactor; iFactor++){
				IFmean_cif[iFactor].col(iTauStore) += iStrata_IFint * factor[iFactor](iNewObs2,iJump);
			  }
			}
		  }
		}

		while((iStrata_tau <= iStrata_tauMax) && ( (((iJump+1)<nJump) && (seqTau[iStrata_tau] < (jump_time[iJump+1]))) || (iJump+1==nJump))){
		  iStrata_tau++;
		  if(iStrata_tau <= iStrata_tauMax){
			if(exportIF || exportSE){
			  IF_cif.slice(iStrata_tau) = IF_cif.slice(iStrata_tau-1);
			}
		  }
		}
	  }

	  // **** IF mean hazard/cumhazard/survival
	  if(exportIFmean && diag == false){
		// <IF>(cif) = E[w * Surv * eXb1] IF_hazard01 + E[w * Surv * eXb1 * X1] * hazard01 * IF_beta1
		//             - \sum_j (E[w * Surv * eXb1 * eXbj] * hazard01 * IF_cumhazard0j + E[w * Surv * eXb1 * eXbj * Xj] * hazard01 * cumhazard0j * IF_betaj)

		if(debug>1){Rcpp::Rcout << " IF mean ";}
		tempo_uvec(0) = iJump;
		
		for(int iFactor=0; iFactor<nFactor; iFactor++){

		  if(isJump_time1(iJump,iStrata)){ // only update for jumps corresponding to the event of interest in the strata
			tempo_uvec2(0) = theCause;
			iStrata_wSeXb1 = factor[iFactor].submat(newdata_index[iStrata],tempo_uvec) % newSurvival.submat(newdata_index[iStrata],tempo_uvec) % neweXb(newdata_index[iStrata],tempo_uvec2);
			iStrata_AIFint.col(iFactor) += iStrata_IFhazard0 * arma::sum(iStrata_wSeXb1);

			if(p(theCause)>0){
			  iStrata_wSeXb1X1 = newX[theCause].rows(newdata_index[iStrata]);
			  iStrata_wSeXb1X1.each_col() %= iStrata_wSeXb1;
			  iStrata_AIFint.col(iFactor) += IFbeta[theCause] * arma::trans(arma::mean(iStrata_wSeXb1X1,0)) * hazard0(iJump,grid_strata(iStrata,theCause));
			}

			if(iJump>0){
			for(int iCause=0; iCause<nCause; iCause++){
			  tempo_uvec2(0) = iCause;
			  iStrata_wSeXb1eXbj = iStrata_wSeXb1 % neweXb(newdata_index[iStrata],tempo_uvec2);
			  iStrata_AIFint.col(iFactor) -= iStrata_IFcumhazard0[iCause] * arma::sum(iStrata_wSeXb1eXbj) * hazard0(iJump,grid_strata(iStrata,theCause));

			  if(p(iCause)>0){
			  	iStrata_wSeXb1eXbjXj = newX[iCause].rows(newdata_index[iStrata]);
				iStrata_wSeXb1eXbjXj.each_col() %= iStrata_wSeXb1eXbj;
				iStrata_AIFint.col(iFactor) -= IFbeta[iCause] * arma::trans(arma::mean(iStrata_wSeXb1eXbjXj,0)) * hazard0(iJump,grid_strata(iStrata,theCause)) * cumhazard0[iCause](iJump-1,grid_strata(iStrata,iCause));
			  }
			}
			}
		  }

		  while((iStrata_tau2 <= iStrata_tauMax) && ( (((iJump+1)<nJump) && (seqTau[iStrata_tau2] < jump_time[iJump+1])) || (iJump+1==nJump))){
			IFmean_cif[iFactor].col(iStrata_tau2) += iStrata_AIFint.col(iFactor);
			iStrata_tau2++;
		  }
		
		} // end IFactor
	  } // end if

	  if(debug>1){Rcpp::Rcout << " end " ;}
      if((iStrata_tau > iStrata_tauMax) || (iStrata_tau2 > iStrata_tauMax)){break;}

	} // end iJump
	if(debug>1){Rcpp::Rcout << std::endl;}
  } // end iStrata

  // ** Post process
  if(debug>0){Rcpp::Rcout << "Post process" << std::endl;}

  if(exportSE){
	for(int iTau=0; iTau<nTau; iTau++){	
	  SE_cif.col(iTau) = arma::trans(sqrt(sum(IF_cif.slice(iTau) % IF_cif.slice(iTau), 0)));
	}
  }

  if(exportIFmean == true){
	for(int iFactor=0; iFactor<nFactor; iFactor++){
	  IFmean_cif[iFactor] /= nNewObs;
	}
  }

  // ** Export
  if(debug>0){Rcpp::Rcout << "Export" << std::endl;}
  return(List::create(Named("IF_cif") = IF_cif,
					  Named("SE_cif") = SE_cif,
					  Named("IFmean_cif") = IFmean_cif));
}


// * calcSeCif2_cpp: compute IF for the absolute risk (method 2)
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




