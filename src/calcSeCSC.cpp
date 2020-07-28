// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


// * calcSeCif_cpp: compute IF for the absolute risk (method 1)
// [[Rcpp::export]]
List calcSeCif_cpp(const arma::vec& seqTau, // horizon time for the predictions
				   const arma::mat& newSurvival,
				   const arma::mat& newHazard0,
				   const std::vector< arma::mat > & newCumHazard0,
				   const std::vector< arma::mat >& newX,
				   const arma::mat& neweXb,
				   const std::vector< arma::mat >& IFbeta,
				   const std::vector< arma::mat >& Ehazard0,
				   const std::vector< std::vector < arma::vec > >& cumEhazard0,
				   const std::vector< arma::vec >& hazard_iS0,
				   const std::vector< std::vector < arma::vec > >& cumhazard_iS0,
				   const std::vector< arma::mat >& delta_iS0,
                   const std::vector< arma::mat >& sample_eXb,
				   const arma::vec& sample_time,
				   const arma::vec& jumpTime,
				   const std::vector< arma::vec>& jumpTheCause,
				   const arma::vec& lastSampleTime,
				   const std::vector< arma::uvec>& newdata_index,
				   const arma::mat& grid_strata,
				   int theCause, 
				   int nTau, int nNewObs, int nSample, int nCause, int nStrata, const IntegerVector& p,
				   bool survtype,
				   bool exportSE, bool exportIF, bool exportIFsum){

  // ** prepare
  Rcpp::Rcout << "start" << std::endl;
  
  // temporary store influence function at time s
  int iStrata_tauMin,iStrata_tauMax,iStrata_tau; 
  NumericVector iStrata_indexJump; // position of the jump times relative to the cause of interest and the current strata among all event times
  int iStrata_nJump; // number of jump times relative to the cause of interest in the current strata
  int iStrata_nNewObs;
  int iiJump;
  int nJump = jumpTime.size();
  arma::mat iStrata_indicator(nSample,nCause);
  arma::colvec iStrata_IFhazard01; // influence function of the baseline hazard (cause of interest)
  arma::colvec iStrata_IFcumhazard012; // influence function of the baseline cumulative hazard (all causes)
  arma::mat iStrata_IFrisk; // influence function of the baseline cumulative hazard (all causes)

  arma::colvec tempo_colvec;
  arma::uvec tempo_uvec;
    
  arma::mat SErisk;
  if(exportSE){
    SErisk.resize(nNewObs, nTau);
    SErisk.fill(0.0);
  }  
  arma::cube IFrisk;
  if(exportIF){
    IFrisk.resize(nNewObs, nSample, nTau);
    IFrisk.fill(0.0);
  }
  arma::mat IFsumrisk;
  if(exportIFsum){
    IFsumrisk.resize(nSample, nTau);
    IFsumrisk.fill(0.0);
  }

  // ** compute IF at interesting times
  for(int iStrata=0; iStrata<nStrata ; iStrata++){
	Rcpp::Rcout << "iStrata " << iStrata << std::endl;
 
	// *** check if any observation in strata
	iStrata_nNewObs = newdata_index[iStrata].size();
	if(iStrata_nNewObs==0){continue;}
	iStrata_IFhazard01.resize(nSample);
	iStrata_IFcumhazard012.resize(nSample);
	iStrata_IFrisk.resize(nSample,iStrata_nNewObs);
	iStrata_IFrisk.fill(0.0);
	  
	// *** jump times
	iStrata_indexJump = jumpTheCause[iStrata];
	iStrata_nJump = iStrata_indexJump.size();
	if(iStrata_nJump==0){continue;}
	
	// *** narrow down prediction times
	iStrata_tauMin=0;
	iStrata_tauMax=nTau-1;
	
	while((iStrata_tauMin < nTau) && jumpTime[iStrata_indexJump[0]]>seqTau[iStrata_tauMin]){ // start at the first event or after 
	  iStrata_tauMin++;    
	}
	while((iStrata_tauMax >= 0) && seqTau[iStrata_tauMax]>lastSampleTime[iStrata]){ // end at the last event or before
	  for(int iNewObs=0; iNewObs<iStrata_nNewObs ; iNewObs++){
		SErisk(newdata_index[iStrata][iNewObs],iStrata_tauMax-1) = NA_REAL;
		IFrisk.slice(iStrata_tauMax-1).col(newdata_index[iStrata][iNewObs]).fill(NA_REAL);
	  }
	  IFsumrisk.col(iStrata_tauMax-1).fill(NA_REAL);
	  iStrata_tauMax--;
	}
	if(iStrata_tauMin >= nTau || iStrata_tauMax < 0){continue;}
	iStrata_tau = iStrata_tauMin;
	
    R_CheckUserInterrupt();
		
	// *** compute IF over jump times
	for(int iJump=0; iJump<iStrata_nJump; iJump++){
	  Rcpp::Rcout << "iJump " << iJump << " / " << nJump << std::endl;
	  iiJump = iStrata_indexJump[iJump]; // iJump is relative to the jump times in the strata and cause, while iiJump is relative to all jump times

	  // **** contribution IF\beta
	  for(int iCause=0; iCause<nCause; iCause++){
		if(p[iCause]>0){
		  for(int iNewObs=0; iNewObs<iStrata_nNewObs ; iNewObs++){
			iStrata_IFrisk.col(iNewObs) -= newSurvival(newdata_index[iStrata][iNewObs],iiJump) * newHazard0(iiJump,iStrata) * neweXb(newdata_index[iStrata][iNewObs],theCause) * newCumHazard0[iCause](iiJump,iStrata) * neweXb(newdata_index[iStrata][iNewObs],iCause) * IFbeta[iCause] * arma::trans(newX[iCause].row(newdata_index[iStrata][iNewObs]));
		  }
		}
	  }
	  if(p[theCause]>0){
		for(int iNewObs=0; iNewObs<iStrata_nNewObs ; iNewObs++){
		  iStrata_IFrisk.col(iNewObs) += newSurvival(newdata_index[iStrata][iNewObs],iiJump) * newHazard0(iiJump,iStrata) * neweXb(newdata_index[iStrata][iNewObs],theCause) * IFbeta[theCause] * arma::trans(newX[theCause].row(newdata_index[iStrata][iNewObs]));
		}
	  }
	  Rcpp::Rcout << "endIFbeta" << std::endl;
	  
	  // **** IF\lambda01
	  iStrata_IFhazard01 = - sample_eXb[0].col(iStrata) * hazard_iS0[iStrata][iiJump] + (jumpTime[iiJump] == sample_time) % delta_iS0[0].col(iStrata);
	  if(p[0]>0){
		iStrata_IFhazard01 -= IFbeta[0] * Ehazard0[iStrata];
	  }

	  // NOTE: strata indicator in eXb and delta_iS0		  
	  for(int iNewObs=0; iNewObs<iStrata_nNewObs ; iNewObs++){
		iStrata_IFrisk.col(iNewObs) = newSurvival(newdata_index[iStrata][iNewObs],iiJump) * iStrata_IFhazard01 * neweXb(newdata_index[iStrata][iNewObs],theCause);
	  }
	  Rcpp::Rcout << "endIFlambda01" << std::endl;

	  // **** compute IF\Lambda01,IF\Lambda02
	  if(iiJump>0){
		for(int iCause=0; iCause<nCause; iCause++){
		  iStrata_IFcumhazard012 = -sample_eXb[iCause].col(iStrata) * cumhazard_iS0[iCause][iStrata][iiJump-1] + (jumpTime[iiJump-1] >= sample_time) % delta_iS0[iCause].col(iStrata);
		  if(p[iCause]>0){
			iStrata_IFcumhazard012 -= IFbeta[iCause] * cumEhazard0[iCause][iStrata];
		  }
		  // NOTE: strata indicator in eXb and delta_iS0		  
		  for(int iNewObs=0; iNewObs<nNewObs ; iNewObs++){
			iStrata_IFrisk.col(iNewObs) -= newSurvival(0,iiJump) * newHazard0(iiJump,iStrata) * neweXb(iNewObs,iCause) * iStrata_IFcumhazard012 * neweXb(iNewObs,iCause);
		  }
		}
	  }
	  Rcpp::Rcout << "endIFLambda012" << std::endl;
	  
	  // **** export
	  Rcpp::Rcout << iStrata_tau << "/" << iStrata_tauMax << " | " << iiJump+1 << "/" << nJump << " | " << seqTau[iStrata_tau] << "/" << jumpTime[iiJump+1] << std::endl;
	  while((iStrata_tau <= iStrata_tauMax) && ( (((iJump+1)<iStrata_nJump) && (seqTau[iStrata_tau] < jumpTime[iiJump+1])) || (iJump+1==iStrata_nJump))){
	  Rcpp::Rcout << "*" << std::endl;

	  	for(int iSample=0; iSample<nSample ; iSample++){
	  	  if(exportSE){
			for(int iNewObs=0; iNewObs<iStrata_nNewObs ; iNewObs++){
			  SErisk(newdata_index[iStrata][iNewObs],iStrata_tau) = arma::mean(iStrata_IFrisk.col(iNewObs) % iStrata_IFrisk.col(iNewObs));
			}
	  	  }
	  	  if(exportIF){
	  		for(int iNewObs=0; iNewObs<iStrata_nNewObs ; iNewObs++){			  
	  		  IFrisk.slice(iStrata_tau).col(newdata_index[iStrata][iNewObs]) = iStrata_IFrisk.col(iNewObs);
	  		}
	  	  }
	  	  if(exportIFsum){
			Rcpp::Rcout << "b" << std::endl;
	  		IFsumrisk.col(iStrata_tau) += arma::mean(iStrata_IFrisk,1);
	  	  }
	  	}
	  	iStrata_tau++;	  
	  }
	}
	
  }
  
  // export
  return(List::create(Named("iid") = IFrisk,
					  Named("iidsum") = IFsumrisk,
					  Named("se") = SErisk));
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

// * calcAIFcif_cpp: compute average IF for the absolute risk (method 3)
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



