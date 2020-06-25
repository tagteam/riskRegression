// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// * calcAIFsurv_cpp: compute average IF for the cumlative hazard / survival
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


// * calcSeCif_cpp: compute IF for the absolute risk (method 1)
// [[Rcpp::export]]
List calcSeCif_cpp(const std::vector<double>& tau,
		   const std::vector<std::vector<arma::mat>>& zipIF_Elambda0,
		   const std::vector<std::vector<arma::mat>>& zipIF_cumElambda0,
		   const std::vector<arma::vec>& zipIF_strata,
		   const std::vector<std::vector<arma::vec>>& zipIF_eXb,
		   const std::vector<std::vector<arma::vec>>& zipIF_time,
		   const std::vector<std::vector<arma::vec>>& zipIF_jump,
		   const std::vector<std::vector<arma::vec>>& zipIF_lambda0_iS0,
		   const std::vector<std::vector<arma::vec>>& zipIF_cumLambda0_iS0,
		   const std::vector<std::vector<arma::vec>>& zipIF_delta_iS0,
		   const std::vector<std::vector<arma::vec>>& zipIF_time1,
		   const std::vector<arma::mat>& ls_IFbeta,
		   const std::vector<arma::mat>& ls_X,
		   const std::vector<arma::mat>& ls_cumhazard,
		   const arma::mat& ls_hazard,
		   const arma::mat& survival,
		   const arma::mat& eXb,
		   const arma::mat& strata,
		   const arma::vec& strataU,
		   const std::vector<arma::uvec>& indexStrata,
		   int theCause, bool hazardType,
		   int nNewObs, int nSample, int nStrata, arma::vec nVar, int nCause, 
		   bool exportSE, bool exportIF, bool exportIFsum, bool diag){

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
      outIF.resize(nNewObs,nSample,1);
    }else{
      outIF.resize(nNewObs,nSample,nTau);
    }
    outIF.fill(0.0);
  }
  arma::mat outIFsum;
  if(exportIFsum){
    if(diag){
      outIFsum.resize(nSample,1);
    }else{
      outIFsum.resize(nSample,nTau);
    }
    outIFsum.fill(0.0);
  }

  // ** initialize variables
  int indexTau_S; // index of the current prediction time

  arma::uvec indexNewObs_S; // prediction observations in the current strata
  int nNewObs_S;
    
  arma::uvec indexOldObs_S;  // training observations in the current strata
  arma::vec oldObsTime_S;
  int nOldObs_S;

  arma::uvec iVec_indexJump_SJ(nCause);
  arma::vec nJump_S(nCause); 
  double currentTime;

  arma::colvec IFhazard_SJ(nSample);
  arma::mat IFcumhazard_SJ(nSample,nCause);

  arma::mat iIFcif_SJ;
  
  // ** loop over strata
  for(int iS = 0; iS<nStrata; iS++){

    // *** prepare
    // position of the first remaining of the prediction time
    indexTau_S = 0; 

    // IF for the CIF at the current jump time
    iIFcif_SJ.resize(nSample,nNewObs_S); 
    iIFcif_SJ.fill(0.0);

    // prediction sample
    indexNewObs_S = indexStrata[iS]; // set of individuals in the strata
    nNewObs_S = indexNewObs_S.size(); // number of individuals in the strata
    
    // training sample from the same strata
    indexOldObs_S = arma::find(zipIF_strata[0] == strata(iS,0));
    oldObsTime_S = zipIF_time[0][[iS]](indexOldObs_S);
    nOldObs_S = indexOldObs_S.size(); 
    
    for(int iC=0; iC<nCause; iC++){
      nJump_S(iC) = zipIF_time1[iC][iS].size();  // number of jumps 
    }
    iVec_indexJump_SJ.fill(-1); // index of the jump for each cause
    
    
    // *** no jump: go to the next strata
    if(nJump_S(0) == 0){
      if(exportSE){
	outSE.rows(nNewObs).fill(0.0);
      }
      if(exportIF){
	outIF.rows(nNewObs).fill(0.0);
      }
      // not need to update IFsum since the contribution is 0
      continue;
    }
    // *** otherwise move to or just before the first jump
    currentTime = zipIF_time1[0][iS](0);
    while(indexTau_S < nTau && tau[indexTau_S] < currentTime){
      if(exportSE){
	if(diag){
	  outSE.rows(nNewObs).col(0).fill(0.0);
	}else{
	  outSE.rows(nNewObs).col(indexTau_S).fill(0.0);
	}
      }
      if(exportIF){
	if(diag){
	  outIF.slice(0).rows(nNewObs).fill(0.0);
	}else{
	  outIF.slice(indexTau_S).rows(nNewObs).fill(0.0);
	}
      }
      // not need to update IFsum since the contribution is 0
      indexTau_S++;
    }

    // *** loop over jump times
    for(int iJ = 0; iS<nJump_S(0); iS++){
      currentTime = zipIF_time1[0][iS](iJ);
      
      // **** compute IFhazard and IFcumhazard
      IFhazard_SJ.fill(0.0);
      IFcumhazard_SJ.fill(0.0);

      // IFcumhazard at t-
      for(int iC=0; iC<nCause; iC++){

	// find jump corresponding to t-
	while(((iVec_indexJump_SJ(iC)+1) < nJump_S(iC)) & (zipIF_time1[iC][iS](iVec_indexJump_SJ(iC)+1) < currentTime)){
	  iVec_indexJump_SJ(iC)++;
	}

	// first term
	if(nVar[iC]>0){
	  if(iC==0){IFhazard_SJ -= ls_IFbeta[0] * zipIF_Elambda0[0][iS].col(iJ);}
	  IFcumhazard_SJ.col(iC) -= ls_IFbeta[iC] * zipIF_cumElambda0[iC][iS].col(iVec_indexJump_SJ(iC));
	}

	// second and third term
	for(int iOldObs=0; iOldObs < nOldObs_S; iOldObs++){
	  if(iC==0 & currentTime <= oldObsTime_S(iOldObs)){IFhazard_SJ(indexOldObs_S(iOldObs)) -= zipIF_lambda0_iS0[0][iS](iJ) * zipIF_eXb[0][iS];}
	  if(iC==0 & currentTime == oldObsTime_S(iOldObs)){IFhazard_SJ(indexOldObs_S(iOldObs)) += zipIF_delta_iS0[0][iS];}
	  
	  iMin = min(zipIF_jump[iC][iS](iIndicator_strata[iObs]),iVec_indexJump_SJ(iC)); 
	  IFcumhazard_SJ(indexOldObs_S(iOldObs),iC) -= zipIF_cumLambda0_iS0[iC][iS](min()) * zipIF_eXb[iC][iS](indexOldObs_S(iOldObs));
	  IFcumhazard_SJ(indexOldObs_S(iOldObs),iC) -= zipIF_delta_iS0[iC][iS] % (currentTime <= zipIF_time[iC][[iS]]);
	}
      }
      // indexOldObs_S = arma::find(zipIF_strata[0] == strata(iS,0));
      // nOldObs_S = indexOldObs_S.size(); 
      
      
      // **** loop over individuals
      for(int iNewObs=0; iNewObs < nNewObs_S; iNewObs++){
	// compute IF
	
      }

      // store
      while((indexTau_S < nTau) && (jump_S[iJ] <= tau[indexTau_S]) && ((iJ < (nTau-1)) && (tau[indexTau_S] <= jump_S[iJ+1])) || ((iJ == nTau-1) && (jump_S[iJ] == tau[indexTau_S])) ){

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

      iIFcif_SJ
   }

  // *** fill remaining columns with NA
    while(iTau < nTau){
      if(exportSE){
	if(diag){
	  outSE.rows(iIndexNewObs).col(0).fill(NA_REAL);
	}else{
	  outSE.rows(iIndexNewObs).col(iTau).fill(NA_REAL);
	}
      }
      if(exportIF){
	if(diag){
	  outIF.slice(0).rows(iIndexNewObs).fill(NA_REAL);
	}else{
	  outIF.slice(iTau).rows(iIndexNewObs).fill(NA_REAL);
	}
      }
      if(exportIFsum){
	if(diag){
	  outIFsum.col(0).fill(NA_REAL);
	}else{
	  outIFsum.col(iTau).fill(NA_REAL);
	}
      }
      iTau++;
    }
  }
  
  // ** export
  return(List::create(Named("se") = outSE,
		      Named("iid") = outIF,
		      Named("average.iid") = outIFsum));
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



