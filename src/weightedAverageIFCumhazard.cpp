// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// J: number of jump times
// n: number of observations in the training set
// N: number of observations in the prediction set
// p: number of regressors
// S: number of strata
// T: number of prediction times
// should be able to speed this part up by using the trick with firsthit
// [[Rcpp::export]]
NumericVector weightedAverageIFCumhazard_cpp(const arma::vec& seqTau, // horizon time for the predictions (T)
                                             const std::vector< arma::vec > & cumhazard0, // baseline cumulative hazard for each strata S:T
                                             const arma::mat& newX, // design matrix (Nxp)
                                             const arma::vec& neweXb, // exponential of the linear predictor (N)
                                             const arma::mat& IFbeta, //influence function for the regression coefficients (nxp)
                                             const std::vector< arma::mat >& cumEhazard0, // to compute the influence function of the baseline hazard S(pxJ)
                                             const std::vector< arma::vec >& cumhazard_iS0, // to compute the influence function of the baseline hazard S:J
                                             const arma::mat& delta_iS0, // to compute the influence function of the baseline hazard (nxS, S because multiplied by the strata indicator)
                                             const arma::mat& sample_eXb, // to compute the influence function of the baseline hazard (nxS, S because multiplied by the strata indicator)
                                             const arma::vec& sample_time, // event times of the training set (n)
                                             const std::vector< arma::uvec>& indexJumpSample_time, // index of the jump time corresponding to the sample times S:J
                                             const std::vector< arma::vec>& jump_time, // jump times for each strata S:J
                                             const std::vector< arma::uvec >& indexJumpTau, // index of the jump time corresponding to the horizon times S:J
                                             const arma::vec& lastSampleTime, // time after which the prediction are NA (S)
                                             const std::vector< arma::uvec>& newdata_index, // index of the observations within strata
                                             int nTau, int nSample, int nStrata, int p, 
                                             bool diag,
                                             int debug,
                                             const arma::vec& weights,
                                             bool isBeforeTau,
                                             double tau){ 
  
  // ** prepare
  if(debug>0){Rcpp::Rcout << "Prepare" << std::endl;}
  int iStrata_tauMin,iStrata_tauMax; 
  arma::vec iStrata_seqTau;
  arma::uvec iStrata_indexJumpTau;
  int iStrata_nNewObs; 
  int iJump; 
  int iTau;
  
  int iNewObs2; 
  arma::colvec iStrata_IFcumhazard0;
  arma::colvec iStrata_IFcumhazard;
  
  arma::uvec index_timestop(nSample);
  
  // ** initialize
  if(debug>0){Rcpp::Rcout << "Initialize" << std::endl;}
  
  arma::vec IF_cumhazard;
  IF_cumhazard.resize(nSample);
  IF_cumhazard.fill(0.0);
  
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
      // we nay have a problem here due to NAs in Brices code
      iStrata_tauMax--;
    }
    if(iStrata_tauMax < 0){continue;}
    
    R_CheckUserInterrupt();
    
    if(debug>1){Rcpp::Rcout << " (tau=" << iStrata_tauMin << "-" << iStrata_tauMax << ") ";}
    
    // *** compute IF/SE/IFmean at each time point
    for(int iTime=iStrata_tauMin; iTime<=iStrata_tauMax; iTime++){
      if(diag){
        iTau = newdata_index[iStrata](iTime);
      }else{
        iTau = iTime;
      }
      
      if (isBeforeTau){
        if (sample_time[iTau] > tau){
          break;
        }
      }
      
      // jump
      iJump = iStrata_indexJumpTau(iTime);
      
      if(debug>1){Rcpp::Rcout << " IF0 " ;}
      
      // **** IF baseline cumulative hazard hazard 
      index_timestop = indexJumpSample_time[iStrata];
      index_timestop.elem(find(index_timestop > iJump)).fill(iJump);
      
      iStrata_IFcumhazard0 = delta_iS0.col(iStrata) % (sample_time <= iStrata_seqTau(iTime)) - sample_eXb.col(iStrata) % cumhazard_iS0[iStrata](index_timestop);
      if(p>0){
        iStrata_IFcumhazard0 -= IFbeta * cumEhazard0[iStrata].col(iTau);
      }
      
      // **** IF/SE hazard/cumhazard/survival
      if(debug>1){Rcpp::Rcout << " IF " ;}
      for(int iNewObs=0; iNewObs<iStrata_nNewObs; iNewObs++){
        if(diag){
          iNewObs2 = newdata_index[iStrata](iTime);
        }else{
          iNewObs2 = newdata_index[iStrata](iNewObs); // the actual index 
        }
        
        if(p>0){
          iStrata_IFcumhazard = neweXb(iNewObs2)*(iStrata_IFcumhazard0 + cumhazard0[iStrata](iTau) * IFbeta * trans(newX.row(iNewObs2)));
        }else{
          iStrata_IFcumhazard = iStrata_IFcumhazard0;
        }
        //Rcout << iStrata_IFcumhazard << "\n";
        
        // store
        IF_cumhazard += iStrata_IFcumhazard*weights[iNewObs2];
      }
    } // end iTime
    if(debug>1){Rcpp::Rcout << std::endl;}
  } // end iStrata
  
  // ** Post process
  if(debug>0){Rcpp::Rcout << "Post process" << std::endl;}
  
  // ** Export
  if(debug>0){Rcpp::Rcout << "Export" << std::endl;}
  return(wrap(IF_cumhazard));
}
