// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


inline double calcIFhazard(double seqTau,
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
			   int p);
  
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
		      const arma::mat& newHazard, // predicted hazard
		      double firstJumpTime, double lastSampleTime,
		      int nTau, int nNewObs, int nSample, int p,
		      bool exportSE, bool exportIF){

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

  // narrow prediction times
  while((iTau0 < nTau) && firstJumpTime>seqTau[iTau0]){ // before the first event    
    iTau0++;    
  }
  while((iTauMax >= 0) && seqTau[iTauMax]>lastSampleTime){ // before the first event    
    iTauMax--;
  }

  // interesting times
  if(iTau0 < nTau && iTauMax >= 0){
    for(int iSample=0; iSample<nSample ; iSample++){
      R_CheckUserInterrupt();
      for(int iNewObs=0; iNewObs<nNewObs ; iNewObs++){

	for(int iTime=iTau0; iTime<=iTauMax; iTime++){
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
				      p);
	  // log transform
	  IF_cumhazard /= newHazard(iNewObs,iTime);

	  if(exportSE){
	    SEcumhazard(iNewObs,iTime) += pow(IF_cumhazard,2);
	  }
	  if(exportIF){
	    IFcumhazard(iNewObs,iTime,iSample) = IF_cumhazard;
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
    }
  }
 
  // SEcumhazard
  return(List::create(Named("iid") = IFcumhazard,
                      Named("se") = SEcumhazard));
}
// }}}

inline double calcIFhazard(double seqTau,
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
			   int p){

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
    if(sampleTime<=seqTau){
      IF_hazard0 = - IFbetaE - sampleEXb * hazard_iS0 + iS0;
    }else{
      IF_hazard0 = - IFbetaE - sampleEXb * hazard_iS0;
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
