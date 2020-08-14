// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat sampleMaxProcess_cpp(int nObject, int nNew, int nSim, const arma::cube& iid, const arma::mat& se);

// * quantileProcess_cpp
// [[Rcpp::export]]
NumericVector quantileProcess_cpp(int nObject, int nNew, int nSim,
                              arma::cube iid,
                              arma::mat se,
			      double confLevel){
  
 
  arma::mat maxTime_sample = sampleMaxProcess_cpp(nObject,
						  nNew,
						  nSim,
						  iid,
						  se);
  
  int indexQuantile = round(nSim * confLevel);
  colvec tempo;
  NumericVector Vquantile(nNew);
    
  for (int iCol = 0; iCol < nNew; iCol++){
    tempo = sort(maxTime_sample.col(iCol));
    Vquantile[iCol] = tempo[indexQuantile];
  }
  
  return(Vquantile);
}

// * sampleMaxProcess_cpp
// [[Rcpp::export]]
arma::mat sampleMaxProcess_cpp(int nObject, int nNew, int nSim,
                              const arma::cube& iid,
                              const arma::mat& se){

  void GetRNGstate(),PutRNGstate(); 
  GetRNGstate();

  colvec G;
  arma::mat iidG;
  arma::mat maxTime_sample(nSim,nNew);
  
  for(int iSim=0; iSim<nSim; iSim++){ 
  G = rnorm(nObject, 0, 1);
  iidG = iid.each_slice() * G;
  iidG = abs(iidG/se);
    
  for (int iCol = 0; iCol < nNew; ++iCol) { // each contrast take the largest statistic
    maxTime_sample(iSim, iCol) = iidG.col(iCol).max();
  }
  }
  
  PutRNGstate();
  
  return(maxTime_sample);
}

// * sampleMaxProcess2_cpp
// nObject: number of observations used to fit the model
// nNew: number of contrasts for which a different max is computed
// nSim: number of simulations
// value: observed value or null hypothesis (nTimes, nContrast)
// iid: influence function (nTimes, nObject, nContrast)
// se: standard error (nTimes, nContrast)
// type: 1 one sided below, 2 one sided above, 3 two sided
// [[Rcpp::export]]
arma::mat sampleMaxProcess2_cpp(int nObject, int nNew, int nSim,
								const arma::mat& value,
								const arma::cube& iid,
								const arma::mat& se,
								int type){

  void GetRNGstate(),PutRNGstate(); 
  GetRNGstate();

  colvec G;
  arma::mat iidG;
  arma::mat maxTime_sample(nSim,nNew);
  
  for(int iSim=0; iSim<nSim; iSim++){ 
	G = rnorm(nObject, 0, 1);
	iidG = iid.each_slice() * G;
	if(type == 1){ // one sided below
	  iidG = -(iidG - value)/se;
	}else if(type == 2){ // one sided above
	  iidG = (iidG - value)se;
	}else if(type == 3){ // two sided
	  iidG = (abs(iidG) - abs(value))/se;
	}
	
	for (int iCol = 0; iCol < nNew; ++iCol) { // each contrast take the largest statistic
	  maxTime_sample(iSim, iCol) = iidG.col(iCol).max();
	}
  }
  
  PutRNGstate();
  
  return(maxTime_sample);
}
