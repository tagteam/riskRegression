// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat sampleMaxProcess_cpp(int nSample, int nContrast, int nSim,
							   const arma::mat& value, const arma::cube& iid, int alternative, int type, bool global);

// * quantileProcess_cpp
// [[Rcpp::export]]
NumericVector quantileProcess_cpp(int nSample, int nContrast, int nSim,
								  arma::cube iid,
								  int alternative,
								  bool global,
								  double confLevel){

  arma::mat matTempo = arma::zeros<arma::mat>(iid.n_slices, nContrast);
  
  arma::mat maxTime_sample = sampleMaxProcess_cpp(nSample,
												  nContrast,
												  nSim,
												  matTempo,
												  iid,
												  alternative,
												  1,
												  global);
  
  int indexQuantile = round(nSim * confLevel);
  colvec tempo;
  NumericVector Vquantile(nContrast);
    
  for (int iCol = 0; iCol < nContrast; iCol++){
    tempo = sort(maxTime_sample.col(iCol));
	if(alternative == 1){
	  Vquantile[iCol] = -tempo[indexQuantile];
	}else{
	  Vquantile[iCol] = tempo[indexQuantile];
	}
  }
  
  return(Vquantile);
}

// * sampleMaxProcess_cpp
// nSample: number of observations used to fit the model
// nContrast: number of contrasts for which a different max is computed
// nSim: number of simulations
// value: observed value or null hypothesis
// iid: influence function (nTimes, nSample, nContrast)
// se: standard error (nTimes, nContrast)
// alternative: 1 one sided below, 2 one sided above, 3 two sided
// type: 1 max test (Kolmogorov-Smirnov type supremum), 2 L2 test (Camer-von-Mises)
// [[Rcpp::export]]
arma::mat sampleMaxProcess_cpp(int nSample, int nContrast, int nSim,
							   const arma::mat& value,
							   const arma::cube& iid,
							   int alternative,
							   int type,
							   bool global){

  void GetRNGstate(),PutRNGstate(); 
  GetRNGstate();

  // ** check arguments
    bool rmValue = abs(value.max())>1e-12;
  if (rmValue == true){
	if(type == 2 && alternative != 3) {         	// log() not defined here
	  throw std::range_error("When argument \'type\' is 2 (CvM test) then argument \'alternative\' should be 3 (two-sided).");
	}
  }

  // ** prepare
  colvec G;
  arma::mat iidG;
  arma::mat maxTime_sample(nSim,nContrast);
  rowvec Svalue;
  if(type==1){
	if(alternative==1){
	  Svalue = min(value,0);
	}else if(alternative==2){
	  Svalue = max(value,0);
	}else if(alternative==3){
	  Svalue = max(abs(value),0);
	}
  }else if(type==2){
	Svalue = sum(value % value, 0);
  }else if(type==3){
	if(alternative==1){
	  Svalue = sum(value,0);
	}else if(alternative==2){
	  Svalue = sum(value,0);
	}else if(alternative==3){
	  Svalue = abs(sum(value,0));
	}
  }

  // ** run
  for(int iSim=0; iSim<nSim; iSim++){ 
	G = rnorm(nSample, 0, 1);
	iidG = iid.each_slice() * G;
 
	  for (int iCol = 0; iCol < nContrast; ++iCol) { // each contrast take the largest statistic
		if(type==1){
		  if(alternative==1){
			maxTime_sample(iSim, iCol) = iidG.col(iCol).min() - Svalue(iCol);
		  }else if(alternative==2){
			maxTime_sample(iSim, iCol) = iidG.col(iCol).max() - Svalue(iCol);
		  }else if(alternative==3){
			maxTime_sample(iSim, iCol) = abs(iidG.col(iCol)).max() - Svalue(iCol);
		  }
		}else if(type==2){
		  maxTime_sample(iSim, iCol) = sum(iidG.col(iCol) % iidG.col(iCol)) - Svalue(iCol);
		}else if(type==3){
		  if(alternative == 1){ // one sided below
			maxTime_sample(iSim, iCol) = - (sum(iidG.col(iCol)) - Svalue(iCol));
		  }else if(alternative == 2){ // one sided above
			maxTime_sample(iSim, iCol) = sum(iidG.col(iCol)) - Svalue(iCol);
		  }else if(alternative == 3){ // two sided
			maxTime_sample(iSim, iCol) = abs(sum(iidG.col(iCol))) - Svalue(iCol);
		  }
		}
	  }
	  
	  if(global){
		maxTime_sample.row(iSim).fill(maxTime_sample.row(iSim).max()); // take the largest statistic over all contrasts
	  }
  }
  
  PutRNGstate();
  
  return(maxTime_sample);
}
