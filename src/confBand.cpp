// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat sampleMaxProcess_cpp(int nSample, int nContrast, int nSim,
							   const arma::colvec& value, const arma::cube& iid, int alternative, int type, bool global);

// * quantileProcess_cpp
// [[Rcpp::export]]
NumericVector quantileProcess_cpp(int nSample, int nContrast, int nSim,
								  arma::cube iid,
								  int alternative,
								  bool global,
								  double confLevel){

  arma::colvec colTempo(iid.n_slices);
  colTempo.fill(0.0);
  
  arma::mat maxTime_sample = sampleMaxProcess_cpp(nSample,
												  nContrast,
												  nSim,
												  colTempo,
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
							   const arma::colvec& value,
							   const arma::cube& iid,
							   int alternative,
							   int type,
							   bool global){

  void GetRNGstate(),PutRNGstate(); 
  GetRNGstate();
   
  colvec G;
  arma::mat iidG;
  arma::mat maxTime_sample(nSim,nContrast);
  bool rmValue = abs(value.max())>1e-12;
  if (rmValue == true){
	if(type==2) {         	// log() not defined here
	  throw std::range_error("When argument \'type\' is 2 then argument \'value\' should be 0");
	}
	if(type==3 && alternative != 2) {         	// log() not defined here
	  throw std::range_error("When argument \'type\' is 3 then argument \'alternative\' should be 2 or argument \'value\' should be 0");
	}
  }
  
  for(int iSim=0; iSim<nSim; iSim++){ 
	G = rnorm(nSample, 0, 1);

	if(alternative == 1){ // one sided below
	  iidG = iid.each_slice() * (-G);
	  if(rmValue){iidG.each_col() += value;}
	}else if(alternative == 2){ // one sided above
	  iidG = iid.each_slice() * G;
	  if(rmValue){iidG.each_col() -= value;}
	}else if(alternative == 3){ // two sided
	  iidG = abs(iid.each_slice() * G);
	  if(rmValue){iidG.each_col() -= abs(value);}
	}
	for (int iCol = 0; iCol < nContrast; ++iCol) { // each contrast take the largest statistic
	  if(type==1){
		maxTime_sample(iSim, iCol) = iidG.col(iCol).max();
	  }else if(type==2){
		maxTime_sample(iSim, iCol) = sum(iidG.col(iCol) % iidG.col(iCol));
	  }else if(type==3){
		maxTime_sample(iSim, iCol) = sum(iidG.col(iCol));
	  }
	}

	if(global){
	  maxTime_sample.row(iSim).fill(maxTime_sample.row(iSim).max()); // take the largest statistic over all contrasts
	}
  }
  
  PutRNGstate();
  
  return(maxTime_sample);
}
