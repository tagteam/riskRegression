// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

void sampleMaxProcess_cpp(int nSample, int nContrast, int nSim,
						  const arma::mat& value, arma::cube& iid, arma::mat& Msample, int alternative, int type, bool global);

// * quantileProcess_cpp
// Compute equicoordinate-quantile using simulations
// nSample: number of observations used to fit the model
// nContrast: number of contrasts for which a different max is computed
// nSim: number of simulations
// iid: influence function (nTimes, nSample, nContrast)
// alternative: 1 one sided below, 2 one sided above, 3 two sided
// global: [logical] should the max be taking over contrasts?
// [[Rcpp::export]]
NumericVector quantileProcess_cpp(int nSample, int nContrast, int nSim,
								  arma::cube& iid,
								  int alternative,
								  bool global,
								  double confLevel){

  void GetRNGstate(),PutRNGstate(); 
  GetRNGstate();

  // ** perform simulation
  arma::colvec G; // individual weights
  arma::mat iidG; // temporary curve
  arma::mat Mstore(nSim, nContrast); // store simulation results

  for(int iSim=0; iSim<nSim; iSim++){ 
	G = rnorm(nSample, 0, 1);
	iidG = iid.each_slice() * G;
 
	for (int iC = 0; iC < nContrast; iC++) { // take the more extreme statistic (over time) for each contrast 
	  if(alternative==1){
		Mstore(iSim,iC) = iidG.col(iC).min();
	  }else if(alternative==2){
		Mstore(iSim,iC) = iidG.col(iC).max();
	  }else if(alternative==3){
		Mstore(iSim,iC) = abs(iidG.col(iC)).max();
	  }
	}
	  
	if(global){ // take the more extreme statistic (over contrasts)
	  if(alternative==1){ 
		Mstore.row(iSim).fill(Mstore.row(iSim).min()); 
	  }else{
		Mstore.row(iSim).fill(Mstore.row(iSim).max()); 
	  }
	}
  }

  // ** compute quantile
  int indexQuantile;
  if(alternative==1){
	indexQuantile = round(nSim * (1-confLevel));
  }else{
	indexQuantile = round(nSim * confLevel);
  }
  arma::colvec tempo;
  NumericVector Vquantile(nContrast);
    
  for (int iCol = 0; iCol < nContrast; iCol++){
    tempo = sort(Mstore.col(iCol));
	Vquantile[iCol] = tempo[indexQuantile];
  }
  
  PutRNGstate();
  return(Vquantile);
}

// * pProcess_cpp
// Compute p-values adjusted for multiple comparisons using simulations
// nSample: number of observations used to fit the model
// nContrast: number of contrasts for which a different max is computed
// nTime: number of timepoints
// nSim: number of simulations
// value: observed value or null hypothesis
// iid: influence function (nTimes, nSample, nContrast)
// alternative: 1 one sided below, 2 one sided above, 3 two sided
// global: [logical] should the max be taking over contrasts?
// [[Rcpp::export]]
arma::mat pProcess_cpp(int nSample, int nContrast, int nTime, int nSim,
						   arma::mat value,
						   arma::cube& iid,
						   int alternative,
						   bool global){

  void GetRNGstate(),PutRNGstate(); 
  GetRNGstate();

  // ** prepare
  arma::mat pmat(nContrast,nTime);
  pmat.fill(0.0);

  arma::colvec G(nSample); // individual weights
  arma::mat iidG(nTime,nContrast); // temporary curve
  if(alternative==3){
	value = abs(value);
  }
  double iEx=NA_REAL;
  
  // ** simulation
  for(int iSim=0; iSim<nSim; iSim++){ 
	G = rnorm(nSample, 0, 1);
	iidG = iid.each_slice() * G;

	if(global==true){
	  if(alternative==1){
		iEx = iidG.min();
	  }else if(alternative==2){
		iEx = iidG.max();
	  }else if(alternative==3){
		iEx = abs(iidG).max();
	  }
	}
	for (int iC = 0; iC < nContrast; iC++) { // take the more extreme statistic (over time) for each contrast
	  if(global==false){
		if(alternative==1){
		  iEx = iidG.col(iC).min();
		}else if(alternative==2){
		  iEx = iidG.col(iC).max();
		}else if(alternative==3){
		  iEx = abs(iidG.col(iC)).max();
		}
	  }
	  for (int iT = 0; iT < nTime; iT++) { // for each time
		if(alternative == 1){
		  if(iEx <= value(iC,iT)){
			pmat(iC,iT)++ ;
		  }
		}else{
		  if(iEx >= value(iC,iT)){
			pmat(iC,iT)++ ;
		  }
		}
	  }
	}
	
  }

  PutRNGstate();
  
  return(pmat / nSim);
}

// * sampleMaxProcess_cpp
// nSample: number of observations used to fit the model
// nContrast: number of contrasts for which a different max is computed
// nSim: number of simulations
// value: observed value or null hypothesis
// iid: influence function (nTimes, nSample, nContrast)
// alternative: 1 one sided below, 2 one sided above, 3 two sided
// type: 1 max test (Kolmogorov-Smirnov type supremum), 2 L2 test (Camer-von-Mises)
// global: [logical] should the max be taking over contrasts?
// [[Rcpp::export]]
arma::mat sampleMaxProcess_cpp(int nSample, int nContrast, int nSim,
							   const arma::mat& value,
							   arma::cube& iid,
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
  arma::colvec G;
  arma::mat iidG;
  arma::mat maxTime_sample(nSim,nContrast);
  arma::rowvec Svalue;
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
