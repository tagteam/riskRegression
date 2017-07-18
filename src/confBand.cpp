// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::mat sampleMaxProcess_cpp(int nObject, int nNew, int nSim, const arma::cube& iid, const arma::mat& se);

// {{{ quantileProcess_cpp
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
// }}}

// {{{ sampleMaxProcess_cpp
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
    
  for (int iCol = 0; iCol < nNew; ++iCol) {
    maxTime_sample(iSim, iCol) = iidG.col(iCol).max();
  }
  }
  
  PutRNGstate();
  
  
  return(maxTime_sample);
}
// }}}
