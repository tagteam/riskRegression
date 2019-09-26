// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// * calcAugmentation
//' @title Compute Augmentation term
//'
//' @description Computes \int_index term(s) ds
//' @param term [matrix] values of the function to be integrated at jump times.
//' @param index [list of vector] the jump times at which to export the value of the integral, for each observation.
//' @param nObs [integer] the number of observations
//' @param nTau [integer] the number of evaluation times.
//'
//' @examples
//'
//'
//'
//'



// [[Rcpp::export]]
arma::mat calcAugmentation_cpp(const arma::mat& term,
							   const arma::mat& index,
							   int nObs, int nTau){

  arma::mat out(nObs,nTau);
  out.fill(NA_REAL);
  int iMaxJump;
  int iTau;
  double iInt;

  // ** compute integral
  for(int iObs=0; iObs<nObs; iObs++){
	iMaxJump = index(iObs,nTau-1);
	iTau = 0;
	iInt = 0.0;

	while((iTau<nTau) && (index(iObs,iTau)<0)){
	  out(iObs,iTau) = 0.0;
	  iTau++;
	}

	for(int iJump=0; iJump<=iMaxJump; iJump++){
	  iInt += term(iObs,iJump);

	  // while there are remaining times to export
	  // AND the prediction time is before the next event
	  //     OR its the last jump, i.e. the prediction time coincide with the last jump (otherwise it would have been removed at the begining)
	  while((iTau<nTau) && (index(iObs,iTau)==iJump)){
		out(iObs,iTau) = iInt;
		iTau++;
	  }
	}
  }

  return(out);  
}
