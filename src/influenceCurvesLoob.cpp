
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix AUCijFun(NumericVector riskCase, NumericVector riskControl){

  int nCase = riskCase.size();
  int nControl = riskControl.size();
  DoubleVector Great(nControl), Equal(nControl);
  NumericMatrix out(nCase,nControl);
  
  for (int x = 0; x<nCase; x++){
    Great = riskCase[x] > riskControl;
    Equal = riskCase[x] == riskControl;
    out(x,_) = Great + Equal*0.5;
  }

  return out;
}

// [[Rcpp::export]]
DoubleVector icCensCC(NumericMatrix icCensC, NumericVector aucIJ){

  int Nc = icCensC.nrow(), N = icCensC.ncol();
  NumericVector aucij(Nc);
  DoubleVector out(N);

  for(int k = 0; k<N; k++){
    for (int i = 0; i<Nc; i++){
      aucij(i) = icCensC(i,k)*aucIJ(i);
    }
    out(k) = sum(aucij);
  }
  
  return out;
}

// [[Rcpp::export]]
DoubleVector icPhi(NumericMatrix icCensC, NumericVector weights){

  int Nc = icCensC.nrow(), N = icCensC.ncol();
  DoubleVector Weight(N), out(N);

    for (int k = 0; k<N; k++){
      for (int i = 0; i<Nc;i++){
	Weight[i] = icCensC(i,k)*weights[i];
      }
      out[k] = sum(Weight);
    }  
  return out;  
}

// [[Rcpp::export]]
NumericVector icWeightSubjectTimesFun(NumericMatrix icCensSubjectTimes, NumericVector icWeightSubjectTimes){

  int N = icCensSubjectTimes.nrow();
  DoubleVector ic1(icWeightSubjectTimes.size()); 
  NumericVector out(N);
  
  for (int k = 0; k<N; k++){
    ic1 = icCensSubjectTimes(k,_)*icWeightSubjectTimes;
    double ic = sum(ic1);
    out(k) = ic;
  }  
  return out;
}

// [[Rcpp::export]]
DoubleVector icWeightTimesFun(NumericVector icCensTimes, NumericVector icWeightTimes){

  int N = icCensTimes.size();
  DoubleVector ic1(icWeightTimes.size()); 
  DoubleVector out(N);
  
  for (int k = 0; k<N; k++){
    ic1 = icCensTimes[k]*icWeightTimes;
    double ic = sum(ic1);
    out(k) = ic;
  }  
  return out;
}


// [[Rcpp::export]]
DoubleVector icTauSubjectTimesFun(NumericMatrix icCensSubjectTimes, NumericVector whichCaseWeights, NumericVector weightsCase){

  int N = icCensSubjectTimes.nrow();
  double icMean, ic;
  DoubleVector icWeightSubjectTimes(whichCaseWeights.size());
  DoubleVector out(N);

    for (int k = 0; k<N; k++){
      icWeightSubjectTimes = icCensSubjectTimes(k,_)*pow(whichCaseWeights, 2);
      icMean = mean(icWeightSubjectTimes);
      ic = weightsCase[k] - icMean;
      out[k] = ic;      
  }
  return out;  
}



// [[Rcpp::export]]
DoubleVector icTauTimesFun(NumericVector icCensTimes, NumericVector whichControlWeights, NumericVector weightsControl){

  int N = icCensTimes.size();
  double icMean, ic;
  DoubleVector icWeightTimes(whichControlWeights.size());
  DoubleVector out(N);

    for (int k = 0; k<N; k++){
      icWeightTimes = icCensTimes[k]*pow(whichControlWeights, 2);
      icMean = mean(icWeightTimes);
      ic = weightsControl[k] - icMean;
      out[k] = ic;      
  }  
  return out;  
}


// [[Rcpp::export]]
DoubleVector icBrierWeightsSubjectTimesFun(NumericVector residuals, NumericMatrix icCensSubjectTimes, NumericVector weightsSubjectTimes){

  int N = icCensSubjectTimes.nrow();
  double ic;
  DoubleVector icWeightSubjectTimes(N);
  DoubleVector out(N);

    for (int k = 0; k<N; k++){
      icWeightSubjectTimes = icCensSubjectTimes(k,_)*weightsSubjectTimes*residuals;
      ic = mean(icWeightSubjectTimes);
      out[k] = ic;      
  }
  return out;  
}

// [[Rcpp::export]]
DoubleVector icBrierWeightsTimesFun(NumericVector residuals, NumericVector icCensTimes, NumericVector weightsTimes){

  int N = icCensTimes.size();
  double ic;
  DoubleVector icWeightTimes(N);
  DoubleVector out(N);

    for (int k = 0; k<N; k++){
      icWeightTimes = icCensTimes(k)*weightsTimes*residuals;
      ic = mean(icWeightTimes);
      out[k] = ic;      
  }
  return out;  
}
