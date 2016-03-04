#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

//// NOTE Adaptation of the Cagsurv5 from the survival package with no weight
// \item[ntimes] number of observations (unique death times)
// \item[nvar] number of covariates
// \item[ndead] number of deaths at that time
// \item[risk] weighted number at risk at the time
// \item[riskDead] sum of weights for the deaths

// [[Rcpp::export]]
NumericVector baseHazEfron_survival_cpp(int ntimes, int nvar, IntegerVector ndead, 
                               NumericVector risk, NumericVector riskDead) {
  double temp;
  int t, j;
  double di;
  NumericVector W(ntimes, 0.0);
  
  for (t=0; t < ntimes; t++) {
    di = ndead[t];
    
    if (di==1){
      
      W[t] = 1/risk[t];
      
    } else {
      
      temp = 1/risk[t];
      
      for (j=0; j < di; j++) {
        
        temp = 1/(risk[t] - riskDead[t]*j/di);
        W[t] += temp/di;
        
      }
      
    }
    
  }
  
  return(W);
}