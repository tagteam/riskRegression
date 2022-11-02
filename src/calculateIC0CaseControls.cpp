#include <Rcpp.h>
using namespace Rcpp;

// To replace the following:
// for (i in 1:n.controls){
//   ic0Case <- ic0Case + (1*(risk.cases > risk.controls[i])+0.5*1*(risk.cases == risk.controls[i]))*w.cases*w.controls[i]
// }
// for (i in 1:n.cases){
//   ic0Control <- ic0Control + (1*(risk.cases[i] > risk.controls)+0.5*1*(risk.cases[i] == risk.controls))*w.cases[i]*w.controls
// }
// [[Rcpp::export]]
void calculateIC0CaseControl(NumericVector ic0Case, NumericVector ic0Control, NumericVector riskCases, NumericVector riskControls, NumericVector wCases, NumericVector wControls) {
  int nCases = ic0Case.size();
  int nControls = ic0Control.size();
  for (int i = 0; i< nControls; i++){
    for (int j = 0; j < nCases; j++){
      if (riskCases[j] > riskControls[i]){
        ic0Case[j] += wCases[j]*wControls[i];
      }
      else if (riskCases[j] == riskControls[i]){
        ic0Case[j] += 0.5*wCases[j]*wControls[i];
      }
    }
  }
  
  for (int i = 0; i< nCases; i++){
    for (int j = 0; j < nControls; j++){
      if (riskCases[i] > riskControls[j]){
        ic0Control[j] += wCases[i]*wControls[j];
      }
      else if (riskCases[i] == riskControls[j]){
        ic0Control[j] += 0.5*wCases[i]*wControls[j];
      }
    }
  }
}
