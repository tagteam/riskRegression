#include <Rcpp.h>
using namespace Rcpp;

// This should calculate
// ic0Case_i = \sum_j I(R_i > R_j)*wControls_j*wCases_i where the sum is taken over controls
// ic0Controls_i = \sum_j I(R_j > R_i)*wCases_j*wControls_i where the sum is taken over cases
// [[Rcpp::export]]
void calculateIC0CaseControl(NumericVector ic0Case, NumericVector ic0Control, NumericVector riskCases, NumericVector riskControls, NumericVector wCases, NumericVector wControls) {
  int nCases = ic0Case.size();
  int nControls = ic0Control.size();
  
  for (int x = 0; x<nCases; x++){
    // LogicalVector log1 = riskCases[x] > riskControls;
    // LogicalVector log2 = riskCases[x] == riskControls;
    // NumericVector vec = (as<NumericVector>(log1)+0.5*as<NumericVector>(log2))*wControls;
    NumericVector vec(nControls);
    for (int y=0;y < nControls;y++){
      if (riskCases[x] > riskControls[y]){
        vec[y] = wControls[y];
      }
      else if (riskCases[x] == riskControls[y]){
        vec[y] = 0.5*wControls[y];
      }
    }
    ic0Case[x] = wCases[x]*sum(vec);
    ic0Control += wCases[x]*vec;
  }
}

// older version
// void calculateIC0CaseControl(NumericVector ic0Case, NumericVector ic0Control, NumericVector riskCases, NumericVector riskControls, NumericVector wCases, NumericVector wControls) {
//   int nCases = ic0Case.size();
//   int nControls = ic0Control.size();
//   for (int i = 0; i< nControls; i++){
//     for (int j = 0; j < nCases; j++){
//       if (riskCases[j] > riskControls[i]){
//         ic0Case[j] += wCases[j]*wControls[i];
//       }
//       else if (riskCases[j] == riskControls[i]){
//         ic0Case[j] += 0.5*wCases[j]*wControls[i];
//       }
//     }
//   }
//   
//   for (int i = 0; i< nCases; i++){
//     for (int j = 0; j < nControls; j++){
//       if (riskCases[i] > riskControls[j]){
//         ic0Control[j] += wCases[i]*wControls[j];
//       }
//       else if (riskCases[i] == riskControls[j]){
//         ic0Control[j] += 0.5*wCases[i]*wControls[j];
//       }
//     }
//   }
// }