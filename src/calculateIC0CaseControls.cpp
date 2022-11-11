#include <Rcpp.h>
using namespace Rcpp;

// This should calculate
// ic0Case_i = sum_i (I(R_k > R_i) + 0.5 I(R_k == R_i)) Weights_i, where the sum is over controls and i any index
// ic0Controls_i = \sum_j W_k * sum_i (I(R_i > R_k) + 0.5 I(R_k == R_i)) Weights_i where the sum is taken over cases and i any index
// [[Rcpp::export]]
void calculateIC0CaseControl(NumericVector ic0Case, NumericVector ic0Control, NumericVector risk, LogicalVector isCase, LogicalVector isControl, NumericVector weight){
  int n = risk.size();
  IntegerVector ordering(n);
  std::iota(ordering.begin(), ordering.end(), 0);
  std::sort(ordering.begin(), ordering.end(),
            [&](int x, int y) { return risk[x] < risk[y]; });
  double valCurr{}, valPrev{};
  int i = n-1;
  while (i >= 0){
    int tieIter = i;
    while (tieIter >= 0 && risk[ordering[tieIter]]==risk[ordering[i]]){
      if (isCase[ordering[tieIter]]){
        valCurr += weight[ordering[tieIter]]; // should set something with valPrev up here
      }
      tieIter--;
    } 
    for (int l = i; l > tieIter;l--){
      if (isCase[ordering[l]]){
        ic0Control[ordering[l]] = 0.5*(valPrev+valCurr - weight[ordering[l]]);   // valPrev+0.5*(valCurr - weight[ordering[l]] - valPrev); //valPrev
      }
      else {
        ic0Control[ordering[l]] = 0.5*(valPrev+valCurr); // valPrev+0.5*(valCurr-valPrev); //valPrev
      }
    }
    i = tieIter;
    valPrev = valCurr;
  }
  
  valCurr = valPrev = 0;
  i = 0;
  while (i < n){
    int tieIter = i;
    while (tieIter < n && risk[ordering[tieIter]]==risk[ordering[i]]){
      if (isControl[ordering[tieIter]]){
        valCurr += weight[ordering[tieIter]];
      }
      tieIter++;
    } 
    for (int l = i; l < tieIter;l++){
      if (isControl[ordering[l]]){
        ic0Case[ordering[l]] = 0.5*(valPrev+valCurr - weight[ordering[l]]);
      }
      else {
        ic0Case[ordering[l]] = 0.5*(valPrev+valCurr);
      }
    }
    i = tieIter;
    valPrev = valCurr;
  }
}
