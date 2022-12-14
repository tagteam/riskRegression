// [[Rcpp::depends(RcppArmadillo)]]
#include "IC-Nelson-Aalen-cens-time.h"

using namespace Rcpp;
using namespace arma;

// part of IFAUC without the influence function from the censoring
// author: Johan Sebastian Ohlendorff
// [[Rcpp::export(rng = false)]]
List getIC0AUC(NumericVector time,
                        NumericVector status,
                        double tau,
                        NumericVector risk,
                        NumericVector GTiminus,
                        NumericVector Gtau,
                        double auc) {
  int n = time.size();
  NumericVector ic0(n), ic0Case(n), ic0Control(n), weights(n);
  LogicalVector cases(n), controls1(n), controls2(n);
  double muCase{}, muControls{}, nu{};
  
  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  auto lower = std::upper_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower) -1;
  if (firsthit == -1){
    firsthit = 0;
  }

  // calculate weights W_t(G;Z_i) = I(status_i != 0, time <= tau) 1/G(Ti-|Xi) + I(time > tau) 1/G(Ti-|Xi) 
  // also calculate muCase = sum_i W_t(G;Z_i) over cases and muControls = sum_i W_t(G;Z_i) over controls
  for (int i = 0; i <= firsthit; i++){
    if(status[i] == 1){ // case
      cases[i] = true;
      weights[i] = 1.0/GTiminus[i];
      muCase += weights[i];
    }
    else if (status[i] == 2){ // control
      controls2[i] = true;
      weights[i] = 1.0/GTiminus[i];
      muControls += weights[i];
    }
  }
  for (int i = firsthit + 1; i < n; i++){ // control
    controls1[i] = true;
    weights[i] = 1.0/Gtau[i];
    muControls += weights[i];
  }
  double mu = muCase*muControls / (double (n*n));
  nu = auc * mu;

  IntegerVector ordering(n);
  std::iota(ordering.begin(), ordering.end(), 0);
  std::sort(ordering.begin(), ordering.end(),
            [&](int x, int y) { return risk[x] < risk[y]; });
  double valCurr{}, valPrev{};
  int i = n-1;
  while (i >= 0){
    int tieIter = i;
    while (tieIter >= 0 && risk[ordering[tieIter]]==risk[ordering[i]]){
      if (cases[ordering[tieIter]]){
        valCurr += weights[ordering[tieIter]]; // should set something with valPrev up here
      }
      tieIter--;
    } 
    for (int l = i; l > tieIter;l--){
      if (cases[ordering[l]]){
        ic0Control[ordering[l]] = weights[ordering[l]]*0.5*(valPrev+valCurr - weights[ordering[l]]);   // valPrev+0.5*(valCurr - weight[ordering[l]] - valPrev); //valPrev
      }
      else {
        ic0Control[ordering[l]] = weights[ordering[l]]*0.5*(valPrev+valCurr); // valPrev+0.5*(valCurr-valPrev); //valPrev
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
      if (controls1[ordering[tieIter]] || controls2[ordering[tieIter]]){
        valCurr += weights[ordering[tieIter]];
      }
      tieIter++;
    } 
    for (int l = i; l < tieIter;l++){
      if (controls1[ordering[l]] || controls2[ordering[l]]){
        ic0Case[ordering[l]] = weights[ordering[l]]*0.5*(valPrev+valCurr - weights[ordering[l]]);
      }
      else {
        ic0Case[ordering[l]] = weights[ordering[l]]*0.5*(valPrev+valCurr);
      }
    }
    i = tieIter;
    valPrev = valCurr;
  }
  // Rcout << "ic0Case: " << ic0Case << "\n";
  // Rcout << "ic0Control: " << ic0Control << "\n";
  
  double IF0num{}, IF0den{};
  for (int i = 0; i < n; i++){
    if (time[i] <= tau && status[i] == 1){ // case
      IF0num = ic0Case[i] / (double (n));
      IF0den = weights[i] / (double (n)) * muControls;
    }
    else if ((time[i] <= tau && status[i] == 2) || (time[i] > tau)){ // control
      IF0num = ic0Control[i] / (double (n));
      IF0den = weights[i] / (double (n)) * muCase;
    }
    else {
      IF0num = IF0den = 0;
    }
    ic0[i] = (IF0num * mu - IF0den * nu)/(mu*mu);
  }
  return(List::create(Named("ic0") = ic0,
                      Named("ic0Case") = ic0Case[cases],
                      Named("ic0Control") = ic0Control[controls1 | controls2],
                      Named("weights") = weights,
                      Named("muCase") = muCase,
                      Named("muControls") = muControls,
                      Named("nu") = nu,
                      Named("firsthit")=firsthit,
                      Named("cases")=cases,
                      Named("controls")=controls1 | controls2,
                      Named("controls1")=controls1,
                      Named("controls2")=controls2));
}

// calculate the term corresponding to KM censoring
// author: Johan Sebastian Ohlendorff
// [[Rcpp::export(rng = false)]]
NumericVector getInfluenceFunctionAUCKMCensoringTerm(NumericVector time,
                                                     NumericVector status,
                                                     double tau,
                                                     NumericVector ic0Case,
                                                     NumericVector ic0Controls,
                                                     NumericVector weights,
                                                     int firsthit,
                                                     double muCase,
                                                     double muControls,
                                                     double nu1,
                                                     double Gtau,
                                                     double auc, 
                                                     int startControls1) {
  // Thomas' code from IC of Nelson-Aalen estimator, i.e. calculate the influence function of the hazard
  // initialize first time point t=0 with data of subject i=0
  int n = time.size();
  NumericVector icpart(n);
  arma::uvec sindex(n,fill::zeros);
  arma::vec utime=unique(time);
  int nu=utime.size();
  arma::vec atrisk(nu);
  arma::vec MC_term2(nu,fill::zeros);
  getInfluenceFunctionKM(time,status,atrisk,MC_term2,sindex,utime);
  double term2numpart1{}, term2denpart1{};
  double term2numpart2 = nu1 * (double (n*n)); //{}, term3num{}, term2denpart2{}, term3den{};
  double term2denpart2 = muCase;
  double term3numpart1{}, term3denpart1{};
  double term3numpart2{};
  double sumPartNum{};
  int nControls = ic0Controls.size();
  for (int i = startControls1; i < nControls; i++){ // control with time > t
    sumPartNum += ic0Controls[i];
  }
  term3numpart2 = nu1 * (double (n*n)) - sumPartNum; //{}, term3num{}, term2denpart2{}, term3den{};
  double sumPartDen = double ((n-(firsthit+1)))/Gtau;
  double term3denpart2 = muControls - sumPartDen;
  double mu = muCase*muControls / (n*n);
  int tieIter{}, tieIterCases{}, tieIterControls{};
  // can do while loops together
  while ((tieIter < n) && (time[tieIter] == time[0])) {
    if ((time[tieIter] <= tau) && (status[tieIter]==1)){
      term2numpart2 -= ic0Case[tieIterCases]; 
      term2denpart2 -= weights[tieIter]; 
      tieIterCases++;
    }
    else if  ((time[tieIter] <= tau) && (status[tieIter]==2)){
      term3numpart2 -= ic0Controls[tieIterControls]; 
      term3denpart2 -= weights[tieIter]; 
      tieIterControls++;
    }
    else if (time[tieIter] > tau){
      tieIterControls++;
    }
    tieIter++;
  }
  int upperTie = tieIter-1;
  double term2num{}, term3num{}, term2den{}, term3den{}, fihattau{};
  for (int i = 0; i<n; i++){
    int const j = i > firsthit ? firsthit : i;
    if (utime[sindex[j]] < time[i]){
      fihattau = - MC_term2[sindex[j]];
    }
    else {
      fihattau =  (1-(status[i] != 0))*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
    }
    term2num = 1.0 / (n*n) * (term2numpart1+fihattau*term2numpart2); // is wrong?
    term2den = 1.0 / (n*n) * (term2denpart1+fihattau*term2denpart2)*muControls; // ok
    term3num = 1.0 / (n*n) * (term3numpart1+fihattau*term3numpart2+fihattau*sumPartNum); // is ok 
    term3den = 1.0 / (n*n) * (term3denpart1+fihattau*term3denpart2+fihattau*sumPartDen)*muCase; // is ok 
    // Rcout << term2num << "\n";
    // Rcout << term2den<< "\n\n";
    // Rcout << term3num<< "\n";
    // Rcout << term3den<< "\n";
    // fast calculation of eq10 and eq17
    if (upperTie == i){
      int tieIter = i+1;
      while ((tieIter < n) && (time[tieIter] == time[i+1])) {
        if ((time[tieIter] <= tau) && (status[tieIter]==1)){
          term2numpart1 -= ic0Case[tieIterCases] * (MC_term2[sindex[i]]); 
          term2denpart1 -= (MC_term2[sindex[i]]) * weights[tieIter]; 
          term2numpart2 -= ic0Case[tieIterCases]; 
          term2denpart2 -= weights[tieIter]; 
          tieIterCases++;
        }
        else if ((time[tieIter] <= tau) && (status[tieIter]==2)){
          term3numpart1 -= ic0Controls[tieIterControls]* MC_term2[sindex[i]]; 
          term3denpart1 -= weights[tieIter] * MC_term2[sindex[i]]; 
          term3numpart2 -= ic0Controls[tieIterControls]; 
          term3denpart2 -= weights[tieIter]; 
          tieIterControls++;
        }
        else if (time[tieIter] > tau){
          tieIterControls++;
        }
        tieIter++;
      }
      upperTie = tieIter-1;
    }
    double IFnum = term2num+term3num;
    double IFden = term2den+term3den;
    
    icpart[i] = (IFnum * mu - IFden * nu1)/(mu*mu);
  }
  return icpart;
}