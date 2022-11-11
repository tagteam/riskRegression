// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"
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
      weights[i] = 1.0/GTiminus[i];
      muCase += weights[i];
    }
    else if (status[i] == 2){ // control
      weights[i] = 1.0/GTiminus[i];
      muControls += weights[i];
    }
  }
  for (int i = firsthit + 1; i < n; i++){ // control
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
      if (time[ordering[tieIter]] <= tau && status[ordering[tieIter]] == 1){
        valCurr += weights[ordering[tieIter]]; // should set something with valPrev up here
      }
      tieIter--;
    } 
    for (int l = i; l > tieIter;l--){
      if (time[ordering[l]] <= tau && status[ordering[l]] == 1){
        ic0Control[ordering[l]] = 0.5*(valPrev+valCurr - weights[ordering[l]]);   // valPrev+0.5*(valCurr - weight[ordering[l]] - valPrev); //valPrev
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
      if ((time[ordering[tieIter]] <= tau && status[ordering[tieIter]] == 2) || time[ordering[tieIter]] > tau){
        valCurr += weights[ordering[tieIter]];
      }
      tieIter++;
    } 
    for (int l = i; l < tieIter;l++){
      if ((time[ordering[l]] <= tau && status[ordering[l]] == 2) || time[ordering[l]] > tau){
        ic0Case[ordering[l]] = 0.5*(valPrev+valCurr - weights[ordering[l]]);
      }
      else {
        ic0Case[ordering[l]] = 0.5*(valPrev+valCurr);
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
      IF0num = weights[i] / (double (n)) * ic0Case[i];
      IF0den = weights[i] / (double (n)) * muControls;
    }
    else if ((time[i] <= tau && status[i] == 2) || (time[i] > tau)){ // control
      IF0num = weights[i] / (double (n)) * ic0Control[i];
      IF0den = weights[i] / (double (n)) * muCase;
    }
    else {
      IF0num = IF0den = 0;
    }
    ic0[i] = (IF0num * mu - IF0den * nu)/(mu*mu);
  }
  return(List::create(Named("ic0") = ic0,
                      Named("ic0Case") = ic0Case,
                      Named("ic0Control") = ic0Control,
                      Named("weights") = weights,
                      Named("muCase") = muCase,
                      Named("muControls") = muControls,
                      Named("nu") = nu,
                      Named("firsthit")=firsthit));
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
                                                     double auc) {
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
  double sumPartNum{};
  for (int i = firsthit + 1; i < n; i++){ // control with time > t
    sumPartNum += ic0Controls[i];
  }
  sumPartNum /= Gtau;
  double sumPartDen = double ((n-(firsthit+1)))/Gtau;
  double term3numpart2 = nu1 * (double (n*n)) - sumPartNum; //{}, term3num{}, term2denpart2{}, term3den{};
  double term3denpart2 = muControls - sumPartDen;
  double mu = muCase*muControls / (n*n);
  int tieIter = 0;
  // can do while loops together
  while ((tieIter < n) && (time[tieIter] == time[0])) {
    if ((time[tieIter] <= tau) && (status[tieIter]==1)){
      term2numpart2 -= ic0Case[tieIter] * weights[tieIter]; 
      term2denpart2 -= weights[tieIter]; 
    }
    else if  ((time[tieIter] <= tau) && (status[tieIter]==2)){
      term3numpart2 -= ic0Controls[tieIter] * weights[tieIter]; 
      term3denpart2 -= weights[tieIter]; 
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
          term2numpart1 -= ic0Case[tieIter] * (MC_term2[sindex[i]]) * weights[tieIter]; 
          term2denpart1 -= (MC_term2[sindex[i]]) * weights[tieIter]; 
          term2numpart2 -= ic0Case[tieIter] * weights[tieIter]; 
          term2denpart2 -= weights[tieIter]; 
        }
        else if ((time[tieIter] <= tau) && (status[tieIter]==2)){
          term3numpart1 -= ic0Controls[tieIter]* MC_term2[sindex[i]]*weights[tieIter]; 
          term3denpart1 -= weights[tieIter] * MC_term2[sindex[i]]; 
          term3numpart2 -= ic0Controls[tieIter]*weights[tieIter]; 
          term3denpart2 -= weights[tieIter]; 
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





// see https://github.com/eestet75/riskRegressionStudy/blob/master/PicsForImplementation/crossvalAUC.png
// [[Rcpp::export(rng = false)]]
NumericVector getInfluenceFunctionAUCKMCensoringCVPart(NumericVector time,
                                                       NumericVector status,
                                                       double tau,
                                                       NumericVector GTiminus,
                                                       double Gtau,
                                                       NumericMatrix aucMat,
                                                       double nu1tauPm) {
  // check for NAs and equal lengths of vectors
  checkNAs(time, GET_VARIABLE_NAME(time));
  checkNAs(status, GET_VARIABLE_NAME(status));
  checkNAs(tau, GET_VARIABLE_NAME(tau));
  checkNAs(GTiminus, GET_VARIABLE_NAME(GTiminus));
  checkNAs(Gtau, GET_VARIABLE_NAME(Gtau));
  
  // should also check matrix for NAs and the IntegerVectors
  compareLengths(time,status);
  compareLengths(status,GTiminus);
  
  int n = time.size();
  NumericVector ic(n);
  arma::uvec sindex(n,fill::zeros);
  arma::vec utime=unique(time);
  int nu=utime.size();
  arma::vec atrisk(nu);
  arma::vec MC_term2(nu,fill::zeros);
  getInfluenceFunctionKM(time,status,atrisk,MC_term2,sindex,utime);
  
  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  auto lower = std::upper_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower) -1;
  if (firsthit == -1){
    firsthit = 0;
  }
  
  // P(tilde{T_i} > tau)
  double Probmu = double ((n-(firsthit+1))) / n;
  //this is the one theta_m(X_i,x), varies over cases
  NumericVector int1 = rowSums(aucMat) / ((double) n); // this is \frac{1}{G(\tilde{T_i}- | Z_i)}\int  \Theta_m(X_i,x') \left(1_{\{t^{\prime}>\tau\}} \frac{1}{G(\tau | z')} + I(delta = 2)/ G(tilde(T_i)-) \right) P(dx')ordered according to cases
  // this is the one with theta_m(x,X_i), varies over controls
  NumericVector int2 = colSums(aucMat) / ((double) n); //  \int  \Theta_m(x,X_i)  1_{\{t \leqslant \tau\}}  \frac{P(dx) \I{\delta=1}}{G(t- | z)}\left( \frac{1_{\{T_i>\tau\}}}{G(\tau | Z_i)} + 1_{\{T_i leqslant \tau\}} \frac{\I{\Delta_i=2}}{G(T_i - | Z_i)} \right) ordered according to controls
  // int1mu = int (1_{t \leq tau} 1(delta = 2)/ hat{G(t-)} + 1_{t > tau} / G(tau)) dP(z) = Q(T_i <= tau, Delta_i = 2 | T_i > tau), also int2mu
  // int2mu = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,1) = Q(T_i <= tau, Delta_i = 1)
  double int1mu = 0, int2mu = 0, eq32part2 = 0, eq33part2 = 0, eq33term1part = 0;
  int numberOfCases = 0;
  int numberOfControls = 0;
  for (int i = 0; i < n; i++){
    if (time[i] <= tau && status[i] == 1){
      int2mu += 1.0/GTiminus[i];
      numberOfCases++;
    }
    else if (time[i] <= tau && status[i] == 2){
      int1mu += 1.0/GTiminus[i];
      eq33part2+=int2[numberOfControls];
      numberOfControls++;
    }
    else if (time[i] > tau){
      eq33term1part+=int2[numberOfControls];
      numberOfControls++;
    }
  }
  int1mu = int1mu / ((double) n) + Probmu / Gtau;
  int2mu = int2mu / ((double) n);
  eq33term1part = eq33term1part / ((double) n);
  eq32part2 = n*n*nu1tauPm;
  eq33part2 *= n;
  // Rcout << "eq32part2 should be " << eq32part2 << "\n";
  // Rcout << "eq33term1part is " << eq33term1part << "\n";
  double mu1 = int1mu*int2mu;
  double nu1 = nu1tauPm;
  // Rcout << "maybe" << n*n*nu1;
  // Rcout << "eq33part2 is now " << eq33part2 << "\n";
  
  double eq32part1{}, eq33part1{}, eq36part1{}, eq37part1{}, eq33term1{}, eq33term2{}, eq37term1{}, eq37term2{};
  double eq36part2 = n*int2mu; //36 was previously 17 // note need to initialize values?
  // Rcout << eq36part2;
  double eq37part2 = n*(int1mu-Probmu / Gtau);
  int tieIter = 0;
  int tieIterCases = 0;
  int tieIterControls = 0;
  // can do while loops together
  while ((tieIter < n) && (time[tieIter] == time[0])) {
    if ((time[tieIter] <= tau) && (status[tieIter]==1)){
      // Rcout << "tieIterCases " << tieIterCases << "\n";
      eq32part2 -= n*int1[tieIterCases];
      eq36part2 -= 1.0 / GTiminus[tieIter];
      tieIterCases++;
    }
    else if ((time[tieIter] <= tau) && (status[tieIter]==2)){ 
      eq33part2 -= n*int2[tieIterControls];
      eq37part2 -= 1.0 / GTiminus[tieIter];
      tieIterControls ++;
    }
    else if (time[tieIter] > tau){
      tieIterControls ++;
    }
    tieIter++;
  }
  
  int upperTie = tieIter-1;
  numberOfCases = 0;
  numberOfControls = 0;
  double eq31{}, eq32{}, eq33{}, eq34{}, eq35{}, eq36{}, eq37{}, eq38{}, fihattau{};
  for (int i = 0; i<n; i++){
    int const j = i > firsthit ? firsthit : i;
    if (utime[sindex[j]] < time[i]){
      fihattau = - MC_term2[sindex[j]];
    }
    else {
      fihattau =  (1-(status[i] != 0))*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
    }
    //Rcout << "fihat is " << fihattau << "\n";
    if (time[i] <= tau && status[i] == 1){
      eq31 = int1[numberOfCases];
      eq34 = 0;
      eq35 = 1.0 / GTiminus[i] * int1mu;
      eq38 = 0;
      numberOfCases++;
    }
    else if ((time[i] <= tau && status[i] == 2) || time[i] > tau){
      eq31 = 0;
      eq34 = int2[numberOfControls];
      eq35 = 0;
      eq38 = (time[i] > tau) ? 1.0 / Gtau * int2mu : 1.0/GTiminus[i] * int2mu;
      numberOfControls++;
    }
    else {
      eq31=eq34=eq35=eq38=0;
    }
    // Rcout << "eq32part1 " << n*n*eq32part1 << "\n";
    // Rcout << "eq32part2 " << n*n*eq32part2 << "\n";
    eq32 = 1.0 / (n*n) * (eq32part1+fihattau*eq32part2);
    eq33term1 = fihattau * eq33term1part;
    // Rcout << "eq33term1 is " << eq33term1 << "\n";
    eq33term2 = 1.0 / (n*n) * (eq33part1+fihattau*eq33part2);
    // Rcout << "eq33term2 is " << eq33term2 << "\n";
    // Rcout << "eq33part1 " << eq33part1 << "\n";
    // Rcout << "eq33part2 " << eq33part2 << "\n";
    eq33 = eq33term1 + eq33term2;
    
    // Rcout << "eq36part1 is " << eq36part1 << "\n";
    // Rcout << "eq36part2 is " << eq36part2 << "\n";
    eq36 = int1mu * 1.0 / (n) * (eq36part1+fihattau*eq36part2);
    eq37term1 = fihattau / Gtau * Probmu * int2mu;
    eq37term2 = 1.0 / (n) * (eq37part1+fihattau*eq37part2) * int2mu;
    eq37 = eq37term1 + eq37term2;
    
    if (upperTie == i){
      int tieIter = i+1;
      while ((tieIter < n) && (time[tieIter] == time[i+1])) {
        if ((time[tieIter] <= tau) && (status[tieIter]==1)){
          // Rcout << "tieIterCases " << tieIterCases << "\n";
          eq32part1 -= n*int1[tieIterCases] * (MC_term2[sindex[i]]);
          eq32part2 -= n*int1[tieIterCases];
          eq36part1 -= (MC_term2[sindex[i]]) / GTiminus[tieIter];
          eq36part2 -= 1.0 / GTiminus[tieIter];
          tieIterCases++;
        }
        else if ((time[tieIter] <= tau) && (status[tieIter]==2)){
          eq33part1 -= n*int2[tieIterControls] * (MC_term2[sindex[i]]);
          eq33part2 -= n*int2[tieIterControls];
          eq37part1 -= (MC_term2[sindex[i]]) / GTiminus[tieIter];
          eq37part2 -= 1.0 / GTiminus[tieIter];
          tieIterControls++;
        }
        else if (time[tieIter] > tau){
          tieIterControls ++;
        }
        tieIter++;
      }
      upperTie = tieIter-1;
    }
    // Rcout << "iteration i = " << i << "\n";
    // // Rcout << "eq31 should be: " << eq31 << "\n";
    // // Rcout << "eq32 should be: " << eq32 << "\n";
    // Rcout << "eq33 should be: " << eq33 << "\n";
    // Rcout << "eq34 should be: " << eq34 << "\n";
    // Rcout << "eq35 should be: " << eq35 << "\n";
    // Rcout << "eq36 should be: " << eq36 << "\n";
    // Rcout << "eq37 should be: " << eq37 << "\n";
    // Rcout << "eq38 should be: " << eq38 << "\n";
    
    double IFnu = eq31+eq32+eq33+eq34;
    double IFmu = eq35+eq36+eq37+eq38;
    
    ic[i] = (IFnu * mu1 - IFmu * nu1)/(mu1*mu1);
  }
  return ic;
}