// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"
#include "IC-Nelson-Aalen-cens-time.h"

using namespace Rcpp;
using namespace arma;

// Calculate influence function for competing risk case/survival case with Nelson-Aalen censoring.
// Equation lines correspond to the document influenceFunctionAUC.pdf in riskRegressionStudy, see
// https://github.com/eestet75/riskRegressionStudy/blob/master/PicsForImplementation/AUCtraintestKM.png
// Author: Johan Sebastian Ohlendorff
// [[Rcpp::export(rng = false)]]
NumericVector getInfluenceFunctionAUCKMCensoring(NumericVector time,
                                                 NumericVector status,
                                                 double tau,
                                                 NumericVector risk,
                                                 NumericVector GTiminus,
                                                 double Gtau,
                                                 double auc,
                                                 bool tiedValues) {
  // check if any of the vectors have NAs and also that the vectors have the same lengths
  checkNAs(time, GET_VARIABLE_NAME(time));
  checkNAs(status, GET_VARIABLE_NAME(status));
  checkNAs(tau, GET_VARIABLE_NAME(tau));
  checkNAs(risk,GET_VARIABLE_NAME(risk));
  checkNAs(GTiminus,GET_VARIABLE_NAME(GTiminus));
  checkNAs(Gtau,GET_VARIABLE_NAME(Gtau));
  checkNAs(auc, GET_VARIABLE_NAME(auc));
  compareLengths(time,status);
  compareLengths(status,risk); 
  compareLengths(risk,GTiminus);
  
  // Thomas' code from IC of Nelson-Aalen estimator, i.e. calculate the influence function of the hazard
  // initialize first time point t=0 with data of subject i=0
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
  // If ties = TRUE, there are equalities for X in the below
  //P(X < X_i, tilde{T_i} > tau)
  NumericVector Probnu(n);
  // int 1*(X_i < x, t <= tau) / G(t-) dP(t,1,x) 
  NumericVector eq1112part(n);
  // int 1*(X_i > x, t <= tau) / G(t-) dP(t,2,x) 
  NumericVector eq1314part(n);
  // risk ordering, i.e. order according to risk 
  IntegerVector ordering(n);
  std::iota(ordering.begin(), ordering.end(), 0);
  std::sort(ordering.begin(), ordering.end(),
            [&](int x, int y) { return risk[x] < risk[y]; });
  // efficient computation of eq1112part and eq1314part by using the ordering of risk
  // the integrals differ if we are trying to compute the part concerning ties in risk
  if (!tiedValues){
    int jCurr = 0, jPrev = 0, i = 0;
    double valCurr = 0, valPrev = 0;
    while (i < n){
      int tieIter = i;
      while (tieIter < n && risk[ordering[tieIter]]==risk[ordering[i]]){
        if (time[ordering[tieIter]] > tau){
          jCurr++;
        }
        else if (time[ordering[tieIter]] <= tau && status[ordering[tieIter]] == 2){
          valCurr += 1.0/(GTiminus[ordering[tieIter]]);
        }
        tieIter++;
      }
      for (int l = i; l < tieIter;l++){
        Probnu[ordering[l]] = ((double) jPrev)/n;
        eq1314part[ordering[l]] = valPrev;
      }
      i = tieIter;
      jPrev = jCurr;
      valPrev = valCurr;
    }
    valCurr = valPrev = 0;
    i = n-1;
    while (i >= 0){
      int tieIter = i;
      while (tieIter >= 0 && risk[ordering[tieIter]]==risk[ordering[i]]){
        if (time[ordering[tieIter]] <= tau && status[ordering[tieIter]] == 1){
          valCurr += 1.0/(GTiminus[ordering[tieIter]]);
        }
        tieIter--;
      } 
      for (int l = i; l > tieIter;l--){
        eq1112part[ordering[l]] = valPrev;
      }
      i = tieIter;
      valPrev = valCurr;
    }
  }
  else {
    int i = 0;
    while (i < n){
      int jProbnu = 0;
      double jPart1112 = 0, jPart1314 = 0;
      int tieIter = i;
      while (tieIter < n && risk[ordering[tieIter]]==risk[ordering[i]]){
        if (time[ordering[tieIter]] > tau){
          jProbnu++;
        }
        else if (time[ordering[tieIter]] <= tau && status[ordering[tieIter]] == 2){
          jPart1314+=1.0/(GTiminus[ordering[tieIter]]);
        }
        else if (time[ordering[tieIter]] <= tau && status[ordering[tieIter]] == 1){
          jPart1112+=1.0/(GTiminus[ordering[tieIter]]);
        }
        tieIter++;
      }
      for (int l = i; l < tieIter;l++){
        int ownWeightProbnu = time[ordering[l]] > tau ? 1 : 0;
        Probnu[ordering[l]] = ((double) jProbnu-ownWeightProbnu)/n;
        double ownWeight1314 = time[ordering[l]] <= tau && status[ordering[l]] == 2 ? 1.0/(((double) n)*GTiminus[ordering[l]]) : 0;
        eq1314part[ordering[l]] = jPart1314 - ownWeight1314;
        double ownWeight1112 = time[ordering[l]] <= tau && status[ordering[l]] == 1 ? 1.0/(((double) n)*GTiminus[ordering[l]]) : 0;
        eq1112part[ordering[l]] = jPart1112 - ownWeight1112;
      }
      i = tieIter;
    }
  }
  // Rcout << "Probnu: " << Probnu << "\n\n";
  // Rcout << "eq1314part: " << eq1314part << "\n\n";
  // Rcout << "eq1112part: " << eq1112part << "\n\n\n";
  
  // F_1(tau) = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,1) = Q(T_i <= tau, Delta_i = 1)
  // F_2(tau) = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,2) = Q(T_i <= tau, Delta_i = 2)
  double F1tau = 0, F2tau = 0, eq9term = 0;
  for (int i = 0; i <= firsthit; i++){
    if (status[i] == 1){
      F1tau += 1.0/GTiminus[i];
      eq9term += Probnu[i]/GTiminus[i];
    }
    else if (status[i] == 2){
      F2tau += 1.0/GTiminus[i];
    }
  }
  eq9term = eq9term / n;
  F1tau = F1tau / ((double) n);
  F2tau = F2tau / ((double) n);
  
  double eq16term = Probmu * F1tau;
  double mu1 = F1tau * Probmu / Gtau + F1tau*F2tau;
  double nu1 = auc * mu1;
  double eq10part1{}, eq12part1{},eq14part1{},eq17part1{},eq19part1{};
  double eq10part2 = n*eq9term;
  double eq12part2 = (nu1-1.0/(Gtau) * eq9term)*n*n;
  double eq14part2 = eq12part2;
  double eq17part2 = n*F1tau;
  double eq19part2 = n*F2tau;
  int tieIter = 0;
  // can do while loops together
  while ((tieIter < n) && (time[tieIter] == time[0])) {
    if ((time[tieIter] <= tau) && (status[tieIter]==1)){
      eq10part2 -= Probnu[tieIter] / GTiminus[tieIter];
      eq14part2 -= eq1314part[tieIter] / GTiminus[tieIter];
      eq17part2 -= 1.0 / GTiminus[tieIter];
    }
    else if  ((time[tieIter] <= tau) && (status[tieIter]==2)){
      eq12part2 -= eq1112part[tieIter] / GTiminus[tieIter];
      eq19part2 -= 1.0 / GTiminus[tieIter];
    }
    tieIter++;
  }
  
  int upperTie = tieIter-1;
  double eq8{}, eq9{}, eq10{}, eq11{},eq12{},eq13{},eq14{}, eq15{}, eq16{},eq17{},eq18{},eq19{},eq20{},eq21{}, fihattau{},eq1721part{};
  for (int i = 0; i<n; i++){
    int const j = i > firsthit ? firsthit : i;
    if (utime[sindex[j]] < time[i]){
      fihattau = - MC_term2[sindex[j]];
    }
    else {
      fihattau =  (1-(status[i] != 0))*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
    }
    eq10 = 1.0 / (Gtau*n) * (eq10part1+fihattau*eq10part2);
    eq12 = 1.0 / (n*n) * (eq12part1+(fihattau-1)*eq12part2);
    eq14 = 1.0 / (n*n) * (eq14part1+(fihattau-1)*eq14part2);
    eq1721part = 1.0 / n * (eq17part1+fihattau*eq17part2);
    eq17 = Probmu / Gtau * eq1721part;
    eq19 = F1tau * (1.0 / n * (eq19part1+fihattau*eq19part2) - F2tau);
    eq21 = F2tau * (eq1721part - F1tau);
    
    // fast calculation of eq10 and eq17
    if (upperTie == i){
      int tieIter = i+1;
      while ((tieIter < n) && (time[tieIter] == time[i+1])) {
        if ((time[tieIter] <= tau) && (status[tieIter]==1)){
          eq10part1 -= Probnu[tieIter] * (MC_term2[sindex[i]]) / GTiminus[tieIter];
          eq14part1 -= eq1314part[tieIter] * (MC_term2[sindex[i]]+1.0) / GTiminus[tieIter];
          eq17part1 -= (MC_term2[sindex[i]]) / GTiminus[tieIter];
          eq10part2 -= Probnu[tieIter] / GTiminus[tieIter];
          eq14part2 -= eq1314part[tieIter] / GTiminus[tieIter];
          eq17part2 -= 1.0 / GTiminus[tieIter];
        }
        else if ((time[tieIter] <= tau) && (status[tieIter]==2)){
          eq12part1 -= eq1112part[tieIter] * (MC_term2[sindex[i]]+1.0) / GTiminus[tieIter];
          eq19part1 -= (MC_term2[sindex[i]]) / GTiminus[tieIter];
          eq12part2 -= eq1112part[tieIter] / GTiminus[tieIter];
          eq19part2 -= 1.0 / GTiminus[tieIter];
        }
        tieIter++;
      }
      upperTie = tieIter-1;
    }
    if ((time[i] <= tau) && (status[i] == 1)){
      eq8 = 1.0/(GTiminus[i]*Gtau)*Probnu[i];
      eq15 = 1.0/(GTiminus[i]*Gtau)*Probmu;
      eq13 = 1.0/ (n*GTiminus[i]) * eq1314part[i];
      eq20 = 1.0/GTiminus[i] * F2tau;
    }
    else {
      eq8 = eq13 = eq15 = eq20 = 0;
    }
    eq9 = (fihattau - 2.0)/Gtau * eq9term;
    eq16 = (fihattau - 2.0)/Gtau * eq16term;
    
    if ((time[i] <= tau) && (status[i] == 2)){
      eq11 = 1.0 / (n*GTiminus[i]) * eq1112part[i];
      eq18 = 1.0 / (GTiminus[i]) * F1tau;
    }
    else if (time[i] > tau){
      eq11 = 1.0 / (n*Gtau) * eq1112part[i];
      eq18 = 1.0 / (Gtau) * F1tau;
    }
    else {
      eq11 = eq18 = 0;
    }
    double IFnu = eq8+eq9+eq10+eq11+eq12+eq13+eq14;
    double IFmu = eq15+eq16+eq17+eq18+eq19+eq20+eq21;
    
    ic[i] = (IFnu * mu1 - IFmu * nu1)/(mu1*mu1);
  }
  return ic;
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