// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Calculate influence function for competing risk case/survival case with Nelson-Aalen censoring.
// Author: Johan Sebastian Ohlendorff
// [[Rcpp::export]]
NumericVector getInfluenceFunctionAUC(NumericVector time,
                                      NumericVector status,
                                      double tau,
                                      NumericVector risk,
                                      NumericVector GTiminus,
                                      double Gtau,
                                      double auc,
                                      bool conservative,
                                      bool tiedValues,
                                      bool survival) {
  // Thomas' code from IC of Nelson-Aalen estimator
  //initialize first time point t=0 with data of subject i=0
  int n = time.size();
  NumericVector ic(n);
  arma::uvec sindex(n,fill::zeros);
  arma::vec utime=unique(time);
  int nu=utime.size();
  arma::vec atrisk(nu);
  arma::vec Cens(nu,fill::zeros);
  arma::vec hazardC(nu,fill::zeros);
  arma::vec MC_term2(nu,fill::zeros);
  int t=0;
  double Y = (double) n;

  if (!conservative){
    atrisk[0]=Y;
    Cens[0]=(1-(status[0] != 0));
    hazardC[0]=Cens[0]/Y;
    MC_term2[0]+=hazardC[0];
    //loop through time points until last subject i=(n-1)
    for (int i=1;i<=n;i++) {
      if (i<n && time[i]==time[i-1]){// these are tied values
        Cens[t] +=(1-(status[i] != 0));
        Y-=1;
        sindex[i]=t;    // index pointer from subject i to unique time point t
      }else{
        utime[t]=time[i-1];
        hazardC[t]=Cens[t]/atrisk[t];
        MC_term2[t]=hazardC[t]*n/atrisk[t];
        //initialize next time point with data of current subject i
        if (i<n){
          t++;
          sindex[i]=t;    // index pointer from subject i to unique time point t
          Y-=1;
          atrisk[t]=Y;
          Cens[t]=(1-(status[i] != 0));
        }
      }
    }
    MC_term2 = arma::cumsum(MC_term2);
  }

  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  auto lower = std::upper_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower) -1;
  LogicalVector logicalIndex(n);
  LogicalVector logicalIndex2(n);
  for (int i = 0; i <= firsthit; i++){
    logicalIndex[i] = (status[i] ==1);
    logicalIndex2[i] = (status[i] ==2);
  }
  // subset relevant vectors
  NumericVector risk2 = risk[logicalIndex];
  NumericVector risk3 = risk[logicalIndex2];
  NumericVector GTiminus2 = GTiminus[logicalIndex];
  NumericVector GTiminus3 = GTiminus[logicalIndex2];
  // get ordering according to risk2
  IntegerVector order(risk2.length());
  IntegerVector order2(risk3.length());
  std::iota(order.begin(), order.end(), 0);
  std::iota(order2.begin(), order2.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int x, int y) { return risk2[x] < risk2[y]; });
  std::sort(order2.begin(), order2.end(),
            [&](int x, int y) { return risk3[x] < risk3[y]; });

  // reorder according to rodering of risk2
  risk2 = risk2[order];
  GTiminus2 = GTiminus2[order];
  risk3 = risk3[order2];
  GTiminus3 = GTiminus3[order2];

  // P(tilde{T_i} > tau)
  double Probmu = double ((n-(firsthit+1))) / n;

  //P(X < X_i, tilde{T_i} > tau)
  NumericVector Probnu(n);
  // copy values from risk[firsthit+1:(n-1)] and sort them
  NumericVector risk1 = risk[Range(firsthit+1,n-1)];
  std::sort(risk1.begin(),risk1.end());
  if (!tiedValues){
    for (int i = 0; i <= firsthit;i++){
      int j = 0;
      while (risk1[j] < risk[i] && j < risk1.length()){
        j++;
      }
      Probnu[i] = ((double) j)/n;
    }
  }
  else {
    for (int i = 0; i <= firsthit;i++){
      int j = 0;
      for (int k = firsthit+1; k<n;k++){
        if (risk[k]==risk[i]){
          j++;
        }
      }
      //P(X = X_i, tilde{T_i} > tau)
      Probnu[i] = ((double) j)/n;
      // if (Probnu[i]!=0){
      //   Rcout << Probnu[i] << "\n";
      // }
    }
  }


  // F_1(tau) = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,1) = Q(T_i <= tau, Delta_i = 1)
  // F_2(tau) = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,2) = Q(T_i <= tau, Delta_i = 2)
  double F1tau, F2tau, eq9term = 0;
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
  NumericVector eq1314part(firsthit+1);
  NumericVector eq1112part(n);
  int begginingIndexEq1112;
  if (!tiedValues){
    if (!survival){
      for (int i = 0; i <= firsthit;i++){
        int j = 0;
        while (j < risk3.length()  && risk[i] > risk3[j]){
          eq1314part[i] += 1.0/GTiminus3[j];
          j++;
        }
      }
      begginingIndexEq1112 = 0;
    }
    else {
      begginingIndexEq1112 = firsthit+1;
    }
  }
  else {
    if (!survival){
      for (int i = 0; i <= firsthit;i++){
        for (int k = 0; k <= firsthit; k++){
          if (risk[k]==risk[i] && k!=i){
            eq1314part[i]+= 1.0/GTiminus[k];
          }
        }
      }
      begginingIndexEq1112 = 0;
    }
    else {
      begginingIndexEq1112 = firsthit+1;
    }
  }
  if (!tiedValues){
    for (int i = begginingIndexEq1112; i < n;i++){
      int j = risk2.length()-1;
      while (j >= 0 && risk[i] < risk2[j]){
        eq1112part[i] += 1.0/GTiminus2[j];
        j--;
      }
    }
  }
  else {
    for (int i = begginingIndexEq1112; i < n;i++){
      for (int j = 0; j <= firsthit; j++){
        if (risk[i]==risk[j] && j!=i){
          eq1112part[i] += 1.0/GTiminus[j];
        }
      }
    }
  }



  double eq10part1, eq12part1,eq14part1,eq17part1,eq19part1 = 0;
  double eq10part2 = n*eq9term;
  double eq12part2 = (nu1-1.0/(Gtau) * eq9term)*n*n;
  double eq14part2 = eq12part2;
  double eq17part2 = n*F1tau;
  double eq19part2 = n*F2tau;
  int tieIter = 0;
  // can do while loops together
  while ((time[tieIter] == time[0]) && (tieIter < n)) {
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
  double eq8, eq9, eq10, eq11,eq12,eq13,eq14, eq15, eq16,eq17,eq18,eq19,eq20,eq21, fihattau,j,eq1721part;
  for (int i = 0; i<n; i++){
    if (!conservative){
      if (i > firsthit){
        j = firsthit;
      }
      else {
        j = i;
      }
      if (utime[sindex[j]] < time[i]){
        fihattau = - MC_term2[sindex[j]];
      }
      else {
        fihattau =  (1-(status[i] != 0))*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
      }
      // fast calculation of eq10 and eq17
      eq10 = 1.0 / (Gtau*n) * (eq10part1+fihattau*eq10part2);
      eq12 = 1.0 / (n*n) * (eq12part1+(fihattau-1)*eq12part2);
      eq14 = 1.0 / (n*n) * (eq14part1+(fihattau-1)*eq14part2);
      eq1721part = 1.0 / n * (eq17part1+fihattau*eq17part2);
      eq17 = Probmu / Gtau * eq1721part;
      eq19 = F1tau * (1.0 / n * (eq19part1+fihattau*eq19part2) - F2tau);
      eq21 = F2tau * (eq1721part - F1tau);
      if (upperTie <= i){
        int tieIter = i+1;
        while ((time[tieIter] == time[i+1]) && (tieIter < n)) {
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
    }
    else {
      fihattau = eq10 = eq17 = 0;
      eq12 = eq14 = 1.0/(Gtau) * eq9term-nu1;
      eq19 = eq21 = -F1tau * F2tau;
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

    if ((time[i] <= tau) & (status[i] == 2)){
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
    // Rcout << "eq8: " << eq8 <<"\n";
    // Rcout << "eq9: " << eq9 <<"\n";
    // Rcout << "eq10: " << eq10 <<"\n";
    // Rcout << "eq11: " << eq11 <<"\n";
    // Rcout << "eq12: " << eq12 <<"\n";
    // Rcout << "eq13: " << eq13 <<"\n";
    // Rcout << "eq14: " << eq14 <<"\n";
    // Rcout << "eq15: " << eq15 <<"\n";
    // Rcout << "eq16: " << eq16 <<"\n";
    // if (IFnu!=0){
    //   Rcout << IFnu;
    // }
    double IFmu = eq15+eq16+eq17+eq18+eq19+eq20+eq21;
    ic[i] = (IFnu * mu1 - IFmu * nu1)/(mu1*mu1);
  }
  return ic;
}

// [[Rcpp::export]]
NumericVector getInfluenceFunctionAUCConservative(NumericVector time,
                                      NumericVector status,
                                      double tau,
                                      NumericVector risk,
                                      NumericVector GTiminus,
                                      double Gtau,
                                      double auc,
                                      bool survival) {
  // Thomas' code from IC of Nelson-Aalen estimator
  //initialize first time point t=0 with data of subject i=0
  int n = time.size();
  NumericVector ic(n);

  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  auto lower = std::upper_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower) -1;
  LogicalVector logicalIndex(n);
  LogicalVector logicalIndex2(n);
  for (int i = 0; i <= firsthit; i++){
    logicalIndex[i] = (status[i] ==1);
    logicalIndex2[i] = (status[i] ==2);
  }
  // subset relevant vectors
  NumericVector risk2 = risk[logicalIndex];
  NumericVector risk3 = risk[logicalIndex2];
  NumericVector GTiminus2 = GTiminus[logicalIndex];
  NumericVector GTiminus3 = GTiminus[logicalIndex2];
  // get ordering according to risk2
  IntegerVector order(risk2.length());
  IntegerVector order2(risk3.length());
  std::iota(order.begin(), order.end(), 0);
  std::iota(order2.begin(), order2.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int x, int y) { return risk2[x] < risk2[y]; });
  std::sort(order2.begin(), order2.end(),
            [&](int x, int y) { return risk3[x] < risk3[y]; });

  // reorder according to rodering of risk2
  risk2 = risk2[order];
  GTiminus2 = GTiminus2[order];
  risk3 = risk3[order2];
  GTiminus3 = GTiminus3[order2];

  // P(tilde{T_i} > tau)
  double Probmu = double ((n-(firsthit+1))) / n;

  //P(X < X_i, tilde{T_i} > tau)
  NumericVector Probnu(n);
  // copy values from risk[firsthit+1:(n-1)] and sort them
  NumericVector risk1 = risk[Range(firsthit+1,n-1)];
  std::sort(risk1.begin(),risk1.end());
  for (int i = 0; i < n;i++){
    int j = 0;
    while (risk1[j] < risk[i] && j < risk1.length()){
      j++;
    }
    Probnu[i] = ((double) j)/n;
  }

  // F_1(tau) = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,1) = Q(T_i <= tau, Delta_i = 1)
  // F_2(tau) = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,2) = Q(T_i <= tau, Delta_i = 2)
  double F1tau, F2tau, eq9term = 0;
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

  double eq8, eq9, eq10, eq11,eq12,eq13,eq14, eq15, eq16,eq17,eq18,eq19,eq20,eq21, fihattau;
  for (int i = 0; i<n; i++){
    fihattau = eq10 = eq17 = 0;
    eq12 = eq14 = 1.0/(Gtau) * eq9term-nu1;
    eq19 = eq21 = -F1tau * F2tau;
    if ((time[i] <= tau) && (status[i] == 1)){
      eq8 = 1.0/(GTiminus[i]*Gtau)*Probnu[i];
      eq15 = 1.0/(GTiminus[i]*Gtau)*Probmu;
      double temp = 0;
      int j = 0;
      while (j < risk3.length()  && risk[i] > risk3[j]){
        temp += 1.0/GTiminus3[j];
        j++;
      }
      eq13 = 1.0/ (n*GTiminus[i]) * temp;
      eq20 = 1.0/GTiminus[i] * F2tau;
    }
    else {
      eq8 = eq13 = eq15 = eq20 = 0;
    }
    eq9 = (fihattau - 2.0)/Gtau * eq9term;
    eq16 = (fihattau - 2.0)/Gtau * eq16term;

    if ((time[i] <= tau) & (status[i] == 2)){
      double temp = 0;
      int j = risk2.length()-1;
      while (j >= 0 && risk[i] < risk2[j]){
        temp += 1.0/GTiminus2[j];
        j--;
      }
      eq11 = 1.0 / (n*GTiminus[i]) * temp;
      eq18 = 1.0 / (GTiminus[i]) * F1tau;
    }
    else if (time[i] > tau){
      double temp = 0;
      int j = risk2.length()-1;
      while (j >= 0 && risk[i] < risk2[j]){
        temp += 1.0/GTiminus2[j];
        j--;
      }
      eq11 = 1.0 / (n*Gtau) * temp;
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


