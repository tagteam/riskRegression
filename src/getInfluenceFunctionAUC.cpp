// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Calculate influence function for survival case with Nelson-Aalen censoring.
// Author: Johan Sebastian Ohlendorff
// [[Rcpp::export]]
NumericVector getInfluenceFunctionAUCSurvival(NumericVector time,
                                              NumericVector status,
                                              double tau,
                                              NumericVector risk,
                                              NumericVector GTiminus,
                                              double Gtau,
                                              double auc) {

  int n = time.size();
  // find first index such that k such that time[k] <= tau but time[k+1] > tau
  auto lower = std::lower_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower)-1;
  // Calculate \hat{mu}_\tau(P),\hat{mu}_1, \hat{nu}_\tau(P) and \hat{nu}_1 (see formulas)
  NumericVector ic(n);
  double mutauP = 0;
  for (int i = 0; i <= firsthit; i++){
    if (status[i] == 1){
      mutauP += 1.0/GTiminus[i];
    }
  }
  double mu2hat2 = mutauP;

  double mu1hat = double ((n-(firsthit+1))) / n;
  mutauP = mu1hat / Gtau * (mutauP / n);
  double nutauP = auc*mutauP;
  double nu3hati2 = nutauP * Gtau*n;

  NumericVector nu1hat(n);

  // copy values from risk[firsthit+1:(n-1)] and sort them
  NumericVector risk1 = risk[Range(firsthit+1,n-1)];
  std::sort(risk1.begin(),risk1.end());

  for (int i = 0; i < n;i++){
    int j = 0;
    while (risk1[j] < risk[i] && j < risk1.length()){
      j++;
    }
    nu1hat[i] = ((double) j)/n;
  }

  LogicalVector logicalIndex(n);
  for (int i = 0; i <= firsthit; i++){
    logicalIndex[i] = (status[i] ==1);
  }
  // subset relevant vectors
  NumericVector risk2 = risk[logicalIndex];
  NumericVector GTiminus2 = GTiminus[logicalIndex];
  // get ordering according to risk2
  IntegerVector order(risk2.length());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](int x, int y) { return risk2[x] < risk2[y]; });
  // reorder according to rodering of risk2
  risk2 = risk2[order];
  GTiminus2 = GTiminus2[order];

  // Thomas code from IC of Nelson-Aalen estimator
  //initialize first time point t=0 with data of subject i=0
  // this can be optimized by saving information from here.
  arma::uvec sindex(n,fill::zeros);
  arma::vec utime=unique(time);
  int nu=utime.size();
  arma::vec atrisk(nu);
  arma::vec Cens(nu,fill::zeros);
  arma::vec hazardC(nu,fill::zeros);
  arma::vec MC_term2(nu,fill::zeros);

  int t=0;
  double Y = (double) n;
  atrisk[0]=Y;
  Cens[0]=(1-status[0]);
  hazardC[0]=Cens[0]/Y;
  MC_term2[0]+=hazardC[0];
  //loop through time points until last subject i=(n-1)
  for (int i=1;i<=n;i++) {
    if (i<n && time[i]==time[i-1]){// these are tied values
      Cens[t] +=(1-status[i]);
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
        Cens[t]=(1-status[i]);
      }
    }
  }
  MC_term2 = arma::cumsum(MC_term2);

  double nu3hati1 = 0;
  double mu2hat1 = 0;
  double nu3hati, mu2hat;
  int tieIter = 0;
  // can do while loops together
  while ((time[tieIter] == time[0]) & (tieIter < n)) {
    if ((time[tieIter] <= tau) & (status[tieIter]==1)){
      nu3hati2 -= nu1hat[tieIter] / GTiminus[tieIter];
      mu2hat2 -= 1.0 / GTiminus[tieIter];
    }
    tieIter++;
  }
  int upperTie = tieIter-1;

  for (int i=0;i<n;i++){
    double firstTermNum, firstTermDen;
    if (time[i] <= tau && status[i] == 1){
      firstTermNum = nu1hat[i] / (GTiminus[i]*Gtau);
      firstTermDen = mu1hat / (GTiminus[i]*Gtau);
    }
    else if (time[i] > tau){
      double nu2hati = 0;
      int j = risk2.length()-1;
      while (j >= 0 && risk[i] < risk2[j]){
        nu2hati += 1.0/GTiminus2[j];
        j--;
      }
      nu2hati = 1.0/n * nu2hati;
      firstTermNum =  nu2hati * 1.0/Gtau;
      firstTermDen =  mutauP / mu1hat;
    }
    else {
      firstTermNum =  0;
      firstTermDen = 0;
    }
    double fihattau = (1-status[i])*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
    nu3hati = 1.0 / n * (nu3hati1+fihattau*nu3hati2);
    mu2hat = 1.0 / n * (mu2hat1+fihattau*mu2hat2);
    if (upperTie <= i){
      int tieIter = i+1;
      while ((time[tieIter] == time[i+1]) & (tieIter < n)) {
        if ((time[tieIter] <= tau) & (status[tieIter]==1)){
          nu3hati1 -= nu1hat[tieIter] * (MC_term2[sindex[i]]) / GTiminus[tieIter];
          mu2hat1 -= (MC_term2[sindex[i]]) / GTiminus[tieIter];
          nu3hati2 -= nu1hat[tieIter] / GTiminus[tieIter];
          mu2hat2 -= 1.0 / GTiminus[tieIter];
        }
        tieIter++;
      }
      upperTie = tieIter-1;
    }

    double icnaTermsNum = fihattau * nutauP + (1.0/Gtau) * nu3hati;
    double icnaTermsDen = fihattau * mutauP + mu1hat/Gtau * mu2hat;
    ic[i] = ((firstTermNum+icnaTermsNum)*mutauP- nutauP*(firstTermDen+icnaTermsDen))/(mutauP*mutauP);
  }
  return ic;
}

// Calculate influence function for competing risk case with Nelson-Aalen censoring.
// Author: Johan Sebastian Ohlendorff
// [[Rcpp::export]]
NumericVector getInfluenceFunctionAUCCompetingRisk(NumericVector time,
                                                   NumericVector status,
                                                   double tau,
                                                   NumericVector risk,
                                                   NumericVector GTiminus,
                                                   double Gtau,
                                                   double auc,
                                                   NumericVector statusCensoring,
                                                   bool conservative) {


  // Thomas code from IC of Nelson-Aalen estimator
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
    Cens[0]=(1-statusCensoring[0]);
    hazardC[0]=Cens[0]/Y;
    MC_term2[0]+=hazardC[0];
    //loop through time points until last subject i=(n-1)
    for (int i=1;i<=n;i++) {
      if (i<n && time[i]==time[i-1]){// these are tied values
        Cens[t] +=(1-statusCensoring[i]);
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
          Cens[t]=(1-statusCensoring[i]);
        }
      }
    }
    MC_term2 = arma::cumsum(MC_term2);
  }

  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  auto lower = std::upper_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower) -1 ; // was -1
  // Rcout << "firsthit " << firsthit << "\n";
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
  // Rcout << "Probmu " <<  Probmu << "\n";

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
  // Rcout << "probnu" << Probnu<< "\n";

  // F_1(tau) = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,1) = Q(T_i <= tau, Delta_i = 1)
  // F_2(tau) = int 1_{t \leq tau} 1/ hat{G(t-)} dP(t,2) = Q(T_i <= tau, Delta_i = 2)
  double F1tau = 0;
  double F2tau = 0;
  double eq9term = 0 ;
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
  // Rcout << "eq9term" << eq9term << "\n";
  F1tau = F1tau / ((double) n);
  F2tau = F2tau / ((double) n);
  double eq16term = Probmu * F1tau;
  // Rcout << "eq16term" << eq16term << "\n";
  double mu1 = F1tau * Probmu / Gtau + F1tau*F2tau;
  double nu1 = auc * mu1;
  // Rcout << "mu1" << mu1 << " nu1 "<< nu1 << "\n";
  // fast calculation of eq10 and eq17 - part 1
  double eq10part2 = 0;
  for (int i = 0; i <= firsthit;i++){
    if (status[i] == 1){
      eq10part2 += Probnu[i] / GTiminus[i];
    }
  }
  double eq10part1 = 0;
  double eq17part1 = 0;
  double eq19part1 = 0;
  double eq17part2 = n*F1tau;
  double eq19part2 = n*F2tau;
  double eq10part, eq17part, eq19part;
  int tieIter = 0;
  // can do while loops together
  while ((time[tieIter] == time[0]) && (tieIter < n)) {
    if ((time[tieIter] <= tau) && (status[tieIter]==1)){
      eq10part2 -= Probnu[tieIter] / GTiminus[tieIter];
      eq17part2 -= 1.0 / GTiminus[tieIter];
    }
    else if  ((time[tieIter] <= tau) && (status[tieIter]==2)){
      eq19part2 -= 1.0 / GTiminus[tieIter];
    }
    tieIter++;
  }
  int upperTie = tieIter-1;
  // Rcout << "MC_term2 is: "<< MC_term2 << "\n";
  // Rcout << "sindex is: " << sindex << "\n";
  // Rcout << "utime is: " << utime << "\n";

  double eq8, eq9, eq10, eq11,eq12,eq13,eq14, eq15, eq16,eq17,eq18,eq19,eq20,eq21,eq24, fihattau, j;
  if (conservative){
    eq24 = nu1-1.0/(Gtau) * eq9term;
  }
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
        fihattau =  (1-statusCensoring[i])*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
      }
      // fast calculation of eq10 and eq17
      eq10part = 1.0 / n * (eq10part1+fihattau*eq10part2);
      eq17part = 1.0 / n * (eq17part1+fihattau*eq17part2);
      eq19part =  1.0 / n * (eq19part1+fihattau*eq19part2);
      eq10 = 1.0 / Gtau * eq10part;
      eq17 = Probmu / Gtau * eq17part;
      if (upperTie <= i){
        int tieIter = i+1;
        while ((time[tieIter] == time[i+1]) && (tieIter < n)) {
          if ((time[tieIter] <= tau) && (status[tieIter]==1)){
            eq10part1 -= Probnu[tieIter] * (MC_term2[sindex[i]]) / GTiminus[tieIter];
            eq17part1 -= (MC_term2[sindex[i]]) / GTiminus[tieIter];
            eq10part2 -= Probnu[tieIter] / GTiminus[tieIter];
            eq17part2 -= 1.0 / GTiminus[tieIter];
          }
          else if ((time[tieIter] <= tau) && (status[tieIter]==2)){
            eq19part1 -= (MC_term2[sindex[i]]) / GTiminus[tieIter];
            eq19part2 -= 1.0 / GTiminus[tieIter];
          }
          tieIter++;
        }
        upperTie = tieIter-1;
      }
      double eq12part = 0;
      double eq14part = 0;
      // one should be able to do this more efficiently
      for (int k = 0; k <= firsthit; k++){
        if ((time[k] <= tau) && (status[k] ==1)){
          double temp1 = 0;
          double temp2 = 0;
          for (int j = 0; j <= firsthit; j++){
            if ((risk[k] > risk[j]) && (status[j] == 2)){
              double fval = 0;
              if (sindex[j] == 0){
                fval = 0;
              }
              if (utime[sindex[j]] <= time[i]){
                fval = - MC_term2[sindex[j]];
              }
              else {
                fval = (1-statusCensoring[i])*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
              }
              temp1 += (fval-1.0)/GTiminus[j];
              temp2 += 1.0/GTiminus[j];
            }
          }
          eq12part += 1.0/GTiminus[k] * temp1;
          double fval = 0;
          if (sindex[k] == 0){
            fval = 0;
          }
          if (utime[sindex[k]] <= time[i]){
            fval = - MC_term2[sindex[k]];
          }
          else {
            fval = (1-statusCensoring[i])*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
          }
          eq14part += (fval-1.0)/GTiminus[k] * temp2;
        }
      }
      eq12 = 1.0 / (double (n*n)) * eq12part;
      eq14 = 1.0 / (double (n*n)) * eq14part;
      eq19 = F1tau * (eq19part - F2tau);
      eq21 = F2tau * (eq17part - F1tau);
    }
    else {
      fihattau = 0;
      eq10 = 0;
      eq17 = 0;
      eq12 = eq24;
      eq14 = eq24;
      eq19 = F1tau * (0 - F2tau);
      eq21 = F2tau * (0 - F1tau);
    }
    // Rcout << "fihattau " << fihattau << "\n";
    if ((time[i] <= tau) && (status[i] == 1)){
      eq8 = 1.0/(GTiminus[i]*Gtau)*Probnu[i];
      eq15 = 1.0/(GTiminus[i]*Gtau)*Probmu;
      double eq13part = 0;
      int j = 0;
      while (j < risk3.length()  && risk[i] > risk3[j]){
        eq13part += 1.0/GTiminus3[j];
        j++;
      }
      eq13 = 1.0/ (n*GTiminus[i]) * eq13part;
      eq20 = 1.0/GTiminus[i] * F2tau;
    }
    else {
      eq8 = 0;
      eq15 = 0;
      eq13 = 0;
      eq20 = 0;
    }
    eq9 = (fihattau - 2.0)/Gtau * eq9term;
    eq16 = (fihattau - 2.0)/Gtau * eq16term;

    if ((time[i] <= tau) & (status[i] == 2)){
      double eq11part = 0;
      int j = risk2.length()-1;
      while (j >= 0 && risk[i] < risk2[j]){
        eq11part += 1.0/GTiminus2[j];
        j--;
      }
      eq11 = 1.0 / (n*GTiminus[i]) * eq11part;
      eq18 = 1.0 / (GTiminus[i]) * F1tau;
    }
    else if (time[i] > tau){
      double eq11part = 0;
      int j = risk2.length()-1;
      while (j >= 0 && risk[i] < risk2[j]){
        eq11part += 1.0/GTiminus2[j];
        j--;
      }
      eq11 = 1.0 / (n*Gtau) * eq11part;
      eq18 = 1.0 / (Gtau) * F1tau;
    }
    else {
      eq11 = 0;
      eq18 = 0;
    }

    // if (i < 19){
    //   Rcout << "i is " << i+1 << "\n";
    //   // Rcout << eq12part;
    //   Rcout << eq8 << " " <<  eq9 << " " << eq10 << " " << eq11 <<" " <<  eq12 <<" " <<  eq13 <<" " <<  eq14 << "\n";
    //   Rcout << eq15 <<" " <<  eq16 <<" " <<  eq17 <<" " <<  eq18 <<" " <<  eq19<<" " <<  eq20 <<" " << eq21 << "\n";
    // }
    double IFnu = eq8+eq9+eq10+eq11+eq12+eq13+eq14;
    double IFmu = eq15+eq16+eq17+eq18+eq19+eq20+eq21;
    ic[i] = (IFnu * mu1 - IFmu * nu1)/(mu1*mu1);
  }
  return ic;
}
