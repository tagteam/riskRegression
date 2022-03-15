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
                                                   double auc) {

  // Thomas code from IC of Nelson-Aalen estimator
  //initialize first time point t=0 with data of subject i=0
  int n = time.size();
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

  // find first index such that k such that tau[k] <= tau but tau[k+1] > tau
  auto lower = std::lower_bound(time.begin(), time.end(), tau);
  int firsthit = std::distance(time.begin(), lower)-1;
  // Calculate \hat{mu}_\tau(P),\hat{mu}_1, \hat{nu}_\tau(P) and \hat{nu}_1 (see formulas)
  NumericVector ic(n);
  // needed as part of another calculation
  double mutauPsurvival2 = 0;
  double mutauPsurvival3 = 0;
  double nu11hat = 0;
  double nu13hat = 0;
  double mu13hat = 0;
  for (int i = 0; i <= firsthit; i++){
    if (status[i] == 1){
      mutauPsurvival2 += 1.0/GTiminus[i];
      mu13hat += (MC_term2[sindex[i]-1]-1.0)/GTiminus[i];
      double counter1 = 0.0;
      for (int j = 0; j <= firsthit; j++){
        if (status[j] == 2 && risk[i] > risk[j]) {
          counter1 += 1.0/GTiminus[j];
        }
      }
      nu11hat += 1.0/GTiminus[i] * counter1;
      // hmm, maybe does not work for ties, but who cares about those?
      nu13hat += (MC_term2[sindex[i]-1]-1.0)/GTiminus[i] * counter1;
    }
    else if (status[i] == 2){
      mutauPsurvival3 += 1.0/GTiminus[i];
    }
  }
  Rcout << "works until here";
  nu11hat = nu11hat / (n*n);
  nu13hat = nu13hat / (n*n);
  mutauPsurvival2 = mutauPsurvival2 / n;
  mutauPsurvival3 = mutauPsurvival3 / n;
  mu13hat = mu13hat / n *mutauPsurvival3;
  double mu1hat = double ((n-(firsthit+1))) / n;
  double mu11hat = mutauPsurvival2 * mutauPsurvival3;
  double mutauPsurvival = mu1hat / Gtau * mutauPsurvival2;
  double mutauP =mutauPsurvival + mu11hat;
  double nutauP = auc*mutauP;
  // recalculate

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

  double nutauPsurvival = 0.0; // = auc*mutauPsurvival previously
  for (int i = 0; i <= firsthit;i++){
    nutauPsurvival += nu1hat[i] / GTiminus[i];
  }
  nutauPsurvival = 1.0/(n*Gtau) * nutauPsurvival;



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

  double nu3hati1 = 0;
  double mu2hat1 = 0;
  double nu3hati2 = 0;
  double mu2hat2 = 0;

  for (int k = 1; k <= firsthit;k++){
    if (status[k] == 1){
      nu3hati2 += nu1hat[k] / GTiminus[k];
      mu2hat2 += 1.0 / GTiminus[k];
    }
  }
  double nu3hati, mu2hat;
  for (int i=0;i<n;i++){
    double firstTermNum, firstTermDen, nu12hat,mu12hat;
    if (time[i] <= tau){
      mu12hat = 1.0/GTiminus[i] * mutauPsurvival3;
      nu12hat = 0;
      for (int j = 0; j <= firsthit; j++){
        if (status[j] == 2 && risk[j] < risk[i]){
          nu12hat += 1.0/GTiminus[j];
        }
      }
      nu12hat = nu12hat / (n*GTiminus[i]);
      if (status[i] == 1){
        firstTermNum = nu1hat[i] / (GTiminus[i]*Gtau);
        firstTermDen = mu1hat / (GTiminus[i]*Gtau);
      }
      else {
        firstTermNum =  0;
        firstTermDen = 0;
      }
    }
    else if (time[i] > tau){
      mu12hat = 0;
      nu12hat = 0;
      double nu2hati = 0;
      int j = risk2.length()-1;
      while (j >= 0 && risk[i] < risk2[j]){
        nu2hati += 1.0/GTiminus2[j];
        j--;
      }
      nu2hati = 1.0/n * nu2hati;
      firstTermNum =  nu2hati * 1.0/Gtau;
      firstTermDen =  mutauPsurvival / mu1hat;
    }
    else {
      mu12hat = 0;
      nu12hat = 0;
      firstTermNum =  0;
      firstTermDen = 0;
    }
    double fihattau = (1-status[i])*n/atrisk[sindex[i]]- MC_term2[sindex[i]];

    if (i==0){
      nu3hati = nu3hati1;
      mu2hat = mu2hat1;
      nu3hati = 1.0 / n * nu3hati;
      mu2hat = 1.0 / n * mu2hat;
    }
    else if (i==1){
      nu3hati = nu3hati1+nu3hati2*fihattau;
      mu2hat = mu2hat1+mu2hat2*fihattau;
      nu3hati = 1.0 / n * nu3hati;
      mu2hat = 1.0 / n * mu2hat;
    }
    if (i>1 && (i-1 <= firsthit)){
      if (status[i-1] == 1){
        nu3hati1 -= nu1hat[i-1] * (MC_term2[sindex[i-1]-1]) / GTiminus[i-1];
        mu2hat1 -= (MC_term2[sindex[i-1]-1]) / GTiminus[i-1];
        nu3hati2 -= nu1hat[i-1] / GTiminus[i-1];
        mu2hat2 -= 1.0 / GTiminus[i-1];
      }
      nu3hati = nu3hati1+nu3hati2*fihattau;
      mu2hat = mu2hat1+mu2hat2*fihattau;
      nu3hati = 1.0 / n * nu3hati;
      mu2hat = 1.0 / n * mu2hat;
    }

    double icnaTermsNum = fihattau * nutauPsurvival + (1.0/Gtau) * nu3hati;
    double icnaTermsDen = fihattau * mutauPsurvival + mu1hat/Gtau * mu2hat;
    ic[i] = ((firstTermNum+icnaTermsNum-2*nutauPsurvival+1*(status[i] == 2 && time[i] > tau)*firstTermNum +(fihattau-1)*nu11hat + nu12hat+nu13hat)*mutauP- nutauP*(firstTermDen+icnaTermsDen-2*mutauPsurvival)+1*(status[i] == 2 && time[i] > tau)*(1.0/Gtau * mutauPsurvival2) +(fihattau-1)*mu11hat + mu12hat+mu13hat)/(nutauP*nutauP);
  }
  return ic;
}

