// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector getInfluenceFunctionAUCSurvival(NumericVector time,
                                              NumericVector status, 
                                              NumericVector tau,
                                              NumericVector risk,
                                              NumericVector cases,
                                              NumericVector controls,
                                              NumericVector GTiminus,
                                              NumericVector Gtau, 
                                              NumericVector auc) {
  int n = time.size();
  arma::vec ic(n);
  double mu1hat = 0;
  double mutauP = 0;
  for (int i = 0; i < n;i++){
    if (time[i] > tau[0]){
      mu1hat += 1.0;
    }
    if (cases[i] == 1){
      mutauP += 1.0/GTiminus[i];
    }
  }
  mu1hat = mu1hat / n;
  mutauP = mu1hat / Gtau[0] * (mutauP / n);
  double nutauP = auc[0]*mutauP;
  arma::vec nu1hat(n);
  // Rcout << "mutauP: " << mutauP << " and nutauP "<< nutauP << "\n";
  // Rcout << "mu1hat: " << mu1hat;
  for (int i = 0; i < n;i++){
    double temp1 = 0;
    for (int j = 0; j < n; j++){
      if (risk[j] < risk[i] && time[j] > tau[0]){
        temp1+=1.0;
      }
    }
    nu1hat[i] = temp1/n;
  }
  // Rcout << "nu1hat is " << nu1hat << "\n";

  // Thomas code from IC of Nelson-Aalen estimator
  arma::uvec sindex(n,fill::zeros);
  arma::vec utime=unique(time);
  int nu=utime.size();
  arma::vec atrisk(nu);
  arma::vec Cens(nu,fill::zeros);
  arma::vec hazardC(nu,fill::zeros);
  arma::vec MC_term2(nu,fill::zeros);
  //initialize first time point t=0 with data of subject i=0
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
  
  // sindex?
  int firsthit = 0;
  // Rcout << tau << "\n";
  // Rcout << time << "\n";
  // Rcout << (firsthit < n-1) << "\n";
  // Rcout << (time[firsthit+1] <= tau[0]) << "\n";
  while (firsthit < n-1 && time[firsthit+1] <= tau[0]) {
    // Rcout << "inside \n";
    firsthit++;
  }
  // Rcout << "firsthit is " << firsthit;
  
  // main loop, here is the problem
  for (int i=0;i<n;i++){
    double firstTermNum, firstTermDen;
    // Rcout << "cases: "<< cases;
    if (cases[i] == 1){
      firstTermNum = nu1hat[i] / (GTiminus[i]*Gtau[0]);
      firstTermDen = mu1hat / (GTiminus[i]*Gtau[0]);
    }
    else if (controls[i] == 1){
      double nu2hati = 0;
      for (int j = 0; j < n;j++){
        if (risk[i] < risk[j] && cases[j] == 1){
          // Rcout << "j is " << j << "\n";
          nu2hati += 1.0/GTiminus[j];
        }
      }
      nu2hati = 1.0/n * nu2hati;
      // Rcout << "nu2hati" << nu2hati << "\n";
      firstTermNum =  nu2hati * 1.0/Gtau[0];
      firstTermDen =  mutauP / mu1hat;
    }
    else {
      firstTermNum =  0;
      firstTermDen = 0;
    }
    // Rcout << "firstTermNum: " << firstTermNum << " and firstTermDen "<< firstTermDen << "\n";

    // forget about ties for now (otherwise we will go up to something less than i)
    arma::vec MC(i+1);
    MC.zeros();
    // something is wrong here
    for (int t=0;t<=i;t++){
      // Rprintf("i=%d\tt=%d\tsindex[i]=%d\tCens[sindex[i]]=%1.2f\tatrisk[sindex[i]]=%1.2f\tMC_term2[sindex[i]]=%1.2f\tMC_term2[t]=%1.2f\t\n",i,t,sindex[i],Cens[sindex[i]],atrisk[sindex[i]],MC_term2[sindex[i]],MC_term2[t]);
      if (utime[t]<time[i]){
        MC[t] = - MC_term2[t];
      } else{
        MC[t] = (1-status[i])*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
      }
    }
    // arma::vec MC(n);
    // MC.zeros();
    // for (int t=0;t<n;t++){
    //   if (utime[t]<time[i]){
    //     MC[t] = - MC_term2[t];
    //   } else{
    //     MC[t] = (1.0-status[i])*n/atrisk[sindex[i]]- MC_term2[sindex[i]];
    //   }
    // }
    // Rcout << "MC is " << MC  << "\n";

    int j;
    if (i > firsthit){
      j = firsthit;
    }
    else {
      j = i;
    }

    double nu3hati = 0;
    double mu2hat = 0;
    for (int k = 0; k < n;k++){
      if (cases[k] == 1){
        if (k>0 && k < i) {
          nu3hati += nu1hat[k] * MC[k-1] / GTiminus[k];
          mu2hat += MC[k-1] / GTiminus[k];
        }
        else if (k > 0 && k >= i) {
          nu3hati += nu1hat[k] * MC[i] / GTiminus[k];
          mu2hat += MC[i] / GTiminus[k];
        }
      }
    }
    
    nu3hati = 1.0 / n * nu3hati;
    mu2hat = 1.0 / n * mu2hat;
    // Rcout << "nu3hati " << nu3hati << "\n";
    // Rcout << "mu2hat " << mu2hat << "\n";
    double fihattau = MC[j];
    // Rcout << "fihattau " << fihattau << "\n";
    // Rcout << "nutauP " << nutauP << "\n";
    double icnaTermsNum = fihattau * nutauP + (1.0/Gtau[0]) * nu3hati;
    double icnaTermsDen = fihattau * mutauP +mu1hat/Gtau[0] * mu2hat;
    // Rcout << "icnaTermsNum: " << icnaTermsNum << "\n"; //<< " and icnaTermsDen"<< icnaTermsDen << "\n";
    
    ic[i] = ((firstTermNum+icnaTermsNum)*mutauP- nutauP*(firstTermDen+icnaTermsDen))/(mutauP*mutauP);
  }
  return wrap(ic);
}


