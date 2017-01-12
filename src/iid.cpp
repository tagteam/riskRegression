// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
List calcE_cpp(const NumericVector& eventtime, 
               const NumericVector& status,
               const NumericVector& eXb,
               const arma::mat& X,
               int p, bool add0){
  
  int nObs = eventtime.size();
  
  // define times
  NumericVector t = eventtime[status>0];
  t = unique(t);
  std::sort(t.begin(),t.end());
  
  if(add0){t.push_back(eventtime[eventtime.size()-1]+1e-12);}
  int nTime = t.size();
  
  // intialisation
  NumericVector resS0(nTime, 0.0);
  NumericMatrix resS1(nTime,p);
  std::fill(resS1.begin(), resS1.end(), 0.0);
  NumericMatrix resE(nTime,p);
  std::fill(resE.begin(), resE.end(), 0.0);
  
  int iTime = nTime-1; 
  while(iTime >= 0 && eventtime[nObs-1]<t[iTime]){
    iTime--;
  }
  double S0=0.0;
  NumericVector S1(p,0.0);
  
  // loop over observations (must be sorted by time)
  for(int iObs=nObs-1;iObs>=0;iObs--){
    
    S0 += eXb[iObs];
    for(int iX=0;iX<p;iX++){
      S1[iX] += X(iObs,iX) * eXb[iObs];
    }
    
    // either the next eventtime is below the time horizon or it is the last event and the eventtime equals the time horizon
    while(iTime >= 0  && ((iObs > 0 && eventtime[iObs-1]<t[iTime]) || (iObs==0 && eventtime[iObs]==t[iTime])) ){
      resS0[iTime] = S0;
      resS1.row(iTime) = S1;
      if(S0>0){
        resE.row(iTime) = S1/S0;
      }// else already initialized at 0
      iTime--;
    }
    
    if(iTime < 0){ break; }
    
  }
  
  return(List::create(Named("E")  = resE,
                      Named("S1")  = resS1,
                      Named("S0")  = resS0,
                      Named("Utime1") = t));
}

// [[Rcpp::export]]
arma::mat ICbeta_cpp(const NumericVector& newT, const NumericVector& neweXb, const arma::mat& newX, const NumericVector& newStatus, const IntegerVector& newIndexJump, 
                     const NumericVector& S01, const arma::mat& E1, const NumericVector& time1, const arma::mat& iInfo,
                     int p){
  
  arma::mat ICbeta;
  int nObs = newT.size();
  
  if(p==0){
    ICbeta.resize(nObs, 1);
    ICbeta.fill(0);
  }else{
    ICbeta.resize(nObs, p);
    ICbeta.fill(NA_REAL);
    
    // initialisation
    int nTime1 = time1.size();
    arma::colvec Score(p);
    double iS0_iter; //  = \sum_tj<tnew delta_j/S0_j
    NumericVector E_iS0_iter(p);// = \sum_tj<tnew E_j delta_j/S0_j
    
    int iObs = 0;
    iS0_iter = 0;
    E_iS0_iter.fill(0);
    while(iObs < nObs && time1[0]>newT[iObs]){ // before the first event
      
      // compute the score
      for(int iX=0;iX<p;iX++){
        Score[iX] = newStatus[iObs] * (newX(iObs,iX)-E1(newIndexJump[iObs],iX) );
      }
      
      ICbeta.row(iObs) = (iInfo * Score).t();
      
      iObs++;
      
    }
    
    for(int iTime1=0 ; iTime1<nTime1 ; iTime1++){
    
      // update the sum
      iS0_iter += 1/S01[iTime1];
      for(int iX=0;iX<p;iX++){
        E_iS0_iter[iX] += E1(iTime1,iX) / S01[iTime1];
      }
      
      // store the value of the influence function
      while(iObs < nObs && ((iTime1 < (nTime1-1) && time1[iTime1+1]>newT[iObs]) || (iTime1==(nTime1-1) && time1[iTime1]==newT[iObs])) ){
        
        // compute the score
        for(int iX=0;iX<p;iX++){
          Score[iX] = newStatus[iObs] * (newX(iObs,iX)-E1(newIndexJump[iObs],iX) ) - neweXb[iObs] * (newX(iObs,iX)*iS0_iter - E_iS0_iter[iX]);
        }
        
        ICbeta.row(iObs) = (iInfo * Score).t();
        iObs++;
      }
      
      if(iObs == nObs){ break; }
      
    }
    
    
  }
  
  return(ICbeta);
}


// [[Rcpp::export]]
arma::mat IClambda0_cpp(const NumericVector& tau, const arma::mat& ICbeta,
                        const NumericVector& newT, const NumericVector& neweXb, const NumericVector& newStatus, const IntegerVector& newStrata, const IntegerVector& newIndexJump, 
                        const NumericVector& S01, const arma::mat& E1, const NumericVector& time1, const NumericVector& lambda0,
                        int p, int strata){
  
  int nObs = newT.size();
  int nTau = tau.size();
  int nTime1 = time1.size();
  arma::mat IClambda0(nObs, nTau);
  IClambda0.fill(NA_REAL);
  
  // Compute delta_iS0
  NumericVector delta_iS0(nObs);
  for(int iObs=0; iObs<nObs ; iObs++){
    delta_iS0[iObs] = newStatus[iObs]/S01[newIndexJump[iObs]];
  } 
  
  // Compute  Elambda0 and lamba0_iS0
  arma::mat Elambda0(p, nTau);
  NumericVector lamba0_iS0(nTime1);
  
  colvec Elambda0_iter(p); Elambda0_iter.fill(0);
  int iTau = 0;
  
  for(int iTime1 = 0; iTime1 < nTime1; iTime1++){
    
    if(p>0){
      for(int iX = 0; iX < p; iX++){
        Elambda0_iter[iX] += E1(iTime1,iX) * lambda0[iTime1];
      }
    }
    
    // update lamba0_iS0
    if(iTime1 == 0){
      lamba0_iS0[0] = lambda0[0]/S01[0];
    }else{
      lamba0_iS0[iTime1] = lamba0_iS0[iTime1-1] + lambda0[iTime1]/S01[iTime1];  
    }
    
    // store in Elambda0 when time is tau
    while(iTau < nTau && ((iTime1 < (nTime1-1) && time1[iTime1+1]>tau[iTau]) || (iTime1==(nTime1-1) && time1[iTime1]==tau[iTau])) ){
      if(p>0){Elambda0.col(iTau) = Elambda0_iter;}
      iTau++;
     }
    
    if(iTau == nTau){ break; }
  }
  
  // main loop
  int iTau0 = 0;
  
  while((iTau0 < nTau) && time1[0]>tau[iTau0]){ // before the first event
    
    for(int iObs=0; iObs<nObs ; iObs++){
      if(newT[iObs]<=tau[iTau0]){
        IClambda0(iObs,iTau0) = delta_iS0[iObs];
      }
    }
    iTau0++;
    
  }
  
  int index_newT_time1; // position of the minimum between t and t_train in lamba0_iS0
  IntegerVector Vindex_tau_time1(nTau);
  for(int iiTau = iTau0 ; iiTau<nTau ; iiTau++){
    Vindex_tau_time1[iiTau] = sum(time1<=tau[iiTau])-1; 
  }
  
  for(int iObs=0; iObs<nObs ; iObs++){ //  first event and after
    
    index_newT_time1 = sum(time1<=newT[iObs])-1;
    
    for(int iiTau = iTau0 ; iiTau<nTau ; iiTau++){
      
      IClambda0(iObs,iiTau) = 0;
      if(p>0){
        for(int iX=0; iX<p; iX++){
          IClambda0(iObs,iiTau) -= ICbeta(iObs,iX) * Elambda0(iX,iiTau);
        }
       }
      
      if(strata == newStrata[iObs]){
        IClambda0(iObs,iiTau) -= neweXb[iObs] * lamba0_iS0[min(Vindex_tau_time1[iiTau],index_newT_time1)];
        if(newT[iObs]<=tau[iiTau]){
          IClambda0(iObs,iiTau) += delta_iS0[iObs];
        }
      }
    }
  }
  
  return(IClambda0);
  
}  
