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
List IClambda0_cpp(const NumericVector& tau, const arma::mat& ICbeta,
                   const NumericVector& newT, const NumericVector& neweXb, const NumericVector& newStatus, const IntegerVector& newStrata, const IntegerVector& newIndexJump, 
                   const NumericVector& S01, const arma::mat& E1, const NumericVector& time1, double lastTime1, const NumericVector& lambda0,
                   int p, int strata){
  
  int nObs = newT.size();
  int nTau = tau.size();
  int nTime1 = time1.size();
  arma::mat IClambda0(nObs, nTau);
  arma::mat ICLambda0(nObs, nTau);
  IClambda0.fill(NA_REAL);
  ICLambda0.fill(NA_REAL);

  // Compute delta_iS0
  NumericVector delta_iS0(nObs);
  for(int iObs=0; iObs<nObs ; iObs++){
    delta_iS0[iObs] = newStatus[iObs]/S01[newIndexJump[iObs]];
  } 
  
  // Compute  Elambda0 and cumLamba0_iS0
  arma::mat Elambda0(p, nTau), cumElambda0(p, nTau);
  Elambda0.fill(0.0); cumElambda0.fill(0.0);
  NumericVector lambda0_iS0(nTime1,0.0), cumLambda0_iS0(nTime1,0.0);
  int iTau = 0;
  
  colvec Elambda0_iter(p), cumElambda0_iter(p); 
  for(int iTime1 = 0; iTime1 < nTime1; iTime1++){
    Elambda0_iter.fill(0);
    
    if(p>0){
      for(int iX = 0; iX < p; iX++){
        Elambda0_iter[iX] = E1(iTime1,iX) * lambda0[iTime1];
        cumElambda0_iter[iX] += Elambda0_iter[iX];
      }
    }
    
    // update (cum)Lambda0_iS0
    lambda0_iS0[iTime1] = lambda0[iTime1]/S01[iTime1];
    if(iTime1 == 0){
      cumLambda0_iS0[0] = lambda0_iS0[0];
    }else{
      cumLambda0_iS0[iTime1] = cumLambda0_iS0[iTime1-1] + lambda0_iS0[iTime1];  
    }
    
    
    // store in Elambda0 when time is tau if
    // (i) there are remaining tau
    // (ii) the next event time is strictly after tau OR it is the last event time and tau is before the last observation time
    if(p>0){
      while(iTau < nTau && ((iTime1 < (nTime1-1) && time1[iTime1+1]>tau[iTau]) || (iTime1==(nTime1-1) && tau[iTau]<= lastTime1) ) ){
        if(time1[iTime1]==tau[iTau]){Elambda0.col(iTau) = Elambda0_iter;}
        cumElambda0.col(iTau) = cumElambda0_iter;
        iTau++;
      }
    }
    
    if(iTau == nTau){ break; }
  }
  
  // main loop
  int iTau0 = 0;
  
  while((iTau0 < nTau) && time1[0]>tau[iTau0]){ // before the first event
    
    IClambda0.col(iTau0).zeros();
    ICLambda0.col(iTau0).zeros();
    iTau0++;
    
  }

  int index_newT_time1; // position of the minimum between t and t_train in cumLamba0_iS0
  int nTau_beforeLast=iTau0; // number of evaluation time before the last event
  IntegerVector Vindex_tau_time1(nTau);
  for(int iiTau = iTau0 ; iiTau<nTau ; iiTau++){
    Vindex_tau_time1[iiTau] = sum(time1<=tau[iiTau])-1; 
    if(tau[iiTau]<=lastTime1){nTau_beforeLast++;}
  }
  
  for(int iObs=0; iObs<nObs ; iObs++){ //  first event and after
    
    index_newT_time1 = sum(time1<=newT[iObs])-1;
    
    for(int iiTau = iTau0 ; iiTau<nTau_beforeLast ; iiTau++){
      
      IClambda0(iObs,iiTau) = 0;
      ICLambda0(iObs,iiTau) = 0;
      
      // first term
      if(p>0){
        for(int iX=0; iX<p; iX++){
          if(tau[iiTau]==time1[Vindex_tau_time1[iiTau]]){IClambda0(iObs,iiTau) -= ICbeta(iObs,iX) * Elambda0(iX,iiTau);}
          ICLambda0(iObs,iiTau) -= ICbeta(iObs,iX) * cumElambda0(iX,iiTau);
        }
      }
      
      if(strata == newStrata[iObs]){
        // // second term
        if(tau[iiTau]==time1[Vindex_tau_time1[iiTau]] && time1[Vindex_tau_time1[iiTau]] <= newT[iObs]){ IClambda0(iObs,iiTau) -= neweXb[iObs] * lambda0_iS0[Vindex_tau_time1[iiTau]]; }
        ICLambda0(iObs,iiTau) -= neweXb[iObs] * cumLambda0_iS0[min(Vindex_tau_time1[iiTau],index_newT_time1)];

        // third term
        if(newT[iObs]<=tau[iiTau]){
          if(newT[iObs]==tau[iiTau]){IClambda0(iObs,iiTau) += delta_iS0[iObs];}
          ICLambda0(iObs,iiTau) += delta_iS0[iObs];
        }
      }
    }
  }
  
  // export
  // return(List::create(Named("cumLambda0_iS0") = cumLambda0_iS0, 
  //                     Named("cumElambda0") = cumElambda0, 
  //                     Named("Elambda0") = Elambda0, 
  //                     Named("hazard") = IClambda0,
  //                     Named("cumhazard") = ICLambda0));
  
  return(List::create(Named("hazard") = IClambda0,
                      Named("cumhazard") = ICLambda0));
  
}  


