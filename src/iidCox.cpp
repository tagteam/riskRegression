// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// * calcE_cpp
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

// * IFbeta_cpp
// [[Rcpp::export]]
arma::mat IFbeta_cpp(const NumericVector& newT, const NumericVector& neweXb, const arma::mat& newX, const NumericVector& newStatus, const IntegerVector& newIndexJump, 
                     const NumericVector& S01, const arma::mat& E1, const NumericVector& time1, const arma::mat& iInfo,
                     int p){
  
  arma::mat IFbeta;
  int nObs = newIndexJump.size();
  
  if(p==0){
    IFbeta.resize(nObs, 1);
    IFbeta.fill(0);
  }else{
    IFbeta.resize(nObs, p);
    IFbeta.fill(NA_REAL);
    
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
      
      IFbeta.row(iObs) = (iInfo * Score).t();
      
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
        
        IFbeta.row(iObs) = (iInfo * Score).t();
        iObs++;
      }
      
      if(iObs == nObs){ break; }
      
    }
    
    
  }
  
  return(IFbeta);
}

// * IFlambda0_cpp
// [[Rcpp::export]]
List IFlambda0_cpp(const NumericVector& tau, const arma::mat& IFbeta,
                   const NumericVector& newT, const NumericVector& neweXb, const NumericVector& newStatus, const IntegerVector& newStrata, const IntegerVector& newIndexJump, 
                   const NumericVector& S01, const arma::mat& E1, const NumericVector& time1, double lastTime1, const NumericVector& lambda0,
                   int p, int strata, bool minimalExport){
  
  int nObs = newT.size();
  int nTau = tau.size();
  int nTime1 = time1.size();

  // ** Prepare output
  arma::mat IFlambda0(nObs, std::max(nTau,1));
  arma::mat IFLambda0(nObs, std::max(nTau,1));
  IFlambda0.fill(NA_REAL);
  IFLambda0.fill(NA_REAL);

  arma::mat Elambda0(p, std::max(nTau,1), fill::zeros);
  arma::mat cumElambda0(p, std::max(nTau,1), fill::zeros);
  NumericVector lambda0_iS0(nTime1,0.0);
  NumericVector cumLambda0_iS0(nTime1,0.0);
  
  // ** Find early prediction times
  // if iTau0 = 5 this means that the first five prediction times are before the first event, i.e IF = 0
  int iTau0 = 0;
  while(iTau0 < nTau && time1[0]>tau[iTau0]){
	iTau0++;
  }

  // ** Compute delta_iS0
  NumericVector delta_iS0(nObs,0.0);
  for(int iObs=0; iObs<nObs ; iObs++){
	if(strata == newStrata[iObs]){
	  delta_iS0[iObs] = newStatus[iObs]/S01[newIndexJump[iObs]];
	}
  } 

  // Exclude case with no Tau
  if(nTau==0){
	if(minimalExport){
	  NumericVector indeXb = clone(neweXb);
	  indeXb[newStrata == strata] = 0;
	  return(List::create(Named("Elambda0") = Elambda0,
						  Named("cumElambda0") = cumElambda0,
						  Named("eXb") = indeXb,
						  Named("lambda0_iS0") = lambda0_iS0,
						  Named("cumLambda0_iS0") = cumLambda0_iS0,
						  Named("delta_iS0") = delta_iS0,
						  Named("time1") = time1
						  ));
	}else{
	  IFlambda0.fill(0.0);
	  IFLambda0.fill(0.0);
	  return(List::create(Named("hazard") = IFlambda0,
						  Named("cumhazard") = IFLambda0));
	}
  }

  // Compute  Elambda0 and cumLamba0_iS0
  int iTau = iTau0;
  
  colvec Elambda0_iter(p), cumElambda0_iter(p); 
  cumElambda0_iter.fill(0);
  
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
  
  if(minimalExport){
	NumericVector indeXb = clone(neweXb);
	indeXb[newStrata != strata] = 0;
	
    return(List::create(Named("delta_iS0") = delta_iS0,
                        Named("Elambda0") = Elambda0,
                        Named("cumElambda0") = cumElambda0,
						Named("eXb") = indeXb,
						Named("lambda0_iS0") = lambda0_iS0,
                        Named("cumLambda0_iS0") = cumLambda0_iS0,
                        Named("time1") = time1));
  }
  
  // main loop
  if(iTau0>0){
    for(int iiTau=0 ; iiTau<iTau0; iiTau++){ // before the first event    
	  IFlambda0.col(iiTau).zeros();
	  IFLambda0.col(iiTau).zeros();    
	}
  }
  
  int index_newT_time1; // position of the minimum between t and t_train in cumLamba0_iS0
  int nTau_beforeLast=iTau0; // number of evaluation time before the last event
  IntegerVector Vindex_tau_time1(nTau);
  for(int iiTau = iTau0 ; iiTau<nTau ; iiTau++){
    Vindex_tau_time1[iiTau] = sum(time1<=tau[iiTau])-1; 
    if(tau[iiTau]<=lastTime1){nTau_beforeLast++;}
  }
  
  for(int iObs=0; iObs<nObs ; iObs++){ //  first event and after
    
    // index_newT_time1 can value -1 when newT is before the first event 
    // in this case cumLambda0 is 0 so the second term can be skipped
    index_newT_time1 = sum(time1<=newT[iObs])-1; 
    
    for(int iiTau = iTau0 ; iiTau<nTau_beforeLast ; iiTau++){
      
      IFlambda0(iObs,iiTau) = 0;
      IFLambda0(iObs,iiTau) = 0;
      
      // first term
      if(p>0){
        for(int iX=0; iX<p; iX++){
          if(tau[iiTau]==time1[Vindex_tau_time1[iiTau]]){IFlambda0(iObs,iiTau) -= IFbeta(iObs,iX) * Elambda0(iX,iiTau);}
          IFLambda0(iObs,iiTau) -= IFbeta(iObs,iX) * cumElambda0(iX,iiTau);
        }
      }
      
      if(strata == newStrata[iObs]){
        // second term
		if(tau[iiTau]==time1[Vindex_tau_time1[iiTau]] && time1[Vindex_tau_time1[iiTau]] <= newT[iObs]){
		  IFlambda0(iObs,iiTau) -= neweXb[iObs] * lambda0_iS0[Vindex_tau_time1[iiTau]];
		}
		IFLambda0(iObs,iiTau) -= neweXb[iObs] * cumLambda0_iS0[min(Vindex_tau_time1[iiTau],index_newT_time1)];

        // third term
        if(newT[iObs]<=tau[iiTau]){
          if(newT[iObs]==tau[iiTau]){
          IFlambda0(iObs,iiTau) += delta_iS0[iObs];}
          IFLambda0(iObs,iiTau) += delta_iS0[iObs];
        }
      }
    }   
  }
  
  return(List::create(Named("hazard") = IFlambda0,
                      Named("cumhazard") = IFLambda0));
  
}  

// }}}
