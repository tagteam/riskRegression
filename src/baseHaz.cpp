// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "sortS.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

struct structExport {
  vector<double> time;
  vector<double> hazard;
  vector<double> cumHazard;
  vector<double> SEhazard;
  vector<double> SEcumHazard;
  arma::mat Xbar;
  arma::mat XbarCumSum;
  int n;
};

structExport baseHazStrata_cpp(const vector<double>& alltimes, const vector<int>& status, const vector<double>& eXb, 
                               bool se, const arma::mat& data, int nVar,
                               int nPatients, double maxtime, int cause, bool Efron);

structExport subset_structExport(const structExport& resAll, const vector<double>& newtimes,
                                 double emaxtimes, int nNew, int nVar, bool se);

//// Data used to fit the Cox model
// alltimes: event times
// status:
// eXb: 
// strata:
// data
// emaxtimes: last event time in each strata
// nPatients
// nStrata
//// For prediction
// times_pred times at which we want to compute the hazard/survival
// [[Rcpp::export]]
List baseHaz_cpp(const NumericVector& alltimes, const IntegerVector& status, const NumericVector& eXb, const IntegerVector& strata,
                 bool se, arma::mat data, int nVar,
                 const std::vector<double>& predtimes, const NumericVector& emaxtimes,
                 int nPatients, int nStrata, int cause, bool Efron){
  // WARNING strata0 must begin at 0
  // WARNING predtimes must be sorted 
  
  vector<int> nObsStrata(nStrata,0);
  vector< vector<double> > alltimes_S(nStrata);
  vector< vector<int> > status_S(nStrata);
  vector< vector<double> > eXb_S(nStrata);
  vector< uvec > index_S(nStrata);
  vector< arma::mat > data_S(nStrata);
  uvec seqVar;
  if(se){
    seqVar = linspace<uvec>(0,nVar-1,nVar);
  }
  int nPredtimes = predtimes.size();
  double maxtime, max_predtimes = 0; // factice intialisation to avoid warning from the compiler
  if(nPredtimes>0){
    max_predtimes = predtimes[nPredtimes-1];
  }
  
  ////// 1- Strata  
  if(nStrata == 1){
    
    nObsStrata[0] = nPatients;
    alltimes_S[0].resize(nPatients);
    status_S[0].resize(nPatients);
    eXb_S[0].resize(nPatients);
    index_S[0].set_size(nPatients);
    if(se){data_S[0].set_size(nPatients,nVar);}
    
    for(int iter_p = 0 ; iter_p < nPatients ; iter_p++){
      alltimes_S[0][iter_p] = alltimes[iter_p];
      status_S[0][iter_p] = status[iter_p];
      eXb_S[0][iter_p] = eXb[iter_p];
      index_S[0][iter_p] = iter_p;
      if(se){data_S[0].row(iter_p) = data.row(iter_p);}
    }
    
  }else{
    
    //// a- Patient per strata
    for(int iter_p = 0 ; iter_p < nPatients ; iter_p++){
      nObsStrata[strata[iter_p]]++;
    }
    
    //// b- Subset the dataset per strata
    for(int iter_s = 0 ; iter_s < nStrata ; iter_s++){
      alltimes_S[iter_s].resize(nObsStrata[iter_s]);
      status_S[iter_s].resize(nObsStrata[iter_s]);
      eXb_S[iter_s].resize(nObsStrata[iter_s]);
      index_S[iter_s].set_size(nObsStrata[iter_s]);
      if(se){data_S[iter_s].set_size(nObsStrata[iter_s],nVar);}
    }
    
    vector<int> index_tempo(nStrata,0); // indicates position of the last observation in each strata
    int strata_tempo;
    
    for(int iter_p = 0 ; iter_p < nPatients ; iter_p++){
      strata_tempo = strata[iter_p];
      alltimes_S[strata_tempo][index_tempo[strata_tempo]] = alltimes[iter_p];
      status_S[strata_tempo][index_tempo[strata_tempo]] = status[iter_p];
      eXb_S[strata_tempo][index_tempo[strata_tempo]] = eXb[iter_p];
      index_S[strata_tempo][index_tempo[strata_tempo]] = index_tempo[strata_tempo];
      if(se){data_S[strata_tempo].row(index_tempo[strata_tempo]) = data.row(iter_p);}
      index_tempo[strata_tempo]++;
    }
  }
  
  ////// 2- Select and compute
  structExport resH;
  vector<double> timeRes(0);
  vector<double> hazardRes(0);
  vector<double> SEhazardRes(0);
  vector<double> cumHazardRes(0);
  vector<double> SEcumHazardRes(0);
  vector<double> strataRes(0);
  arma::mat XbarRes(0,nVar);
  arma::mat XbarCumSumRes(0,nVar);
  
  for(int iter_s = 0 ; iter_s < nStrata ; iter_s++){
    R_CheckUserInterrupt();
    
    // reorder the data
    sortS(alltimes_S[iter_s], status_S[iter_s], eXb_S[iter_s], index_S[iter_s], nObsStrata[iter_s]); // update alltimes, status and eXb
    
    if(se){
      data_S[iter_s] = data_S[iter_s].submat(index_S[iter_s], seqVar);
    }
    
    if(nPredtimes>0){ // set maxtime to the first event after the maximum prediction time
      int i = 0;
      while(i<(nObsStrata[iter_s]-1) && alltimes_S[iter_s][i]<max_predtimes) i++;
      maxtime = alltimes_S[iter_s][i];
    }else{
      maxtime = emaxtimes[iter_s];
    }
    
    // compute the hazard
    resH = baseHazStrata_cpp(alltimes_S[iter_s], status_S[iter_s], eXb_S[iter_s],
                             se, data_S[iter_s], nVar, 
                             nObsStrata[iter_s], maxtime, cause, 
                             Efron);
    
    // subset results according to predtime
    if(nPredtimes>0){
      resH = subset_structExport(resH, predtimes,
                                 emaxtimes[iter_s], nPredtimes, nVar, se);
    }
    
    // store results
    timeRes.insert( timeRes.end(), resH.time.begin(), resH.time.end() );
    hazardRes.insert( hazardRes.end(), resH.hazard.begin(), resH.hazard.end() );
    cumHazardRes.insert( cumHazardRes.end(), resH.cumHazard.begin(), resH.cumHazard.end() );
    if(nStrata > 1){
      strataRes.resize( strataRes.size() + resH.n, iter_s);
    }
    if(se){
      SEhazardRes.insert( SEhazardRes.end(), resH.SEhazard.begin(), resH.SEhazard.end() );
      SEcumHazardRes.insert( SEcumHazardRes.end(), resH.SEcumHazard.begin(), resH.SEcumHazard.end() );
      XbarRes = join_vert(XbarRes, resH.Xbar);
      XbarCumSumRes = join_vert(XbarCumSumRes, resH.XbarCumSum);
    }
  }
  
  
  //// 3- export
  List res = List::create(Named("time")  = timeRes,
                          Named("hazard")  = hazardRes,
                          Named("se.hazard")  = SEhazardRes,
                          Named("cumHazard")  = cumHazardRes,
                          Named("se.cumHazard")  = SEcumHazardRes,
                          Named("strata")  = strataRes,
                          Named("Xbar")  = XbarRes,
                          Named("XbarCumSumRes")  = XbarCumSumRes
  );
  
  return(res);  
  
}


structExport baseHazStrata_cpp(const vector<double>& alltimes, const vector<int>& status, const vector<double>& eXb, 
                               bool se, const arma::mat& data, int nVar,
                               int nPatients, double maxtime, int cause, bool Efron){
  
  //// 1- count the number of events
  size_t nEvents = 1, nEventsLast = 1;
  
  for(int iterPat = 1 ; iterPat < nPatients ; iterPat++){
    
    if(alltimes[iterPat] != alltimes[iterPat-1]){
      nEvents++; // total number of events
      if(maxtime >= alltimes[iterPat]){
        nEventsLast++; // up to maxtime
      }
    }
    
  }
  
  ////// 2- merge the data by event
  // sumEXb : sum(patient still at risk) exp(XB)
  // sumEXb_event : sum(patient having the event at that time) exp(XB)
  // sumEXb_data : sum(patient still at risk) X exp(XB)
  // sumEXb_eventData : sum(patient having the event at that time) X exp(XB)
  
  // temp2 <- rowsum(risk * x, dtime) # at each time the sum of the product between the individual risk and the value of the covariates (E[Xexp(Xb)])
  // xsum <- apply(temp2, 2, rcumsum) # cumulative E[Xexp(XB)]
  // xsum2 <- rowsum((risk * death) * x, dtime) # same as temp2 but the sum is only for people with event 
  
  vector<double> time(nEventsLast);
  vector<double> hazard(nEventsLast, NA_REAL), SEhazard, SEcumHazard;
  vector<double> cumHazard(nEventsLast, NA_REAL);
  arma::mat Xbar, XbarCumSum;
  if(se){
    SEhazard.resize(nEventsLast, NA_REAL);  
    SEcumHazard.resize(nEventsLast, NA_REAL);
    Xbar.set_size(nEvents, nVar); Xbar.fill(0);
  }
  
  vector<int> death(nEvents,0);
  vector<double> sumEXb(nEvents), sumEXb2;
  vector<double> sumEXb_event(nEvents,0.0); // only used if efron correction
  arma::mat sumEXb_data, sumEXb_eventData;
  if(se){
    sumEXb2.resize(nEvents, NA_REAL);
    sumEXb_data.set_size(nEvents, nVar);
    sumEXb_eventData.set_size(nEvents, nVar); sumEXb_eventData.fill(0);
  }
  
  // last observation
  sumEXb[nEvents-1] = eXb[nPatients-1];
  death[nEvents-1] = (status[nPatients-1]==cause);
  if(Efron && (status[nPatients-1]==cause)){
    sumEXb_event[nEvents-1] += eXb[nPatients-1];
    if(se){sumEXb_eventData.row(nEvents-1) = data.row(nPatients-1) * eXb[nPatients-1];}
  }
  if(se){sumEXb_data.row(nEvents-1) = data.row(nPatients-1) * eXb[nPatients-1];}
  
  // remaining observations
  size_t index_tempo = nEvents-1;
  for(int iterPat = nPatients-2 ; iterPat >= 0 ; iterPat--){
    if(alltimes[iterPat] != alltimes[iterPat+1]){
      index_tempo--;
      sumEXb[index_tempo] = sumEXb[index_tempo+1];
      if(se){sumEXb_data.row(index_tempo) = sumEXb_data.row(index_tempo+1);}
    }
    
    death[index_tempo] += (status[iterPat]==cause);
    sumEXb[index_tempo] += eXb[iterPat];
    if(Efron && (status[iterPat]==cause)){
      sumEXb_event[index_tempo] += eXb[iterPat];
      if(se){sumEXb_eventData.row(index_tempo) += data.row(iterPat) * eXb[iterPat];}
    }
    if(se){sumEXb_data.row(index_tempo) += data.row(iterPat) * eXb[iterPat];}
  }
  
  // first observation
  time[0] = alltimes[0];
  
  // remaining observations
  index_tempo = 0;
  for(int iterPat = 1 ; iterPat < nPatients ; iterPat++){
    if(alltimes[iterPat] != alltimes[iterPat-1]){
      if(maxtime < alltimes[iterPat]){break;}      // if after the last time we want to predict
      index_tempo++;
      time[index_tempo] = alltimes[iterPat];
    }
  }
  
  if(se){
    for(size_t iterEvent = 0 ; iterEvent < nEventsLast ; iterEvent++){
      if(death[iterEvent]>0){ // otherwise it will not be used anyway to compute SEhazard
        sumEXb2[iterEvent] = pow(sumEXb[iterEvent],2);
        Xbar.row(iterEvent) = sumEXb_data.row(iterEvent) * pow(1/sumEXb[iterEvent],2);
      }
    }
  }    
  
  //// OPT- Efron correction [from the survival package, function agsurv5]
  if(Efron){
    
    double Wm1_tempo, Wm2_tempo = NA_REAL, sumRi, sumRi_di, sumRi_Efron, di; // it is important that di is a double and not an in for the division in the for loop
    //sumRi is the sum over the patient at risk of exp(Xbeta)
    //sumRi_di is the sum of the number at risk who experience the event at the specific time, of exp(Xbeta)
    //sumRi_Efron is the corrected sum of the number at risk of exp(Xbeta)
    
    for(size_t iterEvent = 0 ; iterEvent < nEventsLast ; iterEvent++){
      
      if (death[iterEvent]>1){
        sumRi = sumEXb[iterEvent];
        sumRi_di = sumEXb_event[iterEvent];
        di = death[iterEvent];
        
        Wm1_tempo = 1/sumEXb[iterEvent];
        if(se){ Wm2_tempo = 1/sumEXb2[iterEvent];}
        
        for(int iterPat = 1; iterPat < di; iterPat++){
          sumRi_Efron = 1/(sumRi - (iterPat/di)*sumRi_di);
          Wm1_tempo += sumRi_Efron;
          if(se){
            Wm2_tempo += pow(sumRi_Efron,2);
            Xbar.row(iterEvent) += (sumEXb_data.row(iterEvent) - (iterPat/di)*sumEXb_eventData.row(iterEvent)) * pow(sumRi_Efron,2);
          }
          
        }
        
        // Make the average over the patient having the event at time i
        sumEXb[iterEvent] = di/Wm1_tempo;
        if(se){sumEXb2[iterEvent] = di/Wm2_tempo;}
        
      }
    }
  }else if(se){ // Breslow se
    for(size_t iterEvent = 0 ; iterEvent < nEventsLast ; iterEvent++){
      Xbar.row(iterEvent) *= death[iterEvent];
    }
  }
  
  if(se){XbarCumSum =  cumsum(Xbar, 0);}  // cumulative sum by column
  
  //// 3- Computation of the hazards
  hazard[0] = death[0] / sumEXb[0];
  cumHazard[0] = hazard[0];
  if(se){
    if(death[0]>0){
      SEhazard[0] = death[0] / sumEXb2[0];
      SEcumHazard[0] = SEhazard[0];
    }else{
      SEhazard[0] = 0;
      SEcumHazard[0] = 0;
    }
  }
  
  for( size_t iterTime = 1 ; iterTime < nEventsLast ; iterTime ++){ // up to last time
    hazard[iterTime] = death[iterTime] / sumEXb[iterTime];
    cumHazard[iterTime] = cumHazard[iterTime-1] + hazard[iterTime];
    
    if(se){
      if(death[iterTime]>0){
        SEhazard[iterTime] = death[iterTime] / sumEXb2[iterTime];
        SEcumHazard[iterTime] = SEcumHazard[iterTime-1] + SEhazard[iterTime];
      }else{
        SEhazard[iterTime] = 0;
        SEcumHazard[iterTime] = SEcumHazard[iterTime-1];
      }
    }
    
  }
  
  //// export
  structExport res;
  res.time = time;
  res.hazard = hazard;
  res.cumHazard = cumHazard;
  res.SEhazard = SEhazard;
  res.SEcumHazard = SEcumHazard;
  res.Xbar = Xbar;
  res.XbarCumSum = XbarCumSum;
  res.n = nEventsLast;
  
  return res;
}

structExport subset_structExport(const structExport& resAll, const vector<double>& newtimes, 
                                 double emaxtimes,
                                 int nNew, int nVar, bool se){
  
  structExport resSubset;
  resSubset.time = newtimes;
  resSubset.hazard.resize(nNew, NA_REAL);
  resSubset.cumHazard.resize(nNew, NA_REAL);
  if(se){
    resSubset.SEhazard.resize(nNew, NA_REAL);
    resSubset.SEcumHazard.resize(nNew, NA_REAL);
    resSubset.Xbar.set_size(nNew, nVar); resSubset.Xbar.fill(NA_REAL);
    resSubset.XbarCumSum.set_size(nNew, nVar); resSubset.XbarCumSum.fill(NA_REAL);
  }
  
  int i = 0;
  
  for (int t=0;t<nNew;t++){
    
    // update index
    while(i<(resAll.n-1) && resAll.time[i+1]<=newtimes[t]){i++;}
    
    // update hazard
    if(newtimes[t]<=emaxtimes){
      
      if(resAll.time[i]==newtimes[t]){ 
        resSubset.hazard[t] = resAll.hazard[i];
        if(se){resSubset.SEhazard[t] = resAll.SEhazard[i];}
      }else{ // if not an event time then 0
        resSubset.hazard[t] = 0;
        if(se){resSubset.SEhazard[t] = 0;} // not clear which value should be put here
      }
      
      if(newtimes[t]>=resAll.time[0]){
        resSubset.cumHazard[t] = resAll.cumHazard[i];
        if(se){resSubset.SEcumHazard[t] = resAll.SEcumHazard[i];}
      }else{  // if before the first event then 0
        resSubset.cumHazard[t] = 0;
        if(se){resSubset.SEcumHazard[t] = 0;} // not clear which value should be put here
      }
      
    }
    
    if(se){
      resSubset.Xbar.row(t) = resAll.Xbar.row(i);
      resSubset.XbarCumSum.row(t) = resAll.XbarCumSum.row(i);
    }
  }
  resSubset.n = nNew;
  
  // export
  return(resSubset);
}

//// NOTE Adaptation of the Cagsurv5 from the survival package with no weight
// ntimes number of observations (unique death times)
// ndead number of deaths at that time
// risk number at risk at the time
// riskDead number at risk at the time for the deaths

// [[Rcpp::export]]
NumericVector baseHazEfron_survival_cpp(int ntimes, IntegerVector ndead, 
                                        NumericVector risk, NumericVector riskDead) {
  double temp;
  int t, j;
  double di;
  NumericVector W(ntimes, 0.0);
  
  for (t=0; t < ntimes; t++) {
    di = ndead[t];
    
    if (di==1){
      
      W[t] = 1/risk[t];
      
    } else if (di>1){
      
      for (j=0; j < di; j++) {
        
        temp = 1/(risk[t] - riskDead[t]*j/di);
        W[t] += temp/di;
        
      }
      
    }
    
  }
  
  return(W);
}
