// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

struct structExport {
  vector<double> time;
  vector<double> hazard;
  vector<double> cumhazard;
  int n;
};

structExport baseHazStrata_cpp(const vector<double>& alltimes,
                               const vector<int>& status,
                               const vector<double>& eXb, 
                               int nPatients,
                               double maxtime,
                               int cause,
                               bool Efron);

structExport subset_structExport(const structExport& resAll,
                                 const vector<double>& newtimes,
                                 double emaxtimes,
                                 int nNew);

//// Data used to fit the Cox model
// alltimes: event times. 
// status: status 0/1
// eXb: linear predictor 
// strata: must begin at 0
// data
// emaxtimes: last event time in each strata
// nPatients
// nStrata
//// For prediction
// times_pred times at which we want to compute the hazard/survival. Must be sorted 
//// WARNING alltimes status eXb and strata must be sorted by strata, alltimes, and status
// [[Rcpp::export]]
List baseHaz_cpp(const NumericVector& alltimes,
                 const IntegerVector& status,
                 const NumericVector& eXb,
                 const IntegerVector& strata,
                 const std::vector<double>& predtimes,
                 const NumericVector& emaxtimes,
                 int nPatients,
                 int nStrata,
                 int cause,
                 bool Efron){
  
  vector<int> nObsStrata(nStrata,0);
  vector< vector<double> > alltimes_S(nStrata);
  vector< vector<int> > status_S(nStrata);
  vector< vector<double> > eXb_S(nStrata);
  vector< uvec > index_S(nStrata);
  uvec seqVar;
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
    
    for(int iter_p = 0 ; iter_p < nPatients ; iter_p++){
      alltimes_S[0][iter_p] = alltimes[iter_p];
      status_S[0][iter_p] = status[iter_p];
      eXb_S[0][iter_p] = eXb[iter_p];
      index_S[0][iter_p] = iter_p;
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
    }
    
    vector<int> index_tempo(nStrata,0); // indicates position of the last observation in each strata
    int strata_tempo;
    
    for(int iter_p = 0 ; iter_p < nPatients ; iter_p++){
      strata_tempo = strata[iter_p];
      alltimes_S[strata_tempo][index_tempo[strata_tempo]] = alltimes[iter_p];
      status_S[strata_tempo][index_tempo[strata_tempo]] = status[iter_p];
      eXb_S[strata_tempo][index_tempo[strata_tempo]] = eXb[iter_p];
      index_S[strata_tempo][index_tempo[strata_tempo]] = index_tempo[strata_tempo];
      index_tempo[strata_tempo]++;
    }
  }
  
  ////// 2- Select and compute
  structExport resH;
  vector<double> timeRes(0);
  vector<double> hazardRes(0);
  vector<double> cumhazardRes(0);
  vector<double> strataRes(0);
  
  for(int iter_s = 0 ; iter_s < nStrata ; iter_s++){
    R_CheckUserInterrupt();
    
    if(nPredtimes>0){ // set maxtime to the first event after the maximum prediction time
      int i = 0;
      while(i<(nObsStrata[iter_s]-1) && alltimes_S[iter_s][i]<max_predtimes) i++;
      maxtime = alltimes_S[iter_s][i];
    }else{
      maxtime = emaxtimes[iter_s];
    }
    
    // compute the hazard
    resH = baseHazStrata_cpp(alltimes_S[iter_s], status_S[iter_s], eXb_S[iter_s],
                             nObsStrata[iter_s], maxtime, cause, 
                             Efron);
    
    // subset results according to predtime
    if(nPredtimes>0){
      resH = subset_structExport(resH,
                                 predtimes,
                                 emaxtimes[iter_s],
                                          nPredtimes);
    }
    
    // store results
    timeRes.insert( timeRes.end(), resH.time.begin(), resH.time.end() );
    hazardRes.insert( hazardRes.end(), resH.hazard.begin(), resH.hazard.end() );
    cumhazardRes.insert( cumhazardRes.end(), resH.cumhazard.begin(), resH.cumhazard.end() );
    if(nStrata > 1){
      strataRes.resize( strataRes.size() + resH.n, iter_s);
    }
  }
  
  
  //// 3- export
  List res = List::create(Named("time")  = timeRes,
                          Named("hazard")  = hazardRes,
                          Named("cumhazard")  = cumhazardRes,
                          Named("strata")  = strataRes);
  
  return(res);  
}


structExport baseHazStrata_cpp(const vector<double>& alltimes,
                               const vector<int>& status,
                               const vector<double>& eXb, 
                               int nPatients,
                               double maxtime,
                               int cause,
                               bool Efron){
  
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
  // temp2 <- rowsum(risk * x, dtime) # at each time the sum of the product between the individual risk and the value of the covariates (E[Xexp(Xb)])
  // xsum <- apply(temp2, 2, rcumsum) # cumulative E[Xexp(XB)]
  // xsum2 <- rowsum((risk * death) * x, dtime) # same as temp2 but the sum is only for people with event 
  
  vector<double> time(nEventsLast);
  vector<double> hazard(nEventsLast, NA_REAL);
  vector<double> cumhazard(nEventsLast, NA_REAL);
  
  vector<int> death(nEvents,0);
  vector<double> sumEXb(nEvents), sumEXb2;
  vector<double> sumEXb_event(nEvents,0.0); // only used if efron correction
  
  // last observation
  sumEXb[nEvents-1] = eXb[nPatients-1];
  death[nEvents-1] = (status[nPatients-1]==cause);
  if(Efron && (status[nPatients-1]==cause)){
    sumEXb_event[nEvents-1] += eXb[nPatients-1];
  }
  
  // remaining observations
  size_t index_tempo = nEvents-1;
  for(int iterPat = nPatients-2 ; iterPat >= 0 ; iterPat--){
    if(alltimes[iterPat] != alltimes[iterPat+1]){
      index_tempo--;
      sumEXb[index_tempo] = sumEXb[index_tempo+1];
    }
    
    death[index_tempo] += (status[iterPat]==cause);
    sumEXb[index_tempo] += eXb[iterPat];
    if(Efron && (status[iterPat]==cause)){
      sumEXb_event[index_tempo] += eXb[iterPat];
    }
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
  
  //// OPT- Efron correction [from the survival package, function agsurv5]
  if(Efron){
    
    double Wm1_tempo, sumRi, sumRi_di, sumRi_Efron, di; // it is important that di is a double and not an in for the division in the for loop
    //sumRi is the sum over the patient at risk of exp(Xbeta)
    //sumRi_di is the sum of the number at risk who experience the event at the specific time, of exp(Xbeta)
    //sumRi_Efron is the corrected sum of the number at risk of exp(Xbeta)
    
    for(size_t iterEvent = 0 ; iterEvent < nEventsLast ; iterEvent++){
      
      if (death[iterEvent]>1){
        sumRi = sumEXb[iterEvent];
        sumRi_di = sumEXb_event[iterEvent];
        di = death[iterEvent];
        
        Wm1_tempo = 1/sumEXb[iterEvent];
        
        for(int iterPat = 1; iterPat < di; iterPat++){
          sumRi_Efron = 1/(sumRi - (iterPat/di)*sumRi_di);
          Wm1_tempo += sumRi_Efron;
        }
        
        // Make the average over the patient having the event at time i
        sumEXb[iterEvent] = di/Wm1_tempo;
        
      }
    }
  }
  
  //// 3- Computation of the hazards
  hazard[0] = death[0] / sumEXb[0];
  cumhazard[0] = hazard[0];
  
  for( size_t iterTime = 1 ; iterTime < nEventsLast ; iterTime ++){ // up to last time
    hazard[iterTime] = death[iterTime] / sumEXb[iterTime];
    cumhazard[iterTime] = cumhazard[iterTime-1] + hazard[iterTime];
  }
  
  //// export
  structExport res;
  res.time = time;
  res.hazard = hazard;
  res.cumhazard = cumhazard;
  res.n = nEventsLast;
  
  return res;
}

structExport subset_structExport(const structExport& resAll, const vector<double>& newtimes, 
                                 double emaxtimes,
                                 int nNew){
  
  structExport resSubset;
  resSubset.time = newtimes;
  resSubset.hazard.resize(nNew, NA_REAL);
  resSubset.cumhazard.resize(nNew, NA_REAL);
  int i = 0;
  
  for (int t=0;t<nNew;t++){
    
    // update index
    while(i<(resAll.n-1) && resAll.time[i+1]<=newtimes[t]){i++;}
    
    // update hazard
    if(newtimes[t]<=emaxtimes){
      
      if(resAll.time[i]==newtimes[t]){ 
        resSubset.hazard[t] = resAll.hazard[i];
      }else{ // if not an event time then 0
        resSubset.hazard[t] = 0;
      }
      
      if(newtimes[t]>=resAll.time[0]){
        resSubset.cumhazard[t] = resAll.cumhazard[i];
      }else{  // if before the first event then 0
        resSubset.cumhazard[t] = 0;
      }
      
    }
    
  }
  resSubset.n = nNew;
  
  // export
  return(resSubset);
}
