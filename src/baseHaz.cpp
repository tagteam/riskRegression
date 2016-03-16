#include <Rcpp.h>
#include "sortS.h"

using namespace Rcpp;
using namespace std;

const double CST_EPSILON = 1e-10;

vector< vector<double> > BaseHaz_cpp(const vector<double>& alltimes, const vector<int>& status, const vector<double>& Xb, 
                                     int nPatients, double lasttime, int cause, bool Efron, bool addFirst, bool addLast);

// [[Rcpp::export]]
List BaseHazStrata_cpp(const NumericVector& alltimes, const IntegerVector& status, const NumericVector& Xb, const IntegerVector& strata,
                       int nPatients, int nStrata, double lasttime, int cause, bool Efron, bool addFirst, bool addLast){
  // WARNING strata0 must begin at 0
  
  vector<int> nObsStrata(nStrata,0);
  vector< vector<double> > alltimes_S(nStrata);
  vector< vector<int> > status_S(nStrata);
  vector< vector<double> > Xb_S(nStrata);
  
  ////// 1- Strata  
  if(nStrata == 1){
    
    nObsStrata[0] = nPatients;
    alltimes_S[0].resize(nPatients);
    status_S[0].resize(nPatients);
    Xb_S[0].resize(nPatients);
    
    for(int iter_p = 0 ; iter_p < nPatients ; iter_p++){
      alltimes_S[0][iter_p] = alltimes[iter_p];
      status_S[0][iter_p] = status[iter_p];
      Xb_S[0][iter_p] = Xb[iter_p];
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
      Xb_S[iter_s].resize(nObsStrata[iter_s]);
    }
    
    vector<int> index_tempo(nStrata,0);
    int strata_tempo;
    
    for(int iter_p = 0 ; iter_p < nPatients ; iter_p++){
      strata_tempo = strata[iter_p];
      alltimes_S[strata_tempo][index_tempo[strata_tempo]] = alltimes[iter_p];
      status_S[strata_tempo][index_tempo[strata_tempo]] = status[iter_p];
      Xb_S[strata_tempo][index_tempo[strata_tempo]] = Xb[iter_p];
      index_tempo[strata_tempo]++;
    }
  }
  
  ////// 2- Select and compute
  bool test;
  vector<vector <double> > resH;
  vector<double> timeRes(0);
  vector<double> hazardRes(0);
  vector<double> cumHazardRes(0);
  vector<double> strataRes(0);
  
  for(int iter_s = 0 ; iter_s < nStrata ; iter_s++){
    
    // reorder the data
    test = sortS(alltimes_S[iter_s], status_S[iter_s], Xb_S[iter_s], nObsStrata[iter_s]); // update alltimes, status and Xb
    
    // deal with specific cases
    if(test == false){
      return(List::create(Named("time") = NA_REAL,
                          Named("strata") = NA_REAL,
                          Named("hazard") = NA_REAL,
                          Named("cumHazard") = NA_REAL,
                          Named("Pb") = "Could not allocate memory in sortS")
      );
    }
    if(lasttime < alltimes_S[iter_s][0]){
      if(addFirst){
        timeRes.push_back(CST_EPSILON);
        hazardRes.push_back(0);
        cumHazardRes.push_back(0);
        if(nStrata > 1){
          strataRes.push_back(iter_s);
        }
      }
      if(addLast){
        timeRes.push_back(alltimes_S[iter_s][alltimes_S[iter_s].size()-1]+CST_EPSILON);
        hazardRes.push_back(NA_REAL);
        cumHazardRes.push_back(NA_REAL);
        if(nStrata > 1){
        strataRes.push_back(iter_s);
        }
      }
      continue;
    }
    
    // compute the hazard
    resH = BaseHaz_cpp(alltimes_S[iter_s], status_S[iter_s], Xb_S[iter_s], 
                       nObsStrata[iter_s], lasttime, cause, 
                       Efron, addFirst, addLast);
    
    // store results
    timeRes.insert( timeRes.end(), resH[0].begin(), resH[0].end() );
    hazardRes.insert( hazardRes.end(), resH[1].begin(), resH[1].end() );
    cumHazardRes.insert( cumHazardRes.end(), resH[2].begin(), resH[2].end() );
    if(nStrata > 1){
      strataRes.resize( strataRes.size() + resH[0].size(), iter_s);
    }
  }
  
  
  //// 3- export
  return(List::create(Named("time")  = timeRes,
                      Named("strata") = strataRes,
                      Named("hazard")  = hazardRes,
                      Named("cumHazard")  = cumHazardRes)
  );
}


vector< vector<double> > BaseHaz_cpp(const vector<double>& alltimes, const vector<int>& status, const vector<double>& Xb, 
                                     int nPatients, double lasttime, int cause, bool Efron, bool addFirst, bool addLast){
  
  //// 1- count the number of events
  size_t nEvents = 1, nEventsLast = 1;
  
  for(int iterPat = 1 ; iterPat < nPatients ; iterPat++){
    
    if(alltimes[iterPat] != alltimes[iterPat-1]){
      nEvents++; // total number of events
      if(lasttime >= alltimes[iterPat]){
        nEventsLast++; // up to lasttime
      }
    }
    
  }
  
  ////// 2- merge the data by event
  vector<double> time(nEventsLast);
  vector<int> death(nEvents,0);
  vector<double> sumEXb(nEvents);
  vector<double> sumEXb_event(nEvents,0.0); // only used if efron correction
  
  // last observation
  sumEXb[nEvents-1] = exp(Xb[nPatients-1]);
  death[nEvents-1] = (status[nPatients-1]==cause);
  if(Efron && (status[nPatients-1]==cause)){
    sumEXb_event[nEvents-1] += exp(Xb[nPatients-1]);
  }
  
  // remaining observations
  size_t index_tempo = nEvents-1;
  for(int iterPat = nPatients-2 ; iterPat >= 0 ; iterPat--){
    if(alltimes[iterPat] != alltimes[iterPat+1]){
      index_tempo--;
      sumEXb[index_tempo] = sumEXb[index_tempo+1];
    }
    
    death[index_tempo] += (status[iterPat]==cause);
    sumEXb[index_tempo] += exp(Xb[iterPat]);
    if(Efron && (status[iterPat]==cause)){
      sumEXb_event[index_tempo] += exp(Xb[iterPat]);
    }
  }
  
  // first observation
  time[0] = alltimes[0];
  
  // remaining observations
  index_tempo = 0;
  for(int iterPat = 1 ; iterPat < nPatients ; iterPat++){
    if(alltimes[iterPat] != alltimes[iterPat-1]){
      if(lasttime < alltimes[iterPat]){break;}      // if after the last time we want to predict
      index_tempo++;
      time[index_tempo] = alltimes[iterPat];
    }
  }
  
 //// OPT- Efron correction

    if(Efron){
      
      double Wm1_tempo, sumRi, sumRi_di, di; // it is important that di is a double and not an in for the division in the for loop
      
      for(size_t iterEvent = 0 ; iterEvent < nEventsLast ; iterEvent++){
        
        if (death[iterEvent]>1){
          sumRi = sumEXb[iterEvent];
          sumRi_di = sumEXb_event[iterEvent];
          di = death[iterEvent];
            
          Wm1_tempo = 1/sumRi;
          
          for(int iterPat = 1; iterPat < di; iterPat++){
            Wm1_tempo += 1/(sumRi - sumRi_di*iterPat/di);
          }
         
          sumEXb[iterEvent] = di/Wm1_tempo;
        }
        
      }
      
    }
  
  //// 3- Computation of the hazards
  vector<double> hazard(nEventsLast, NA_REAL);
  vector<double> cumHazard(nEventsLast, NA_REAL);
  
  hazard[0] = death[0] / sumEXb[0];
  cumHazard[0] = hazard[0];
  
  for( size_t iterTime = 1 ; iterTime < nEventsLast ; iterTime ++){ // up to last time
    
    hazard[iterTime] = death[iterTime] / sumEXb[iterTime];
    cumHazard[iterTime] = cumHazard[iterTime-1] + hazard[iterTime];
  }
  
  //// 4- before and after the events
  if(addLast){
    time.push_back(alltimes[nPatients-1] + CST_EPSILON);
    hazard.push_back(NA_REAL);
    cumHazard.push_back(NA_REAL);
  }
  if(addFirst){
    time.insert(time.begin(),-CST_EPSILON);
    hazard.insert(hazard.begin(),0);
    cumHazard.insert(cumHazard.begin(),0);
  }
  
  //// export
  vector< vector<double> > res(3);
  res[0] = time;
  res[1] = hazard;
  res[2] = cumHazard;
  
  return res;
  
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
