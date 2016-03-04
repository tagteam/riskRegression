// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "sortS.h"

using namespace Rcpp;
using namespace std;

// Q1 is that better to do 1 and 2 (two loops one to get the number of patient by strata then create the vector for each strata) or use push.back
// Q2 export from BaseHaz_cpp may be improved

vector< vector<double> > BaseHaz_cpp(const vector<double>& alltimes, const vector<int>& status, const vector<double>& Xb, 
                                     int nPatients, double lasttime, int cause, bool Efron, bool addFirst, bool addLast);

//' @export
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
  vector<int> death(nEventsLast,0); // could be avoided by only computing the number of deaths at the current time in 3
  vector<int> event(nEventsLast,0); // only for Efron correction
  vector<double> sumEXb(nEvents);
  
  // first observation
  time[0] = alltimes[0];
  event[0]++;  // only for Efron correction
  death[0] = (status[0]==cause);
  sumEXb[nEvents-1] = exp(Xb[nPatients-1]);
  
  // remaining observations
  size_t index_tempo = nEvents-1;
  for(int iterPat = nPatients-2 ; iterPat >= 0 ; iterPat--){
    if(alltimes[iterPat] != alltimes[iterPat+1]){
      index_tempo--;
      sumEXb[index_tempo] = sumEXb[index_tempo+1];
    }
    sumEXb[index_tempo] += exp(Xb[iterPat]);
  }
  
  index_tempo = 0;
  for(int iterPat = 1 ; iterPat < nPatients ; iterPat++){
    if(alltimes[iterPat] != alltimes[iterPat-1]){
      if(lasttime < alltimes[iterPat]){break;}      // if after the last time we want to predict
      index_tempo++;
      time[index_tempo] = alltimes[iterPat];
    }
    death[index_tempo] += (status[iterPat]==cause);
    event[index_tempo]++;  // only for Efron correction
  }
  
//   //// OPT- Efron correction
//   if(Efron){
//     
//     double correction, Wtempo, Wbase, Wdiff;
//     int countEvent= 0; 
//     sumEXb.push_back(0);
//     
//     index_tempo = 0;
//     Wbase = sumEXb[index_tempo];
//     Wtempo = 1;
//     Wdiff = sumEXb[index_tempo] - sumEXb[index_tempo+1];
//     
//     countEvent ++;
//     correction = Wdiff * (countEvent-1)/((double) event[index_tempo]);
//     Wtempo *= pow(Wbase - 0*correction, 1/((double) event[index_tempo]));
//     
//     for(int iterPat = 0 ; iterPat < nPatients ; iterPat++){ 
//       
//       if(alltimes[iterPat] != alltimes[iterPat-1]){
//         sumEXb[index_tempo] = Wbase;// Wtempo;
//         index_tempo++;
//         Wbase = sumEXb[index_tempo];
//         Wtempo = 1;
//         Wdiff = sumEXb[index_tempo] - sumEXb[index_tempo+1];
//         countEvent = 0;
//       }
//       countEvent ++;
//       correction = Wdiff * (countEvent-1)/((double) event[index_tempo]);
//      
//       Wtempo *= pow(Wbase - 0*correction, 1/((double) event[index_tempo]));
//       Rcout << nEvents << " " <<  index_tempo << " "  << event[index_tempo]  << " " << 1/((double) event[index_tempo]) << " "  << Wbase << " " << Wtempo << endl;
//     }
//     
//     sumEXb[index_tempo] = Wtempo;
//   }           
      
//       for(size_t iterTime = 0 ; iterTime < nEventsLast-1 ; iterTime ++){ 
//         
//         correction = event[iterTime] * (event[iterTime]-1)/((double)2 * event[nEventsLast-1]);
//         sumEXb[iterTime] = sumEXb[iterTime] - (sumEXb[iterTime] - sumEXb[iterTime+1]) * correction;
//         
//       }
      
//       correction = event[nEventsLast-1] * (event[nEventsLast-1]-1)/((double)2 * event[nEventsLast-1]);
//       sumEXb[nEventsLast-1] = sumEXb[nEventsLast-1] - sumEXb[nEventsLast-1] * correction;
     

    
    

//   if(Efron){
//     index_tempo = nEvents-1;
//     int countEvent= 1;
//     double Wtempo = 0;
//     
//     for(int iterPat = nPatients-2 ; iterPat >= 0 ; iterPat--){
//       
//       if(alltimes[iterPat] != alltimes[iterPat+1]){
//         Wtempo = ((countEvent+1)/2.0 - 1)/((double)countEvent);
//         // Rcout <<  countEvent << " " << Wtempo << endl;
//         sumEXb[index_tempo] = sumEXb[index_tempo] - (sumEXb[index_tempo] - sumEXb[index_tempo+1]) * Wtempo;
//         
//         index_tempo--;
//         countEvent = 1;
//         
//       }else{
//         countEvent++;
//       }
//     }
//     
//     if(alltimes[0] == alltimes[1]){
//       Wtempo = ((countEvent+1)/2.0 - 1)/((double)countEvent); 
//       //Rcout <<  countEvent << " " << Wtempo << endl;
//       sumEXb[index_tempo] =  sumEXb[index_tempo] - (sumEXb[index_tempo] - sumEXb[index_tempo+1]) * Wtempo;//pow(Wtempo,1/(double)death[index_tempo]);
//     }
//   }
  
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
    time.push_back(time[nEventsLast-1]);
    hazard.push_back(NA_REAL);
    cumHazard.push_back(NA_REAL);
  }
  if(addFirst){
    time.insert(time.begin(),0);
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


