getInfluenceCurve.AUC <- function(t,
                                  time,
                                  event,
                                  WTi,
                                  Wt,
                                  risk,
                                  MC,
                                  auc,
                                  nth.times,
                                  conservative,
                                  cens.model){
    ## assert that time is sorted in ascending order with stop
    if (is.unsorted(time)){
        stop("Internal error. Time is not sorted in ascending order. ")
    }
    
    conservativeIFcalculation <- getIC0AUC(time,event,t,risk,WTi,Wt,auc)
    if (conservative[[1]] || cens.model[[1]] == "none"){
        conservativeIFcalculation[["ic0"]]
    }
    else {
        conservativeIFcalculation[["ic0"]]+getInfluenceFunction.AUC.censoring.term(time = time,
                                                                                   event = event,
                                                                                   t = t,
                                                                                   IFcalculationList = conservativeIFcalculation,
                                                                                   MC = MC,
                                                                                   cens.model = cens.model, 
                                                                                   Wt = Wt,
                                                                                   auc = auc, 
                                                                                   nth.times = nth.times)
    }
}


