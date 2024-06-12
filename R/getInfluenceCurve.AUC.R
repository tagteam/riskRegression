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
    if (t >= max(time)){conservative <- TRUE}
    # call c++
    conservativeIFcalculation <- getIC0AUC(time = time,
                                           status = event,
                                           tau = t,
                                           risk = risk,
                                           GTiminus = WTi,
                                           Gtau = Wt,
                                           auc = auc)
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


