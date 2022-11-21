predictCoxWeights <- function(object,
                              times,
                              newdata,
                              diag = FALSE,
                              weights, 
                              isBeforeTau = FALSE,
                              tau = -1){
  
  call <- match.call()
  ## centering
  centered <- TRUE
  if(inherits(centered,"data.frame")){
    df.reference <- centered
    centered2 <- TRUE ## for the linear predictor of the hazard
  }else{
    df.reference <- NULL
    centered2 <- centered ## for the linear predictor of the hazard
  }
  centered <- TRUE ## for the linear predictor of the baseline hazard
  
  ## ** Extract elements from object
  if (missing(times)) {
    nTimes <- 0
    times <- numeric(0)
  }else{
    nTimes <- length(times)
  }
  times.sorted <- times
  order.times <- 1:nTimes
  oorder.times <- 1:nTimes
  
  object.n <- coxN(object)  
  object.modelFrame <- coxModelFrame(object)
  infoVar <- coxVariableName(object, model.frame = object.modelFrame)
  object.baseEstimator <- coxBaseEstimator(object)
  
  ## ease access
  is.strata <- infoVar$is.strata
  object.levelStrata <- levels(object.modelFrame$strata) ## levels of the strata variable
  nStrata <- length(object.levelStrata) ## number of strata
  nVar.lp <- length(infoVar$lpvars) ## number of variables in the linear predictor
  
  ## ** normalize model frame
  ## convert strata to numeric
  object.modelFrame[,c("strata.num") := as.numeric(.SD$strata) - 1]
  
  ## linear predictor
  ## if we predict the hazard for newdata then there is no need to center the covariates
  object.modelFrame[,c("eXb") := exp(coxLP(object, data = NULL, center = if(is.null(newdata)){centered}else{FALSE}))]
  ## add linear predictor and remove useless columns
  rm.name <- setdiff(names(object.modelFrame),c("start","stop","status","eXb","strata","strata.num"))
  if(length(rm.name)>0){
    object.modelFrame[,c(rm.name) := NULL]
  }
  
  ## sort the data
  object.modelFrame[, c("statusM1") := 1-.SD$status] ## sort by statusM1 such that deaths appear first and then censored events
  object.modelFrame[, c("XXXindexXXX") := 1:.N] ## keep track of the initial positions (useful when calling calcSeCox)
  data.table::setkeyv(object.modelFrame, c("strata.num","stop","start","statusM1"))
  
  ## last event time in each strata
  if(!is.null(attr(times,"etimes.max"))){ ## specified by the user
    etimes.max <- attr(times,"etimes.max")
    attr(times,"etimes.max") <- NULL
    attr(times.sorted,"etimes.max") <- etimes.max
  }else if(is.strata){ ## based on the data
    if(nVar.lp==0){
      iDTtempo <- object.modelFrame[, .SD[which.max(.SD$stop)], by = "strata.num"]
      etimes.max <- iDTtempo[,if(.SD$status==1){1e12}else{.SD$stop}, by = "strata.num"][[2]]
    }else{
      etimes.max <- object.modelFrame[, max(.SD$stop), by = "strata.num"][[2]]
    }
  }else{
    if(nVar.lp==0 && (utils::tail(object.modelFrame$status,1)==1)){ ## no covariates and ends by a death
      etimes.max <- 1e12
    }else{
      etimes.max <- max(object.modelFrame[["stop"]])
    }
  }
  
  ## ** checks
  ## check user inputs 
  if(nTimes[1]>0 && any(is.na(times))){
    stop("Missing (NA) values in argument \'times\' are not allowed.\n")
  }
  ## predictCox is not compatible with all coxph/cph object (i.e. only handle only simple cox models)
  if(!is.null(object$weights) && !all(object$weights==1)){
    stop("predictCox does not know how to handle Cox models fitted with weights \n")
  }
  if(!is.null(object$naive.var)){
    stop("predictCox does not know how to handle frailty.") 
  }
  if(any(object.modelFrame[["start"]]!=0)){
    warning("The current version of predictCox was not designed to handle left censoring \n",
            "The function may be used on own risks \n") 
  }    
  if(object.baseEstimator == "exact"){
    stop("Prediction with exact handling of ties is not implemented.\n")
  }
  if(!is.null(object$call$tt)){
    stop("predictCox does not know how to handle time varying effects.\n") 
  }
  ## convergence issue
  if(!is.null(coef(object)) && any(is.na(coef(object)))){
    print(coef(object))
    stop("One or several parameters of the regression model have no value, i.e., a value 'NA'.\n")
  }
  ## prediction 
  if (missing(newdata)){
    stop("Argument 'newdata' is missing. Cannot compute standard errors in this case.")
  }
  if(!is.null(newdata)){
    if(nTimes[1]==0 && !identical(as.character("cumhazard"),"lp")){
      stop("Time points at which to evaluate the predictions are missing \n")
    }
    if(!is.vector(times)){
      stop("Argument \'times\' must be a vector \n")
    }
    name.regressor <- c(infoVar$lpvars.original, infoVar$stratavars.original)
    if(length(name.regressor) > 0 && any(name.regressor %in% names(newdata) == FALSE)){
      stop("Missing variables in argument \'newdata\': \"",
           paste0(setdiff(name.regressor,names(newdata)), collapse = "\" \""),
           "\"\n")
    }
  }
  ## diag argument
  if(!is.logical(diag)){
    stop("Argument \'diag\' must be of type logical \n")
  }
  if(diag){
    if(NROW(newdata)!=length(times)){
      stop("When argument \'diag\' is TRUE, the number of rows in \'newdata\' must equal the length of \'times\' \n")
    }
  }
  
  ## ** baseline hazard
  
  ## compute the baseline hazard
  Lambda0 <- baseHaz_cpp(starttimes = object.modelFrame$start,
                         stoptimes = object.modelFrame$stop,
                         status = object.modelFrame$status,
                         eXb = object.modelFrame$eXb,
                         strata = object.modelFrame$strata.num,
                         nPatients = object.n,
                         nStrata = nStrata,
                         emaxtimes = etimes.max,
                         predtimes = times.sorted,
                         cause = 1,
                         Efron = (object.baseEstimator == "efron"))
  ## restaure strata levels
  Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
  
  ## ** compute cumlative hazard and survival
  ## *** reformat newdata (compute linear predictor and strata)
  new.n <- NROW(newdata)
  ## newdata <- copy(newdata)
  setDT(newdata)
  
  Xb <- coxLP(object, data = newdata, center = FALSE)
  lp.iid <- FALSE
  new.eXb <- exp(Xb)
  new.strata <- coxStrata(object, data = newdata, 
                          sterms = infoVar$strata.sterms, 
                          strata.vars = infoVar$stratavars, 
                          strata.levels = infoVar$strata.levels)
  
  new.levelStrata <- levels(droplevels(new.strata))
  
  if(nVar.lp > 0){
    ## get the (new) design matrix
    new.LPdata <- model.matrix(object, data = newdata)
    if(NROW(new.LPdata)!=NROW(newdata)){
      stop("NROW of the design matrix and newdata differ. \n",
           "Maybe because newdata contains NA values \n")
    }
  }else{
    new.LPdata <- matrix(0, ncol = 1, nrow = new.n)
  }
  
  ## restaure original ordering
  data.table::setkeyv(object.modelFrame,"XXXindexXXX")
  if(diag){
    Lambda0$oorder.times <- oorder.times
  }else{
    Lambda0$oorder.times <- 1:nTimes
  }
  
  ## Computation of the influence function and/or the standard error
  if(diag){
    times2 <- times
  }else{
    times2 <- times.sorted
  }
  attr(times2,"etimes.max") <- attr(times.sorted,"etimes.max")
  times = times2
  object.time = object.modelFrame$stop
  object.eXb = object.modelFrame$eXb
  object.strata =  object.modelFrame$strata
  new.survival = NULL
  export = "iid"
  ##** Computation of the influence function
  iid.object <- iidCox(object, tau.hazard = times, store.iid = "minimal", return.object = FALSE)
  
  ## ** Prepare arguments
  if(diag){
    nTimes <- 1
  }
  new.strata <- as.numeric(new.strata)
  
  Lambda0$strata <- as.numeric(Lambda0$strata)    
  Lambda0$cumhazard <- lapply(1:nStrata,function(s){Lambda0$cumhazard[Lambda0$strata==s][Lambda0$oorder.times]})
  
  ## ** hazard / cumulative hazard / survival
  c(weightedAverageIFCumhazard_cpp(seqTau = times,
                                   cumhazard0 = Lambda0$cumhazard,
                                   newX = new.LPdata,
                                   neweXb = new.eXb,
                                   IFbeta = iid.object$IFbeta,
                                   cumEhazard0 = iid.object$calcIFhazard$cumElambda0,
                                   cumhazard_iS0 = iid.object$calcIFhazard$cumLambda0_iS0,
                                   delta_iS0 = iid.object$calcIFhazard$delta_iS0,
                                   sample_eXb = iid.object$calcIFhazard$eXb,
                                   sample_time = iid.object$obstime,
                                   indexJumpSample_time = lapply(iid.object$calcIFhazard$time1,
                                                                 function(iTime){prodlim::sindex(jump.times = iTime, eval.times = iid.object$obstime)-1}), 
                                   jump_time = iid.object$calcIFhazard$time1, 
                                   indexJumpTau = lapply(iid.object$calcIFhazard$time1,
                                                         function(iTime){prodlim::sindex(jump.times = iTime, eval.times = times)-1}), 
                                   lastSampleTime = iid.object$etime.max,
                                   newdata_index = lapply(1:nStrata,
                                                          function(iS){which(new.strata == iS)-1}),
                                   nTau = nTimes, nSample = object.n, nStrata = nStrata, p = nVar.lp,
                                   diag = diag,
                                   debug = 0, weights = weights, isBeforeTau = isBeforeTau, tau = tau))*new.n
}