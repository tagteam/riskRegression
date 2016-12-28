#' @title Predicting hazard or cumulative hazard
#' @rdname predict.CauseSpecificCox
#' @aliases predict.CauseSpecificCox
#' @aliases predictBig.CauseSpecificCox
#'
#' @description  Apply formula to combine two or more Cox models into absolute risk (cumulative incidence function)
#' 
#' @param object The fitted cause specific Cox model
#' @param newdata A data frame containing the values of the variables
#'     in the right hand side of 'coxph' for each patient.
#' @param times Vector of times at which to return the estimated
#'     hazards/survival
#' @param cause Identifies the cause of interest among the competing events.
#' @param t0 the starting time for the conditional survival.
#' @param colnames Logical. If \code{TRUE} name the columns of the matrix containing the absolute risk with times
#' @param se Logical. If \code{TRUE} add the standard errors corresponding to the output. Experimental !!
#' @param ... not used
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds
#'     tag@@biostat.ku.dk
#' @details Note: for Cox regression models with time varying
#'     covariates it does not make sense to use this function, because
#'     the predicted risk has to be a measurable function of the data
#'     available at the time origin.
#' 
#' @return A data table containing the predictions for each patient (in rows) and each time (in columns). 
#' 
#' @examples 
#' 
#' d <- SimCompRisk(1e2)
#' d$time <- round(d$time,1)
#' ttt <- sample(x = unique(sort(d$time)), size = 10) 
#'
#' #### coxph function
#' CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow")
#' 
#' predCSC <- predict(CSC.fit, newdata = d, cause = 2, times = ttt)
#' 
#' #### cph function
#' CSC.cph <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow", fitter = "cph")
#' 
#' predcph <- predict(CSC.cph, newdata = d, cause = 2, times = ttt)
#'
#' #### conditional survival
#' T0 <- 1
#' predCSC_afterT0 <- predict(CSC.fit, newdata = d, cause = 2, times = ttt, t0 = T0)
#'
#' @method predict CauseSpecificCox
#' @export
predict.CauseSpecificCox <- function(object,newdata, times, cause, t0 = NA, colnames = TRUE, se  = FALSE, ...){
    survtype <- object$survtype
    if (length(cause) > 1){
        stop(paste0("Can only predict one cause. Provided are: ", 
                    paste(cause, collapse = ", "), sep = ""))
    }
    if (missing(cause)) {
        cause <- object$theCause
    }
    causes <- object$causes
    eTimes <- object$eventTimes
    if (any(match(as.character(cause), causes, nomatch = 0)==0L))
        stop(paste0("Requested cause ",as.character(cause)," does not match fitted causes which are:\n ",paste0("- ",causes,collapse="\n")))
    ## stopifnot(match(as.character(cause), causes, nomatch = 0) != 
    ## 0)
    if(se){
        stop("Not yet implemented \n")
    }
    if (survtype == "survival") {
        if (object$theCause != cause) 
            stop("Object can be used to predict cause ", object$theCause, 
                 " but not ", cause, ".\nNote: the cause can be specified in CSC(...,cause=).")
    }
    if(any(is.na(times))){
        stop("NA values in argument \'times\' \n")
    }
    if(length(t0)!=1){
        stop("\'t0\' must have length one \n")
    }
  
  # relevant event times to use
  eventTimes <- eTimes[which(eTimes <= max(times))] 
  if(length(eventTimes) == 0){eventTimes <- min(times)} # at least the first event
  
  # predict cumulative cause specific hazards
  nData <- NROW(newdata)
  nTimes <- length(eventTimes)
  if (survtype == "hazard") {
    nCause <- length(causes)
    
    ls.hazard <- vector(mode = "list", length = nCause)
    ls.cumHazard <- vector(mode = "list", length = nCause)
    M.eXb_h <- matrix(NA, nrow = nData, ncol = nCause)
    M.strata <- matrix(NA, nrow = nData, ncol = nCause)
    M.etimes.max <- matrix(NA, nrow = nData, ncol = nCause)
    
    for(iterC in 1:nCause){
        infoVar <- CoxVariableName(object$models[[iterC]])
      
      ## baseline hazard from the Cox model
      causeBaseline <- predictCox(object$models[[iterC]],
                                  times = eventTimes, newdata = NULL,
                                  type = c("hazard","cumHazard"), 
                                  keep.strata = TRUE, keep.lastEventTime = TRUE, keep.times = TRUE,
                                  se = FALSE, format = "list")
      
      ls.hazard[[iterC]] <- matrix(causeBaseline$hazard, byrow = FALSE, nrow = nTimes)
      ls.cumHazard[[iterC]] <- matrix(causeBaseline$cumHazard, byrow = FALSE, nrow = nTimes) 
      
      ## linear predictor for the new observations
      newdata.eXb <- exp(CoxLP(object$models[[iterC]], data = newdata, center = TRUE))
      
      ## strata for the new observations
      newdata.strata <- CoxStrata(object$models[[iterC]], data = newdata, 
                                  sterms = infoVar$sterms, 
                                  stratavars = infoVar$stratavars, 
                                  levels = levels(causeBaseline$strata), 
                                  stratalevels = infoVar$stratalevels)
      
      M.strata[,iterC] <- as.numeric(newdata.strata)-1
      M.eXb_h[,iterC] <- newdata.eXb
      
      ## last time by strata
      M.etimes.max[,iterC] <- causeBaseline$lastEventTime[M.strata[,iterC]+1]
    }
    
    M.eXb_cumH <- M.eXb_h
    
  }else{
    nCause <- 1
    tdiff <- 0*min(diff(eTimes))/2 # TO MATCH test-CauseSpecificCoxRegresion.R but will not match pec
    
    #### cause ####
    infoVar_Cause <- CoxVariableName(object$models[[paste("Cause",cause)]])
    
    ## baseline hazard from the Cox model
    causeBaseline <- predictCox(object$models[[paste("Cause",cause)]],
                                times = eventTimes, newdata = NULL,
                                type = c("hazard","cumHazard"), 
                                keep.strata = TRUE, keep.lastEventTime = TRUE, keep.times = TRUE,
                                se = FALSE, format = "list")
    
    ls.hazard <- list(matrix(causeBaseline$hazard, byrow = FALSE, nrow = nTimes))
    
    ## linear predictor for the new observations
    newdata.eXb_cause <- exp(CoxLP(object$models[[paste("Cause",cause)]], data = newdata, center = TRUE))
    M.eXb_h <- cbind(newdata.eXb_cause)
    
    ## strata for the new observations
    newdata.strata <- CoxStrata(object$models[[paste("Cause",cause)]], data = newdata, 
                                sterms = infoVar_Cause$sterms, 
                                stratavars = infoVar_Cause$stratavars, 
                                levels = levels(causeBaseline$strata), 
                                stratalevels = infoVar_Cause$stratalevels)
    
    M.strata <- cbind(as.numeric(newdata.strata)-1)
    M.etimes.max <- cbind(causeBaseline$lastEventTime[M.strata+1]) # last time by strata
    
    #### overall ####

    ## baseline
    overallBaseline <- predictCox(object$models[["OverallSurvival"]],
                                  times = eventTimes-tdiff, newdata = NULL,
                                  type = "cumHazard", 
                                  keep.strata = TRUE, keep.lastEventTime = TRUE, keep.times = TRUE,
                                  se = FALSE, format = "list")
    ls.cumHazard <- list(matrix(overallBaseline$cumHazard, byrow = FALSE, nrow = nTimes))
    
    ## linear predictor for the new observations
    newdata.eXb_All <- exp(CoxLP(object$models[["OverallSurvival"]], data = newdata, center = TRUE))
    M.eXb_cumH <- cbind(newdata.eXb_All)
    
  }
  CIF <- predictCIF_cpp(hazard = ls.hazard, 
                        cumHazard = ls.cumHazard, 
                        eXb_h = M.eXb_h, 
                        eXb_cumH = M.eXb_cumH, 
                        strata = M.strata,
                        newtimes = sort(times), 
                        etimes = eventTimes, 
                        etimeMax = apply(M.etimes.max,1,min), 
                        t0 = t0,
                        nTimes = nTimes,
                        nNewTimes = length(times), 
                        nData = nData,
                        cause = which(causes == cause) - 1, 
                        nCause = nCause)
  
  #### export ###
  if(any(order(times) != 1:length(times))){# reorder times
    CIF <- CIF[,order(order(times))]
  }
  if(colnames){
    colnames(CIF) <- times
  }
  return(CIF)
}

seCSC <- function(indexTimes,
                  hazardCause1, survivalCause1, survivalCause2,
                  hazardCause1.se, survivalCause1.se, survivalCause2.se){
  
  integrand <- (hazardCause1*survivalCause1*survivalCause2.se)^2 + (hazardCause1*survivalCause1.se*survivalCause2)^2 + (hazardCause1.se*survivalCause1*survivalCause2.se)^2
  out <- c(0, rowCumSum(integrand))[,indexTimes + 1,drop = FALSE] # it is unclear to which value should be set the variance of the absolute risk before the first event. 0???
  return(out)
}

