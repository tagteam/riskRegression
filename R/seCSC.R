
#' @title Standard error of the absolute risk predicted from cause-specific Cox models
#' @rdname seCSC
#'
#' @description  Standard error of the absolute risk predicted from cause-specific Cox models.
#' 
#' @param hazard list containing the baseline hazard for each cause in a matrix form. Columns correspond to the strata.
#' @param cumHazard list containing the cumulative baseline hazard for each cause in a matrix form. Columns correspond to the strata.
#' @param object.time a vector containing all the events regardless to the cause.
#' @param object.maxtime a matrix containing the latest event in the strata of the observation for each cause.
#' @param iid the value of the influcence function for each cause 
#' @param eXb_h a matrix containing the exponential of the linear predictor evaluated for the new observations (rows) for each cause (columns)
#' @param eXb_cumH same as before except when considering \code{survtype == "survival"}
#' @param new.LPdata a list of design matrices for the new observations for each cause.
#' @param new.strata a matrix containing the strata indicator for each observation and each cause.
#' @param times the time points at which to evaluate the predictions.  
#' @param new.n the number of new observations.
#' @param cause the cause of interest.
#' @param nCause the number of causes.
#'    
#' @examples 
#' \dontrun{
#' set.seed(10)
#' d <- SimCompRisk(2e1)
#' d$time <- round(d$time,1)
#' ttt <- unique(sort(d$time))#sort(sample(x = unique(sort(d$time)), size = 10))
#'
#' #### coxph function
#' CSC.fit <- CSC(Hist(time,event)~ X1+X2,data=d, method = "breslow")
#' 
#' predCSC <- predict(CSC.fit, newdata = d[1,,drop=FALSE], cause = 2, times = ttt, se = TRUE)
#' 
#' predC <- predictCox(CSC.fit$model[[2]],
#'                     newdata = d[1,,drop=FALSE], times = ttt, type = "hazard", se = FALSE)
#' 
#' 
#' predC <- predictCox(CSC.fit$model[[2]], newdata = d[1,],
#'                     times = ttt, type = "hazard", se = FALSE)
#' predC <- predictCox(CSC.fit$model[[2]],
#'                     newdata = d[1,], times = ttt,
#'                     type = c("hazard","cumHazard"), se = FALSE)
#' 
#' attr(predCSC,"se")-predC$hazard.se
#' predC$hazard.se
#' }
seCSC <- function(hazard, cumHazard, object.time, object.maxtime, iid,
                  eXb_h, eXb_cumH, new.LPdata, new.strata, times,
                  new.n, cause, nCause){
   
  nEtimes <- length(object.time)
  object.n <- NROW(iid[[1]]$ICbeta)
  CIF.se <- matrix(NA, nrow = new.n, ncol = length(times))
  
  for(iObs in 1:new.n){
    
    iStrata <- new.strata[iObs,]
    
    iHazard1 <- hazard[[cause]][,iStrata[cause]+1]
    iCumHazard <- rep(0, nEtimes)
    iIChazard1 <- NULL
    iICcumHazard <- matrix(0, nrow = object.n, ncol = nEtimes)
    
    for(iCause in 1:nCause){
      hazard_tempo <- hazard[[iCause]][,iStrata[iCause]+1]
      cumHazard_tempo <- cumHazard[[iCause]][,iStrata[iCause]+1]
      ICbeta_tempo <- iid[[iCause]]$ICbeta
      IChazard0_tempo <- iid[[iCause]]$IChazard[[iStrata[iCause]+1]]
      ICcumHazard0_tempo <- iid[[iCause]]$ICcumHazard[[iStrata[iCause]+1]]
      eXb_h_tempo <- eXb_h[iObs,iCause]
      eXb_cumH_tempo <- eXb_cumH[iObs,iCause]
      new.LPdata_tempo <- new.LPdata[[iCause]][iObs,,drop=FALSE]
      #
      iCumHazard <- iCumHazard + cumHazard_tempo
      
      #
      X_ICbeta <- ICbeta_tempo %*% t(new.LPdata_tempo)
      if(cause == iCause){iIChazard1 <- eXb_h_tempo*(IChazard0_tempo + X_ICbeta %*% hazard_tempo)}
      iICcumHazard <- iICcumHazard + eXb_cumH_tempo*(ICcumHazard0_tempo + X_ICbeta %*% cumHazard_tempo)
    }
    
   
    CIF.se_tempo <- rowCumSum(rowMultiply_cpp(iIChazard1 - rowMultiply_cpp(iICcumHazard, scale = iHazard1), scale = exp(-iCumHazard)))
    CIF.se_tempo <- cbind(0,CIF.se_tempo)[,prodlim::sindex(object.time, eval.times = times)+1]
    CIF.se[iObs,] <- sqrt(apply(CIF.se_tempo^2,2,sum))
    
  }
  
   return(CIF.se)
}
