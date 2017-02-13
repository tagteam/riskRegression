#' @title Extract i.i.d. decomposition from a Cox model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iid
#' 
#' @param object object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata Optional new data at which to do i.i.d. decomposition 
#' @param tauHazard the vector of times at which the i.i.d decomposition of the baseline hazard will be computed
#' @param keep.times Logical. If \code{TRUE} add the evaluation times to the output.
#' @param center.result Temporary argument. Should the IF be rescale to match timereg results.
#' @param ... additional arguments
#'
#' @details If there is no event in a strata, the influence function for the baseline hazard is set to 0.
#'
#' @examples
#' library(survival)
#' library(data.table)
#' set.seed(10)
#' d <- sampleData(7e1, outcome = "survival")[,.(eventtime,event,X1,X6)]
#' setkey(d, eventtime)
#' 
#' library(timereg)
#' mGS.cox <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6),
#'                      data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
#' ICbeta_GS <- mGS.cox$gamma.iid
#' IClambda_GS <- t(as.data.table(mGS.cox$B.iid))
#'  
#' m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
#' system.time(IC.cox <- iidCox(m.cox))
#' 
#' IC.cox <- iidCox(m.cox, tauHazard = 7)
#'  
#' 

#' @rdname iid
#' @export
iidCox <- function(object, newdata = NULL, tauHazard = NULL, 
                   keep.times = TRUE, center.result = TRUE){
  
  center.eXb <- TRUE # Temporary argument. Should the linear predictor be centered on the exponential scale.
  
  #### extract elements from object ####
  infoVar <- CoxVariableName(object)
  iInfo <- CoxVarCov(object)
  object.design <- CoxDesign(object)
  
  object.status <- object.design[,"status"]
  object.time <- object.design[,"stop"]
  object.strata <- CoxStrata(object, stratavars = infoVar$stratavars)
  object.levelStrata <- levels(object.strata)
  object.eXb <- exp(CoxLP(object, data = NULL, center = center.eXb))
  object.LPdata <- as.matrix(object.design[,infoVar$lpvars,drop = FALSE])
  nStrata <- length(levels(object.strata))
  
  object.center <- CoxCenter(object)
  
  #### Extract new observations ####
  if(!is.null(newdata)){
    
    if("data.frame" %in% class(newdata) == FALSE){
      stop("class of \'newdata\' must inherit from data.frame \n")
    }
    new.status <- newdata[[infoVar$status]]
    new.time <- newdata[[infoVar$time]]
    new.strata <- CoxStrata(object, data = newdata, 
                            sterms = infoVar$sterms, stratavars = infoVar$stratavars, levels = object.levelStrata, stratalevels = infoVar$stratalevels)
    new.eXb <- exp(CoxLP(object, data = newdata, center = center.eXb))
    new.LPdata <- model.matrix(object, newdata)
    
  }else{
    
    new.status <-  object.status
    new.time <-  object.time
    new.strata <-  object.strata
    new.eXb <- object.eXb
    new.LPdata <- object.LPdata
    
  }
  
  #### tests ####
  ## time at which the influence function is evaluated
  if(is.list(tauHazard) && length(tauHazard)!=nStrata){
    stop("argument \"tauHazard\" must be a list with ",nStrata," elements \n",
         "each element being the vector of times for each strata \n")
  }
  
  #### Compute quantities of interest ####
  p <- NCOL(object.LPdata)
  
    ## baseline hazard
  lambda0 <- predictCox(object, type = "hazard", centered = TRUE, keep.strata = TRUE)
  
  ## resale factor
  if(center.result == TRUE && !is.null(object.center)){
    scalingFactor <- exp(-as.double(coef(object) %*% object.center))
  }
  
  ## S0, E, jump times
  object.index_strata <- list() 
  object.order_strata <- list()
  
  object.eXb_strata <- list()
  object.LPdata_strata <- list()
  object.status_strata <- list()
  object.time_strata <- list()
  
  new.index_strata <- list()
  new.order_strata <- list()
  
  new.eXb_strata <- list()
  new.LPdata_strata <- list()
  new.status_strata <- list()
  new.time_strata <- list()
  
  Ecpp <- list()
  new.indexJump <- list()
  new.order <- NULL
  
  for(iStrata in 1:nStrata){
    
   ## reorder object data
    object.index_strata[[iStrata]] <- which(object.strata == object.levelStrata[iStrata])
    object.order_strata[[iStrata]] <- order(object.time[object.index_strata[[iStrata]]])
    
    object.eXb_strata[[iStrata]] <- object.eXb[object.index_strata[[iStrata]][object.order_strata[[iStrata]]]]
    object.LPdata_strata[[iStrata]] <- object.LPdata[object.index_strata[[iStrata]][object.order_strata[[iStrata]]],,drop = FALSE]
    object.status_strata[[iStrata]] <- object.status[object.index_strata[[iStrata]][object.order_strata[[iStrata]]]]
    object.time_strata[[iStrata]] <- object.time[object.index_strata[[iStrata]][object.order_strata[[iStrata]]]]
    
    ## reorder new data
    if(!is.null(newdata)){
      new.index_strata[[iStrata]] <- which(new.strata == object.levelStrata[iStrata])
      new.order_strata[[iStrata]] <- order(new.time[new.index_strata[[iStrata]]])
      
      new.eXb_strata[[iStrata]] <- new.eXb[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
      new.LPdata_strata[[iStrata]] <- new.LPdata[new.index_strata[[iStrata]][new.order_strata[[iStrata]]],,drop = FALSE]
      new.status_strata[[iStrata]] <- new.status[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
      new.time_strata[[iStrata]] <- new.time[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
    }else{
      new.index_strata[[iStrata]] <- object.index_strata[[iStrata]]
      new.order_strata[[iStrata]] <- object.order_strata[[iStrata]]
      
      new.eXb_strata[[iStrata]] <- object.eXb_strata[[iStrata]]
      new.LPdata_strata[[iStrata]] <- object.LPdata_strata[[iStrata]]
      new.status_strata[[iStrata]] <- object.status_strata[[iStrata]]
      new.time_strata[[iStrata]] <- object.time_strata[[iStrata]]
    }
    
    ## E
    Ecpp[[iStrata]] <-  calcE_cpp(status = object.status_strata[[iStrata]], 
                                  eventtime = object.time_strata[[iStrata]],
                                  eXb = object.eXb_strata[[iStrata]],
                                  X = object.LPdata_strata[[iStrata]],
                                  p = p, add0 = TRUE)
    
    new.indexJump[[iStrata]] <- prodlim::sindex(Ecpp[[iStrata]]$Utime1, new.time) - 1
    # if event/censoring is before the first event in the training dataset 
    # then sindex return 0 thus indexJump is -1
    # the following 3 lines convert -1 to 0
    if(any(new.indexJump[[iStrata]]<0)){
      new.indexJump[[iStrata]][new.indexJump[[iStrata]]<0] <- 0
    }
    
    ## store order
    if(length(new.order>0)){
      new.order <- c(new.order, new.index_strata[[iStrata]][new.order_strata[[iStrata]]])
    }else{
      new.order <- new.index_strata[[iStrata]][new.order_strata[[iStrata]]]
    }
    
  }
  
  #### Computation of the influence function ####
  ICbeta <- NULL
  ICcumhazard <- NULL
  ls.Utime1 <- NULL
  
  #### beta
  for(iStrata in 1:nStrata){
    
    new.indexJump_strata <- new.indexJump[[iStrata]][new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
    
    ## IF
    if(p>0){
      ICbeta_tempo <- ICbeta_cpp(newT = new.time_strata[[iStrata]], neweXb = new.eXb_strata[[iStrata]], newX = new.LPdata_strata[[iStrata]], newStatus = new.status_strata[[iStrata]], 
                                 newIndexJump = new.indexJump_strata, 
                                 S01 = Ecpp[[iStrata]]$S0, E1 = Ecpp[[iStrata]]$E, time1 = Ecpp[[iStrata]]$Utime1, iInfo = iInfo,
                                 p = p)    
    }else{
      ICbeta_tempo <- matrix(NA, ncol = 1, nrow = length(new.index_strata[[iStrata]]))
    }
    
    ## output
    ICbeta <- rbind(ICbeta, ICbeta_tempo)
    
  }
  
  ## set original order
  ICbeta <- ICbeta[order(new.order),,drop=FALSE]
  
  
  #### lambda
  for(iStrata in 1:nStrata){
    
    ## hazard
    if(nStrata==1){
      lambda0_strata <- lambda0
    }else{
      lambda0_strata <- lambda0[lambda0$strata == object.levelStrata[iStrata],, drop = FALSE]
    }
    
    ## tauHazard
    if(is.null(tauHazard)){
      tauHazard_strata <- object.time_strata[[iStrata]][object.status_strata[[iStrata]] == 1]
    }else if(is.list(tauHazard)){
      tauHazard_strata <- tauHazard[[nStrata]]
    }else{
      tauHazard_strata <- tauHazard
    }
    
    ## E
    nUtime1_strata <- length(Ecpp[[iStrata]]$Utime1)
    if(p>0){
      Etempo <- Ecpp[[iStrata]]$E[-NROW(Ecpp[[iStrata]]$E),,drop = FALSE]
    }else{
      Etempo <- matrix(0, ncol = 1, nrow = nUtime1_strata-1)
    }
    
    ## IF
    if(any(new.status_strata[[iStrata]]>0)){
      ICcumhazard_tempo <- IClambda0_cpp(tau = tauHazard_strata,
                                       ICbeta = ICbeta,
                                       newT = new.time, neweXb = new.eXb, newStatus = new.status, newIndexJump = new.indexJump[[iStrata]], newStrata = as.numeric(new.strata),
                                       S01 = Ecpp[[iStrata]]$S0[1:(nUtime1_strata-1)],
                                       E1 = Etempo,
                                       time1 = Ecpp[[iStrata]]$Utime1[1:(nUtime1_strata-1)],
                                       lambda0 = lambda0_strata[match(Ecpp[[iStrata]]$Utime1[-nUtime1_strata],lambda0_strata[,"time"]),"hazard"],
                                       p = p, strata = iStrata)
      
      # rescale
      if(center.result == TRUE && !is.null(CoxCenter(object))){
        ICcumhazard_tempo <- ICcumhazard_tempo * scalingFactor
      }
      
    }else{
      ICcumhazard_tempo <- matrix(0, ncol = 1, nrow = length(new.index_strata[[iStrata]]))
      if(length(tauHazard_strata)==0){tauHazard_strata <- NA}
    }
    
    # output 
    ls.Utime1 <- c(ls.Utime1, list(tauHazard_strata))
    if(keep.times){
      colnames(ICcumhazard_tempo) <- tauHazard_strata
    }
    ICcumhazard <- c(ICcumhazard, list(ICcumhazard_tempo))
    
  }
  
  #### export
  return(list(ICbeta = ICbeta,  
              ICcumhazard = ICcumhazard,
              time = ls.Utime1,
              indexObs = new.order
  ))
}