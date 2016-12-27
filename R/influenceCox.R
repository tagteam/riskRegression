#' @title Extract i.i.d. decomposition from a Cox model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iid
#' 
#' @param object object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata Optional new data at which to do i.i.d. decomposition 
#' @param tauLambda the vector of times at which the i.i.d decomposition of the baseline hazard will be computed
#' @param keep.times Logical. If \code{TRUE} add the evaluation times to the output.
#' @param keep.indexObs Logical. If \code{TRUE} set the names of the rows in the output to the position of the observations in the original dataset.  
#' @param keep.originalOrder Temporary argument. Should the linear predictor be centered when computing the baseline hazard (equivalent to centered in basehaz).
#' @param center.result Temporary argument. Should the IF be rescale to match timereg results.
#' @param ... additional arguments
#'
#' @examples
#' library(survival)
#' library(data.table)
#' set.seed(10)
#' d <- sampleData(7e1, outcome = "survival")[,.(eventtime,event,X1,X6)]
#' setkey(d, eventtime)
#' 
#' library(timereg)
#' mGS.cox <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d, resample.iid = TRUE, max.timepoint.sim=NULL)
#' ICbeta_GS <- mGS.cox$gamma.iid
#' IClambda_GS <- t(as.data.table(mGS.cox$B.iid))
#'  
#' m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE)
#' system.time(IC.cox <- iidCox(m.cox))
#' 
#' IC.cox <- iidCox(m.cox, tauLambda = 7)
#'  
#' 

#' @rdname iid
#' @export
iidCox <- function(object, newdata = NULL, tauLambda = NULL, 
                   keep.times = TRUE, keep.indexObs = TRUE, keep.originalOrder = FALSE,
                  center.result = TRUE){
  
  center.eXb <- TRUE # Temporary argument. Should the linear predictor be centered on the exponential scale.
  center.LPdata <- FALSE # Temporary argument. Should the covariates be centered.
  
  #### extract elements from object ####
  infoVar <- CoxStrataVar(object)
  iInfo <- CoxVarCov(object)
  is.strata <- infoVar$is.strata
  object.designVar <- colnames(iInfo)
  
  object.status <- CoxStatus(object)
  object.time <- CoxEventtime(object)
  object.strata <- CoxStrata(object, stratavars = infoVar$stratavars)
  object.levelStrata <- levels(object.strata)
  object.eXb <- exp(CoxLP(object, data = NULL, center = center.eXb))
  object.LPdata <- as.matrix(CoxDesign(object, data = CoxData(object), 
                                       lpvars = infoVar$lpvars, stratavars = infoVar$stratavars,
                                       rm.intercept = TRUE, center = center.LPdata))
  nStrata <- length(levels(object.strata))
  
  # for factor variables: retain only the relevant column from the design matrix, 
  # i.e. remove the column(s) corresponding to the reference level(s)
  object.LPdata <- object.LPdata[, object.designVar, drop = FALSE]
  
  #### Extract new observations ####
  if(!is.null(newdata)){
    info2 <- CoxResponseVar(object)
    new.time <- newdata[[info2$time]]
    new.strata <- CoxStrata(object, data = newdata, 
                            sterms = infoVar$sterms, stratavars = infoVar$stratavars, levels = object.levelStrata, stratalevels = infoVar$stratalevels)
    new.status <- newdata[[info2$event]]
    new.eXb <- exp(CoxLP(object, data = newdata, center = center.eXb))
    new.LPdata <- as.matrix(CoxDesign(object, data = newdata,  
                                      lpvars = infoVar$lpvars, stratavars = infoVar$stratavars,
                                      rm.intercept = TRUE, center = center.LPdata))
    
    # for factor variables: retain only the relevant column from the design matrix, 
    new.LPdata <- new.LPdata[, object.designVar, drop = FALSE]
  }
  
  
  #### Compute quantities of interest only for event times ####
  lastEventtime <- object.time[length(object.time)]
  p <- NCOL(object.LPdata)
  
  ## baseline hazard
  lambda0 <- predictCox(object, type = "hazard", centered = TRUE, keep.strata = TRUE)
  
  ## resale factor
  if(center.result == TRUE && !is.null(CoxCenter(object))){
    scalingFactor <- exp(-as.double(coef(object) %*% CoxCenter(object)))
  }
  
  ## time at which the influence function is evaluated
  if(is.list(tauLambda) && length(tauLambda)!=nStrata){
    stop("argument \"tauLambda\" must be a list with ",nStrata," elements \n",
         "each element being the vector of times for each strata \n")
  }
  
  #### Computation of the influence function ####
  ICbeta <- NULL
  ICLambda0 <- NULL
  new.order <- NULL
  ls.Utime1 <- NULL
  
  for(iStrata in 1:nStrata){
  
    ## reorder object data
    object.index_strata <- which(object.strata == object.levelStrata[iStrata])
    object.order_strata <- order(object.time[object.index_strata])
  
    object.eXb_strata <- object.eXb[object.index_strata[object.order_strata]]
    object.LPdata_strata <- object.LPdata[object.index_strata[object.order_strata],,drop = FALSE]
    object.status_strata <- object.status[object.index_strata[object.order_strata]]
    object.time_strata <- object.time[object.index_strata[object.order_strata]]
    
    ## reorder new data
    if(is.null(newdata)){
      new.index_strata <- object.index_strata
      new.order_strata <- object.order_strata
      
      new.time_strata <- object.time_strata
      new.status_strata <- object.status_strata
      new.eXb_strata <- object.eXb_strata
      new.LPdata_strata <- object.LPdata_strata
    }else{
      new.index_strata <- which(new.strata == object.levelStrata[iStrata])
      new.order_strata <- order(new.time[new.index_strata])
      
      new.eXb_strata <- new.eXb[new.index_strata[new.order_strata]]
      new.LPdata_strata <- new.LPdata[new.index_strata[new.order_strata],,drop = FALSE]
      new.status_strata <- new.status[new.index_strata[new.order_strata]]
      new.time_strata <- new.time[new.index_strata[new.order_strata]]
    }
    new.n_strata <- length(new.index_strata)
    
    ## hazard
    if(nStrata==1){
      lambda0_strata <- lambda0
    }else{
      lambda0_strata <- lambda0[lambda0$strata == object.levelStrata[iStrata],, drop = FALSE]
    }
    
    ## tauLambda
    if(is.null(tauLambda)){
      tauLambda_strata <- object.time_strata[object.status_strata == 1]
    }else if(is.list(tauLambda)){
      tauLambda_strata <- tauLambda[[nStrata]]
    }else{
      tauLambda_strata <- tauLambda
    }
    
    ## E
    Ecpp <-  calcE_cpp(status = object.status_strata, 
                       eventtime = object.time_strata,
                       eXb = object.eXb_strata,
                       X = object.LPdata_strata,
                       p = p, add0 = TRUE)
    
    new.indexJump_strata <- prodlim::sindex(Ecpp$Utime1, new.time_strata) - 1
    
    nUtime1_strata <- length(Ecpp$Utime1)
    if(p>0){
      Etempo <- Ecpp$E[-NROW(Ecpp$E),,drop = FALSE]
    }else{
      Etempo <- matrix(0, ncol = 1, nrow = nUtime1_strata-1)
    }
    
    ## IC BETA
    if(p>0){
      ICbeta_tempo <- ICbeta_cpp(newT = new.time_strata, neweXb = new.eXb_strata, newX = new.LPdata_strata, newStatus = new.status_strata, 
                                 newIndexJump = new.indexJump_strata, 
                                 S01 = Ecpp$S0, E1 = Ecpp$E, time1 = Ecpp$Utime1, iInfo = iInfo,
                                 p = p)    
    }else{
      ICbeta_tempo <- matrix(NA, ncol = 1, nrow = new.n_strata)
    }
    
    ## IC LAMBDA
    IClambda0_tempo <- IClambda0_cpp(tau = tauLambda_strata, 
                                     ICbeta = ICbeta_tempo,
                                     newT = new.time_strata, neweXb = new.eXb_strata, newStatus = new.status_strata, newIndexJump = new.indexJump_strata, 
                                     S01 = Ecpp$S0[1:(length(Ecpp$Utime1)-1)], 
                                     E1 = Etempo, 
                                     time1 = Ecpp$Utime1[1:(nUtime1_strata-1)], 
                                     lambda0 = lambda0_strata[match(Ecpp$Utime1[-nUtime1_strata],lambda0_strata[,"time"]),"hazard"], 
                                     p = p)
   
    
    # rescale
    if(center.result == TRUE && !is.null(CoxCenter(object))){
      IClambda0_tempo <- IClambda0_tempo * scalingFactor
    }
    
    # store order
      if(length(new.order>0)){
        new.order <- c(new.order, new.index_strata[new.order_strata])
      }else{
        new.order <- new.index_strata[new.order_strata]
      }
      ls.Utime1 <- c(ls.Utime1, list(Ecpp$Utime1[1:(nUtime1_strata-1)]))
    
    # output 
    if(keep.times){
      colnames(IClambda0_tempo) <- tauLambda_strata
    }
    if(keep.indexObs){
      rownames(IClambda0_tempo) <- new.index_strata[new.order_strata]
      rownames(ICbeta_tempo) <- new.index_strata[new.order_strata]
    }
    ICbeta <- rbind(ICbeta, ICbeta_tempo)
    ICLambda0 <- c(ICLambda0, list(IClambda0_tempo))
  }

  #### set original order
  if(keep.originalOrder){
    original.order <- order(new.order)
    
    ICbeta <- ICbeta[original.order,,drop=FALSE]
    
    AllTau <- sort(unlist(ls.Utime1))
    
    ICLambda0_all  <- lapply(1:nStrata, function(iStrata){
      jump <- prodlim::sindex(eval.times = AllTau, jump.times = c(0,ls.Utime1[[iStrata]],tail(ls.Utime1[[iStrata]],1)+1e-12))
      cbind(0,ICLambda0[[iStrata]],NA)[,jump, drop = FALSE]
    })
    ICLambda0_all <- do.call(rbind,ICLambda0_all)
    
    if(keep.times){
      colnames(ICLambda0_all) <- unlist(ls.Utime1)
    }
 
    ICLambda0 <- ICLambda0_all[original.order,,drop=FALSE]
  }
  
  #### export
  return(list(ICbeta = ICbeta,  # restaure original ordering
              ICLambda0 = ICLambda0,
              time = ls.Utime1,
              indexObs = new.order
  ))
}