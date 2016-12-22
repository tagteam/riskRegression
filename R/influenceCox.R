#' @title Extract i.i.d. decomposition from a Cox model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iid
#' 
#' @param object object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata Optional new data at which to do i.i.d. decomposition 
#' @param tauLambda the vector of times at which the i.i.d decomposition of the baseline hazard will be computed
#' @param tauBeta the time at which the i.i.d decomposition of the beta will be computed
#' @param keep.times Logical. If \code{TRUE} add the evaluation times to the output. 
#' @param center.eXb Temporary argument. Should the linear predictor be centered on the exponential scale.
#' @param center.LPdata Temporary argument. Should the covariates be centered.
#' @param center.lambda0 Temporary argument. Should the linear predictor be centered when computing the baseline hazard (equivalent to centered in basehaz).
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
iidCox <- function(object, newdata = NULL, tauLambda = NULL, tauBeta = NULL, keep.times = TRUE, 
                   center.eXb = TRUE, center.LPdata = FALSE, center.lambda0 = TRUE, center.result = TRUE){
  ## TODO
  # to be adapted for stratified Cox models
  # check if it is ok for categorical variables
  
  #### extract elements from object ####
  infoVar <- CoxStrataVar(object)
  iInfo <- CoxVarCov(object)
  is.strata <- infoVar$is.strata
  object.designVar <- colnames(iInfo)
  
  object.n <- CoxN(object)
  object.status <- CoxStatus(object)
  object.time <- CoxEventtime(object)
  object.strata <- CoxStrata(object, stratavars = infoVar$stratavars)
  object.levelStrata <- levels(object.strata)
  object.eXb <- exp(CoxLP(object, data = NULL, center = center.eXb))
  object.LPdata <- as.matrix(CoxDesign(object, data = CoxData(object), 
                                       lpvars = infoVar$lpvars, stratavars = infoVar$stratavars,
                                       rm.intercept = TRUE, center = center.LPdata))
  nStrata <- length(levels(object.strata))
  
  # reorder data
  dtOrder <- data.table(strata = object.strata, time = object.time)
  dtOrder[, index := 1:.N]
  setkeyv(dtOrder, cols = c("strata","time"))
  object.order <- dtOrder$index
  
  object.eXb <- object.eXb[object.order]
  object.LPdata <- object.LPdata[object.order,,drop = FALSE]
  object.status <- object.status[object.order]
  object.strata <- object.strata[object.order]
  object.time <- object.time[object.order]
  
  object.indexStrata <- lapply(object.levelStrata, function(level){which(object.strata==level)-1})
  
  if(is.null(tauLambda)){
    dtTempo <- data.table(status = object.status, time = object.time, strata = object.strata)
    tauLambda <- dtTempo[,.(lsTime = list(time[status == 1])),by = strata][["lsTime"]]
  }else{
    if(!is.list(tauLambda)){
      tauLambda <- tauLambda[[1]]
      tauLambda <- lapply(1:nStrata, function(i){tauLambda})
    }else if(length(tauLambda)!=nStrata){
      stop("argument \"tauLambda\" must be a list with ",nStrata," elements \n",
           "each element being the vector of times for each strata \n")
    }
  }
  if(is.null(tauBeta)){tauBeta <- tail(object.time[object.status>0],1)}

  if(any(unlist(tauLambda) > tauBeta)){
    stop("Elements in tauLambda must not exceed tauBeta \n")
  }
  
  # for factor variables: retain only the relevant column from the design matrix, 
  # i.e. remove the column(s) corresponding to the reference level(s)
  object.LPdata <- object.LPdata[, object.designVar, drop = FALSE]
  
  #### Extract new observations
  
  if(is.null(newdata)){
    new.n <- object.n
    new.time <- object.time
    new.strata <- object.strata
    new.status <- object.status
    new.eXb <- object.eXb
    new.LPdata <- object.LPdata
    new.indexStrata <- object.indexStrata
    # reorder data (already done)
    new.order <- object.order
  }else{
    info2 <- CoxResponseVar(object)
    new.n <- NROW(newdata)
    new.time <- newdata[[info2$time]]
    new.strata <- CoxStrata(object, data = newdata, 
                            sterms = infoVar$sterms, stratavars = infoVar$stratavars, levels = object.levelStrata, stratalevels = infoVar$stratalevels)
    new.status <- newdata[[info2$event]]
    new.eXb <- exp(CoxLP(object, data = newdata, center = center.eXb))
    new.LPdata <- as.matrix(CoxDesign(object, data = newdata,  
                                      lpvars = infoVar$lpvars, stratavars = infoVar$stratavars,
                                      rm.intercept = TRUE, center = center.LPdata))
    
    # reorder data
    dtOrder <- data.table(strata = new.strata, time = new.time)
    dtOrder[, index := 1:.N]
    setkeyv(dtOrder, cols = c("strata","time"))
    object.order <- dtOrder$index
    
    new.eXb <- new.eXb[new.order]
    new.LPdata <- new.LPdata[new.order,,drop = FALSE]
    new.status <- new.status[new.order]
    new.strata <- new.strata[new.order]
    new.time <- new.time[new.order]
    
    new.indexStrata <-lapply(object.levelStrata, function(level){which(new.strata==level)-1})
    
    # for factor variables: retain only the relevant column from the design matrix, 
    new.LPdata <- new.LPdata[, object.designVar, drop = FALSE]
  }
  
  
  #### Compute quantities of interest only for event times
  lastEventtime <- object.time[length(object.time)]
  p <- NCOL(object.LPdata)
  
  ## baseline hazard
  lambda0 <- predictCox(object, type = "hazard", centered = center.lambda0, keep.strata = TRUE)
  
  ## E
  Ecpp <- calcEstrata_cpp(tau = tauBeta,
                          indexStrata = object.indexStrata,
                          status = object.status,
                          p = p, nStrata = nStrata,
                          eventtime = object.time,
                          eXb = object.eXb,
                          X = object.LPdata,
                          add0 = TRUE)
  
  ## resale
  if(center.result == TRUE && !is.null(CoxCenter(object))){
    scalingFactor <- exp(-as.double(coef(object) %*% CoxCenter(object)))
  }
  
  #### Computation of the influence function
  ICbeta <- NULL
  ICLambda0 <- NULL
  
  for(iStrata in 1:nStrata){
  
    ## IC BETA
    ICbeta_tempo <- iidBeta(tau = tauBeta,
                            newT = new.time[new.indexStrata[[iStrata]]+1],
                            neweXb = new.eXb[new.indexStrata[[iStrata]]+1],
                            newX = new.LPdata[new.indexStrata[[iStrata]]+1,,drop = FALSE],
                            newStatus = new.status[new.indexStrata[[iStrata]]+1],
                            S01 = Ecpp$S0[[iStrata]],
                            E1 = matrix(Ecpp$E[[iStrata]], ncol = p, byrow = FALSE),
                            time1 = Ecpp$Utime1[[iStrata]],
                            n = Ecpp$n[iStrata], p = p, iInfo = iInfo)
    
    ICbeta <- rbind(ICbeta, ICbeta_tempo)
    
    ## IC LAMBDA
    if(p==0){
      E1tempo <- matrix(0, ncol = 1, nrow = Ecpp$n_Utime1[iStrata]-1)
      Xtempo <-matrix(0, ncol = 1, nrow = Ecpp$n_Utime1[iStrata]-1)
    }else{
      E1tempo <- matrix(Ecpp$E[[iStrata]], ncol = p, byrow = FALSE)[1:(Ecpp$n_Utime1[iStrata]-1),,drop = FALSE]
      Xtempo <- new.LPdata[new.indexStrata[[iStrata]]+1,,drop = FALSE]
    }
    
    IClambda0_tempo <- iidLambda0(tau = tauLambda[[iStrata]], max.time = lastEventtime,
                                  newT = new.time[new.indexStrata[[iStrata]]+1],
                                  neweXb = new.eXb[new.indexStrata[[iStrata]]+1],
                                  newX = Xtempo,
                                  newStatus = new.status[new.indexStrata[[iStrata]]+1],
                                  ICbeta = ICbeta_tempo,
                                  S01 = Ecpp$S0[[iStrata]][1:(Ecpp$n_Utime1[iStrata]-1)],
                                  E1 = E1tempo,
                                  time1 = Ecpp$Utime1[[iStrata]][1:(Ecpp$n_Utime1[iStrata]-1)],
                                  lambda0 = lambda0[match(Ecpp$Utime1[[iStrata]][1:(Ecpp$n_Utime1[iStrata]-1)],lambda0[,"time"]),"hazard"],
                                  n = Ecpp$n[iStrata], p = p)
   
    
    # rescale
    if(center.result == TRUE && !is.null(CoxCenter(object))){
      IClambda0_tempo <- IClambda0_tempo * scalingFactor
    }
    
    # reorder
    indexTempo <- new.order[new.strata==levels(new.strata)[iStrata]]
    IClambda0_tempo <- IClambda0_tempo[indexTempo-(indexTempo[1]-1),,drop=FALSE]
    
    # output 
    if(keep.times){
       colnames(IClambda0_tempo) <- tauLambda[[iStrata]]
    }
    
    ICLambda0 <- c(ICLambda0, list(IClambda0_tempo))
  }

  ## rename 
  if(keep.times){
    attr(ICbeta,"time") <- tauBeta
  }
 
  #### export
  return(list(ICbeta = ICbeta[new.order,,drop=FALSE],  # restaure original ordering
              ICLambda0 = ICLambda0
  ))
}

iidBeta <- function(tau, newT, neweXb, newX, newStatus,
                    S01, E1, time1, 
                    n, p, iInfo){
  
  if(is.null(iInfo)){return(matrix(0, nrow = n, ncol = 1))}
  IC <- matrix(NA, nrow = n, ncol = p)
  
  for(iterI in 1:n){
    if(newT[iterI]<=tau){
      term1 <- calcU_cpp(newX = newX[iterI,,drop = FALSE], newStatus = newStatus[iterI], newN = 1,
                         IndexNewT = 0, ENewT = E1[prodlim::sindex(time1, newT[iterI]),,drop = FALSE]*newStatus[iterI],
                         p = p, aggregate = FALSE)
    }else{
      term1 <- matrix(0, nrow = 1, ncol = p)
    }
    
    term2 <- -colScale_cpp(rowCenter_cpp(E1[time1<=min(newT[iterI],tau),,drop = FALSE], center = newX[iterI,,drop = FALSE]), 
                           scale = S01[time1<=min(newT[iterI],tau)])
    
    IC[iterI,] <- as.double(iInfo %*% t(term1 - neweXb[iterI] * colSums(term2)))
  }
  
  return(IC)
}

iidLambda0 <- function(tau, max.time,
                       newT, neweXb, newX, newStatus, ICbeta,
                       S01, E1, time1, lambda0,
                       n, p){
  n.tau <- length(tau)
  ICLamda0 <- matrix(NA, nrow = n, ncol = n.tau)
  
  res <- 0
  
  for(iterTau in 1:n.tau){
    
    if(tau[iterTau]>max.time){
      return(ICLamda0)
    }
    
    for(iterObs in 1:n){
      sum1 <- colMultiply_cpp(E1, scale = lambda0*(time1<=tau[iterTau]))
      sum2 <- lambda0/S01*(time1<=newT[iterObs])*(time1<=tau[iterTau])
      S0_tempo <- S01[prodlim::sindex(jump.times = time1, eval.times = newT[iterObs])]
      ICLamda0[iterObs,iterTau] <- - ICbeta[iterObs,,drop= FALSE] %*% colSums(sum1) - neweXb[iterObs] * sum(sum2) + (newT[iterObs]<=tau[iterTau]) * newStatus[iterObs]/S0_tempo
      res <- res +  (newT[iterObs]<=tau[iterTau]) * newStatus[iterObs]/S0_tempo
    }
    
  }
  
  return(ICLamda0)
}