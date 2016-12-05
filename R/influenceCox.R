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
#' IC.cox <- iidCox(m.cox)
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
  
  # reorder data
  object.order <- order(object.time)
  object.time <- object.time[object.order]
  object.status <- object.status[object.order]
  object.eXb <- object.eXb[object.order]
  object.LPdata <- object.LPdata[object.order,,drop = FALSE]
  object.indexStrata <- lapply(object.levelStrata, function(level){which(object.strata==level)-1})
  
  if(is.null(tauLambda)){tauLambda <- object.time[object.status>0]}
  if(is.null(tauBeta)){tauBeta <- tail(object.time[object.status>0],1)}
  if(any(tauLambda > tauBeta)){
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
    new.order <- order(new.time)
    new.time <- new.time[new.order]
    new.strata <- new.strata[new.order]
    new.status <- new.status[new.order]
    new.eXb <- new.eXb[new.order]
    new.LPdata <- new.LPdata[new.order,,drop = FALSE]
    new.indexStrata <-lapply(object.levelStrata, function(level){which(new.strata==level)-1})
    
    # for factor variables: retain only the relevant column from the design matrix, 
    new.LPdata <- new.LPdata[, object.designVar, drop = FALSE]
  }
  
  
  #### Compute quantities of interest only for event times
  lastEventtime <- object.time[length(object.time)]
  p <- NCOL(object.LPdata)
  nStrata <- length(levels(object.strata))
  
  ## baseline hazard
  lambda0 <- predictCox(object, type = "hazard", centered = center.lambda0, keep.strata = TRUE)
  
  ### v1
  # Ecpp <- calcEstrata_cpp(tau = tauBeta,
  #                         indexStrata = object.indexStrata,
  #                         status = object.status,
  #                         p = p, nStrata = nStrata,
  #                         eventtime = object.time,
  #                         eXb = object.eXb,
  #                         X = object.LPdata,
  #                         add0 = TRUE)
  
  ### v2
  U1.time <- c(unique(object.time[object.status>0]),tail(object.time,1)+1e-12)
  nU1.time <- length(U1.time)
  indexTime1 <- rep(nU1.time, object.n)
  indexTime1[object.status>0] <- match(object.time[object.status>0], U1.time)-1

  resE <- lapply(U1.time, FUN = function(t){
    return(calcE_cpp(t = t, eventtime = object.time, eXb = object.eXb, X = object.LPdata, n = object.n, p = p))
  })

  S0_U1times <- unlist(lapply(resE,"[[","S0"))
  S1_U1times <- matrix(unlist(lapply(resE,"[[","S1")), nrow = nU1.time, ncol = p, byrow = TRUE)
  if(NCOL(S1_U1times)==0){S1_U1times <- cbind(S1_U1times,0)}
  E_U1times <- matrix(unlist(lapply(resE,"[[","E")), nrow = nU1.time, ncol = p, byrow = TRUE)
  if(NCOL(E_U1times)==0){E_U1times <- cbind(E_U1times,0)}
  
  #### Computation of the influence function
  ICbeta <- iidBeta(tau = tauBeta,
                    newT = new.time,
                    neweXb = new.eXb,
                    newX = new.LPdata,
                    newStatus = new.status,
                    S01 = S0_U1times,
                    E1 = E_U1times,
                    time1 = U1.time,
                    n = new.n, p = p, iInfo = iInfo)

  
  ICLambda0 <- iidLambda0(tau = tauLambda, max.time = lastEventtime,
                          newT = new.time,
                          neweXb = new.eXb,
                          newX = new.LPdata,
                          newStatus = new.status,
                          ICbeta = ICbeta,
                          S01 = S0_U1times[-nU1.time],
                          E1 = E_U1times[-nU1.time,,drop = FALSE],
                          time1 = U1.time[-nU1.time],
                          lambda0 = lambda0[match(U1.time[-nU1.time],lambda0[,"time"]),"hazard"],
                          n = object.n, p = p)

  if(keep.times){
    attr(ICbeta,"time") <- tauBeta
    attr(ICLambda0,"time") <- tauLambda
  }
  
  
  # ls <- list(tau = tauLambda, max.time = lastEventtime,
  #            newT = new.time, neweXb = new.eXb, newX = new.LPdata, newStatus = new.status, ICbeta = ICbeta,
  #            S01 = S0_U1times[-nU1.time], E1 = E_U1times[-nU1.time,,drop = FALSE], indexTime1 = indexTime1[-nU1.time], time1 = U1.time[-nU1.time], 
  #            lambda0 = lambda0[match(U1.time[-nU1.time],lambda0[,"time"]),"hazard"],
  #            n = object.n, p = p)
  
  
  ## rescale
  if(center.result == TRUE && !is.null(CoxCenter(object))){
    ICLambda0 <- ICLambda0 * exp(-as.double(coef(object) %*% CoxCenter(object)))
  }
  
  #### export
  return(list(ICbeta = ICbeta[new.order,,drop=FALSE],  # restaure original ordering
              ICLambda0 = ICLambda0[new.order,,drop=FALSE]
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
    
    term2 <- -sweep(sweep(E1[time1<=min(newT[iterI],tau),,drop = FALSE], MARGIN = 2, FUN = "-", STATS = newX[iterI,,drop = FALSE]),
                    MARGIN = 1, FUN = "/", STATS = S01[time1<=min(newT[iterI],tau)]) 
    
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
  
  for(iterTau in 1:n.tau){
    
    if(tau[iterTau]>max.time){
      return(ICLamda0)
    }
    
    for(iterObs in 1:n){
      sum1 <- sweep(E1, MARGIN = 1, FUN = "*", STATS = lambda0*(time1<=tau[iterTau]))
      sum2 <- lambda0/S01*(time1<=newT[iterObs])*(time1<=tau[iterTau])
      S0_tempo <- S01[prodlim::sindex(jump.times = time1, eval.times = newT[iterObs])]
      ICLamda0[iterObs,iterTau] <- - ICbeta[iterObs,,drop= FALSE] %*% colSums(sum1) - neweXb[iterObs] * sum(sum2) + (newT[iterObs]<=tau[iterTau]) * newStatus[iterObs]/S0_tempo
    }
    
  }
  
  return(ICLamda0)
}