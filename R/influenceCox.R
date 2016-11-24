#' @title Extract i.i.d. decomposition from a Cox model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iid
#' 
#' @param x object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
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
#' IC.cox <- iid(m.cox)
#' 
#' range(IC.cox$ICbeta-ICbeta_GS)
#' crossprod(ICbeta_GS)
#' m.cox$var
#' mGS.cox$robvar.gamma
#' crossprod(IC.cox$ICbeta)
#'
#' range(IC.cox$ICLambda0/IClambda_GS[,NCOL(IClambda_GS)])
#' range(IC.cox$ICLambda0 - IC.cox$iid_timereg$cumHazard[,NCOL(IC.cox$iid_timereg$cumHazard)])
#' 
#' range(IC.cox$iid_timereg$cumHazard / IClambda_GS[,-1])
#' 
#' IC.cox <- iid(m.cox, tauLambda = 7)
#'  
#' IC.cox$iid_timereg$cumHazard /  IClambda_GS[,-1]
#'  
#' sweep(IC.cox$iid_timereg$cumHazard,IC.cox$ICLambda0,FUN = "/",MARGIN = 1)
#' 

#' @rdname iid
#' @export
`iid` <- function(x, tau, ...) UseMethod("iid")

#' @rdname iid
#' @method iid coxph
#' @export
iid.coxph <- function(x, newdata = NULL, tauLambda = NULL, tauBeta = NULL, keep.times = TRUE, add.timeRegIID = FALSE, ...){
  # to be adapted for stratified Cox models
  
  #### extract information from the Cox model
  resInfo <- getCoxInfo(x, design = FALSE)
  data <- eval(x$call$data)
  eventtime <- x$y[,"time"]
  status <- x$y[,"status"]
  eXb <- exp(getLP(x, data = data, lpVars = resInfo$lpvars, stratavars = resInfo$stratavars, center = FALSE))
  X <- as.matrix(getData(x, data = eval(x$call$data), lpVars = resInfo$lpvars, stratavars = resInfo$stratavar, center = FALSE))
  lastEventtime <- eventtime[length(eventtime)]
  
  iInfo <- x$var
  
  # reorder data
  order <- order(eventtime)
  data <- data[order,,drop = FALSE]
  eventtime <- eventtime[order]
  status <- status[order]
  eXb <- eXb[order]
  X <- X[order,,drop = FALSE]
  
  #### Extract new observations
  if(is.null(newdata)){
    newT <- eventtime
    newStatus <- status
    neweXb <- eXb
    newX <- X
  }else{
    newT <- newdata[[resInfo$timevars]]
    newStatus <- newdata[[resInfo$eventvars]]
    neweXb <- exp(getLP(x, data = newdata, lpVars = resInfo$lpvars, stratavars = resInfo$stratavars, center = FALSE))
    newX <- as.matrix(getData(x, data = newdata, lpVars = resInfo$lpvars, stratavars = resInfo$stratavars, center = FALSE))
  }
  
  # reorder data
  neworder <- order(newT)
  newT <- newT[neworder]
  newStatus <- newStatus[neworder]
  neweXb <- neweXb[neworder]
  newX <- newX[neworder,,drop = FALSE]
  
  #### Compute quantities of interest only for event times
  n.obs <- NROW(X)
  p.X <- NCOL(X)
  
  U1.time <- c(unique(eventtime[status>0]),tail(eventtime,1)+1e-12)
  nU1.time <- length(U1.time)
  indexTime1 <- rep(nU1.time, n.obs)
  indexTime1[status>0] <- match(eventtime[status>0], U1.time)-1
  
  resE <- lapply(U1.time, FUN = function(t){
    return(calcE_cpp(t = t, eventtime = eventtime, eXb = eXb, X = X, n = n.obs, p = p.X))
  })

  S0_U1times <- unlist(lapply(resE,"[[","S0"))
  S1_U1times <- matrix(unlist(lapply(resE,"[[","S1")), nrow = nU1.time, ncol = p.X, byrow = TRUE)
  E_U1times <- matrix(unlist(lapply(resE,"[[","E")), nrow = nU1.time, ncol = p.X, byrow = TRUE)
  
  if(is.null(tauLambda)){tauLambda <- U1.time[-nU1.time]}
  if(is.null(tauBeta)){tauBeta <- U1.time[nU1.time-1]}
  
  #### Influence function for the beta
  ICbeta <- iidBeta(tau = tauBeta, newT = newT, neweXb = neweXb, newX = newX, newStatus = newStatus,
                    S01 = S0_U1times, E1 = E_U1times, indexTime1 = indexTime1, time1 = U1.time, 
                    n = n.obs, p = p.X, iInfo = iInfo)
  if(keep.times){attr(ICbeta,"time") <- tauBeta}
  
  #### Influence function for the baseline hazard
  # this function is slow for the moment
  lambda0 <- predictCox(x, type = "hazard")
  
  ICLambda0 <- iidLambda0(tau = tauLambda, max.time = lastEventtime,
                          newT = newT, neweXb = neweXb, newX = newX, newStatus = newStatus, ICbeta = ICbeta,
                          S01 = S0_U1times[-nU1.time], E1 = E_U1times[-nU1.time,,drop = FALSE], indexTime1 = indexTime1[-nU1.time], time1 = U1.time[-nU1.time], 
                          lambda0 = lambda0[match(U1.time[-nU1.time],lambda0[,"time"]),"hazard"],
                          n = n.obs, p = p.X)
  if(keep.times){attr(ICLambda0,"time") <- tauLambda}
  
  ## rescale
  ICLambda0 <- ICLambda0 * exp(-as.double(coef(x) %*% x$means))
  
  #### Timereg iid
  iid_timereg <- list()
  if(add.timeRegIID){
    
    ## beta
    epsilon1 <- calcEpsilon1(x = x, data = data, eXb = eXb, X = X, status = status, time = eventtime, n = n.obs, p = p.X)
    iid_timereg$beta <-  t(apply(epsilon1, 1, function(row){iInfo %*% row}))
    
    ## lambda
    epsilon2 <- matrix(0, nrow = n.obs, ncol = nU1.time)
    
    for(iterT in 1:(nU1.time-1)){
      Htempo <- S1_U1times[iterT,,drop = FALSE]/S0_U1times[iterT]^2
      #              event  at jump                                                                                                                        still at risk
      dmartinRes <- status*(eventtime==U1.time[iterT]) - predictCox(x, newdata = data, times = U1.time[iterT], type = "hazard", keep.times = FALSE)[[1]]*(eventtime>=U1.time[iterT])
      epsilon3 <- dmartinRes[,1]/S0_U1times[iterT]
      
      epsilon2[,iterT+1] <- epsilon2[,iterT] + epsilon3 - iid_timereg$beta %*% t(Htempo)
    }
    iid_timereg$cumHazard <- epsilon2[,-1,drop = FALSE]# cbind(epsilon2,0)[,prodlim::sindex(jump.times = c(0,U1.time), eval.times = unique(eventtime))]
    
    attr(iid_timereg$cumHazard,"cumHazard") <- U1.time[-length(U1.time)]
  }
  
  #### export
  return(list(ICbeta = ICbeta[order[neworder],,drop=FALSE],  # restaure original ordering
              ICLambda0 = ICLambda0[order[neworder],,drop=FALSE],
              iid_timereg = iid_timereg
  ))
}


iidBeta <- function(tau, newT, neweXb, newX, newStatus,
                    S01, E1, indexTime1, time1, 
                    n, p, iInfo){
  
  IC <- matrix(NA, nrow = n, ncol = p)
  
  for(iterI in 1:n){
    if(newT[iterI]<=tau){
      term1 <- calcU_cpp(newX = newX[iterI,,drop = FALSE], newStatus= newStatus[iterI], newN = 1,
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
                       S01, E1, indexTime1, time1, lambda0,
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


calcEpsilon1 <- function(x, data, eXb, X, status, time, n, p){
  
  epsilon1 <- matrix(0, nrow = NROW(X), ncol = NCOL(X))
  d_martinRes <- calcMartinRes(x = x, data = data, status = status, time = time, type = "hazard") # one value for each unique time and each patient
  
  Utime <- unique(time)
  for(iterT in 1:length(Utime)){
    Etempo <- calcE_cpp(t = Utime[iterT], eventtime = time, eXb = eXb, X = X, n = n, p = p)$E # one value for each patient
    epsilon1 <- epsilon1 + sweep(sweep(X, MARGIN = 2, STATS = Etempo, FUN = "-"), 
                                 MARGIN = 1, STATS = d_martinRes[,iterT], FUN = "*")
  }
  
  return(epsilon1)
}  

calcMartinRes <- function(x, data, status, time, type, browser = FALSE){ # dNi-dLambdai
  Utime <- unique(time)
  martinRes <- -predictCox(x, newdata = data, times = Utime, type = type, keep.times = FALSE)[[1]]
  for(iterI in 1:NROW(data)){
    martinRes[iterI,time[iterI]<unique(time)] <- 0 # no more at risk => Lambda = 0
    if(type == "cumHazard"){ # add dN
      martinRes[iterI,time[iterI]<=Utime] <- status[iterI] + martinRes[iterI,time[iterI]<=Utime] 
    }else if(type == "hazard"){
      if(browser) browser()
      martinRes[iterI,time[iterI]==Utime] <- status[iterI] + martinRes[iterI,time[iterI]==Utime]
      martinRes[iterI,] <- martinRes[iterI,]*(time<=data$eventtime[iterI]) # only residual when at risk
    }
  }
  
  return(martinRes)
}