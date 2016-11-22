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
#' d <- sampleData(1e2, outcome = "survival")[,.(eventtime,event,X1,X6)]
#' setkey(d, eventtime)
#' 
#' library(timereg)
#' mGS.cox <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d, resample.iid = TRUE)
#' ICbeta_GS <- mGS.cox$gamma.iid
#' IClambda_GS <- t(as.data.table(mGS.cox$B.iid))
#'  
#' m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE)
#' IC.cox <- iid.coxph(m.cox, center = FALSE)
#' 
#' range(IC.cox$ICbeta-ICbeta_GS)
#' range(IC.cox$ICLambda0/IClambda_GS[,NCOL(IClambda_GS)])
#' 
#' var(ICbeta_GS)
#' m.cox[["var"]]
#' var(ICbeta_GS)/mGS.cox$var.gamma
#' 
#' 
#' var(IC.cox$IC)/m.cox$var
#'
#' var(mGS.cox$gamma.iid)
#' m.cox$var
#' 
#' var(mGS.cox$gamma.iid)/m.cox$var
#' 
#' 
#'  
#' IC.cox$epsilon11/mGS.cox$gamma.iid
#' 

#' @rdname iid
#' @export
`iid` <- function(x,...) UseMethod("idd")

#' @rdname iid
#' @export
iid.coxph <- function(x, tau = NULL, ...){
  # to be adapted for stratified Cox models
  
  #### extract information from the Cox model
  resInfo <- getCoxInfo(x, design = FALSE, center = FALSE)
  data <- eval(x$call$data)
  eventtime <- x$y[,"time"]
  status <- x$y[,"status"]
  eXb <- exp(lpCox(x, data = data, xvars = resInfo$xvars, stratavars = resInfo$stratavars))
  X <- centerData(x, data = eval(x$call$data), stratavars = resInfo$stratavar, center = FALSE)
  
  if(is.null(tau)){tau <- max(eventtime)[1]}
  iInfo <- x$var
  
  # reorder data
  order <- order(eventtime)
  data <- data[order,,drop = FALSE]
  eventtime <- eventtime[order]
  status <- status[order]
  eXb <- eXb[order]
  X <- X[order,,drop = FALSE]
  
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
  
  #### Influence function for the beta
  # print(calcU(newT = eventtime, newX = X, newStatus = status, newN = NROW(X),
  #             eventtime = eventtime, eXb = eXb, X = X, p = NCOL(X), aggregate = TRUE))
  # print(calcU_cpp(newT = eventtime, newX = X, newStatus = status, newN = NROW(X),
  #       eventtime = eventtime, eXb = eXb, X = X, p = NCOL(X), aggregate = FALSE))
  ICbeta <- iidBeta(newT = eventtime, neweXb = eXb, newX = X, newStatus = status,
                    S01 = S0_U1times, E1 = E_U1times, indexTime1 = indexTime1, time1 = U1.time, 
                    n = n.obs, p = p.X, iInfo = iInfo)
  
 
  #### Influence function for the lambda0
  lambda0 <- predictCox(x, type = "hazard")
  
  ICLambda0 <- iidLambda0(tau = tau,
                          newT = eventtime, neweXb = eXb, newX = X, newStatus = status, ICbeta = ICbeta,
                          eventtime = eventtime, eXb = eXb, X = X, status = status, lambda0 = lambda0$hazard)
  
  #### Timereg iid
  iid_timereg <- list()
  Ueventtime <- unique(eventtime)
  n.Ueventtime <- length(Ueventtime)
  n.data <- NROW(data) 
  
  ## beta
  epsilon1 <- calcEpsilon1(x = x, data = data, eXb = eXb, X = X, status = status, time = eventtime)
  iid_timereg$beta <-  t(apply(epsilon1, 1, function(row){iInfo %*% row}))
  
  ## lambda
  iid_timereg$cumHazard <- matrix(NA, nrow = n.data, ncol = n.Ueventtime)
  
  for(iterT in 1:n.Ueventtime){
    epsilon3 <- calcEpsilon3(x = x, data = data, eXb = eXb, X = X, status = status, time = Ueventtime[iterT])
    #calcMartinRes(x =x ,data = data,status = status,time = eventtime[1],type = "hazard")/calcS0(t = eventtime[1], eventtime = eventtime, eXb = eXb)
    H <-  calcH(t = Ueventtime[iterT], eventtime = eventtime, eXb = eXb, X = X, status = status)
    #calcS1(t = eventtime[1], eventtime = eventtime, eXb = eXb, X = X)/calcS0(t = eventtime[1], eventtime = eventtime, eXb = eXb)^2
    
    epsilon2 <- epsilon3 + iid_timereg$beta %*% H
    iid_timereg$cumHazard[,iterT] <- epsilon2
  }
  
  #### export
  return(list(ICbeta = ICbeta, 
              ICLambda0 = ICLambda0, 
              iid_timereg = iid_timereg
  ))
  
}


iidBeta <- function(newT, neweXb, newX, newStatus,
                    S01, E1, indexTime1, time1, 
                    n, p, iInfo){
  
  #### Compute IC
  IC <- matrix(NA, nrow = n, ncol = p)
  
  for(iterI in 1:n){
    term1 <- calcU_cpp(newX = newX[iterI,,drop = FALSE], newStatus= newStatus[iterI], newN = 1,
                       IndexNewT = 0, ENewT = E1[prodlim::sindex(time1, newT[iterI]),,drop = FALSE]*newStatus[iterI],
                       p = p, aggregate = FALSE)
    
    term2 <- -sweep(sweep(E1[time1<=newT[iterI],,drop = FALSE], MARGIN = 2, FUN = "-", STATS = newX[iterI,,drop = FALSE]),
                    MARGIN = 1, FUN = "/", STATS = S01[time1<=newT[iterI]])
    
    IC[iterI,] <- as.double(iInfo %*% t(term1 - neweXb[iterI] * colSums(term2)))
  }
  
  return(IC)
}

iidLambda0 <- function(tau,
                       newT, neweXb, newX, newStatus, ICbeta,
                       eventtime, eXb, X, status, lambda0){
  
  n.time <- length(eventtime)
  n.obs <- NROW(X)
  p.X <- NCOL(X)
  
  #### Compute E
  res <- lapply(newT, FUN = function(t){
    return(calcE_cpp(t = t, eventtime = eventtime, eXb = eXb, X = X, n = n.obs, p = p.X))
  })
  
  S0_U1times <- unlist(lapply(res,"[[","S0"))
  S1_U1times <- matrix(unlist(lapply(res,"[[","S1")), nrow = n.time, ncol = p.X, byrow = TRUE)
  E_U1times <- matrix(unlist(lapply(res,"[[","E")), nrow = n.time, ncol = p.X, byrow = TRUE)
  
  #### Compute IC
  ICLamda0 <- rep(NA, n.obs)
  
  for(iterObs in 1:n.obs){
    sum1 <- sweep(E_U1times, MARGIN = 1, FUN = "*", STATS = lambda0*(eventtime<=tau))
    sum2 <- status*lambda0/S0_U1times*(eventtime<=newT[iterObs])*(eventtime<=tau)
    S0_tempo <- S0_U1times[prodlim::sindex(jump.times = eventtime, eval.times = newT[iterObs])]
    
    ICLamda0[iterObs] <- - ICbeta[iterObs,,drop= FALSE] %*% colSums(sum1) - neweXb[iterObs] * sum(sum2) + (newT[iterObs]<=tau) * newStatus[iterObs]/S0_tempo
  }
  
  return(ICLamda0)
}

calcS0 <- function(t, eventtime, eXb){
  res <- (t<=eventtime) * eXb 
  return(sum(res))
}
calcS1 <- function(t, eventtime, eXb, X){
  res <- sweep(X, MARGIN = 1, STATS = (t<=eventtime) * eXb, FUN = "*")
  return(colSums(res))
}
calcE <- function(t, eventtime, eXb, X){
  calcS1(t, eventtime, eXb, X)/calcS0(t, eventtime, eXb)
}
calcU <- function(newT, newX, newStatus, newN,
                  eventtime, eXb, X, 
                  p, E, aggregate){
  
  ## compute E
  index1 <- which(newStatus>0)
  U1.time <- unique(newT[index1])
  E_U1times <- t(sapply(U1.time, function(t){calcE(t = t,  eventtime = eventtime, eXb = eXb, X = X)}))
  
  ## compute the score
  score <- matrix(0, nrow = newN, ncol = p)
  
  for(iterPat in 1:newN){
    if(newStatus[iterPat]==0){next}
    score[iterPat,] <- newX[iterPat, ,drop = FALSE] - E_U1times[U1.time == newT[iterPat],,drop = FALSE]
  }
  if(aggregate){
    score <- colSums(score)
  }
  return(score)
}

calcS2 <- function(t, eventtime, eXb, X){
  vecH <- apply(X, 1, function(row){ tcrossprod(row)})
  res <- sweep(vecH, MARGIN = 2, STATS = (t<=eventtime) * eXb, FUN = "*")
  res <- rowSums(res)
  
  return(matrix(res, nrow = NCOL(X), ncol = NCOL(X)))
}

calcI <- function(t, eventtime, eXb, X, status, aggregate){
  indexEvent <- which(status==1)
  evaltime <- eventtime[indexEvent]
  n.evaltime <- sum(status)
  info <- array(0, dim = c(NCOL(X), NCOL(X), NROW(X)))
  
  for(iterT in 1:n.evaltime){
    Etempo <- calcE(t = evaltime[iterT], eventtime = eventtime, eXb, X)
    S0tempo <- calcS0(t = evaltime[iterT], eventtime = eventtime, eXb)
    S2tempo <- calcS2(t = evaltime[iterT], eventtime = eventtime, eXb, X)
    
    # handle ties
    indexPat <- which(eventtime[indexEvent] == evaltime[iterT])
    info[,,indexEvent[indexPat]] <- S2tempo/S0tempo-tcrossprod(Etempo)
  }
  
  if(aggregate){
    info <- apply(info,1:2,sum)
  }
  
  return(info)
}

calcU0 <- function(t, X, eventtime, status){
  # here t should be min(t, t_obs) where t is the time at which the IC is computed and t_obs is the event time of the observation
  status[eventtime>t] <- 0 # artificial censoring of the observations after t
  evaltime <- eventtime[status == 1]
  n.evaltime <- sum(status)
  
  score0 <- rep(0,NCOL(X))
  
  for(iterT in 1:n.evaltime){
    Etempo <- calcE(t = evaltime[iterT])
    S0tempo <- calcS0(t = evaltime[iterT])
    
    score0 <- score0 + (X - Etempo)/S0tempo
  }
  
  return(score0)
}

calcMartinRes <- function(x, data, status, time, type){ # dNi-dLambdai
  Utime <- unique(time)
  martinRes <- -predictCox(x, newdata = data, times = Utime, type = type, keep.times = FALSE)[[1]]
  
  for(iterI in 1:NROW(data)){
    
    martinRes[iterI,time[iterI]<unique(time)] <- 0 # no more at risk => Lambda = 0
    
    if(type == "cumHazard"){ # add dN
      martinRes[iterI,time[iterI]<=Utime] <- status[iterI] + martinRes[iterI,time[iterI]<=Utime] 
    }else if(type == "hazard"){
      martinRes[iterI,time[iterI]==Utime] <- status[iterI] + martinRes[iterI,time[iterI]==Utime]
    }
    
  }
  return(martinRes)
}

calcEpsilon1 <- function(x, data, eXb, X, status, time){
  
  epsilon1 <- matrix(0, nrow = NROW(X), ncol = NCOL(X))
  d_martinRes <- calcMartinRes(x = x, data = data, status = status, time = time, type = "hazard") # one value for each unique time and each patient
  
  Utime <- unique(time)
  for(iterT in 1:length(Utime)){
    Etempo <- calcE(t = Utime[iterT], eventtime = time, eXb = eXb, X = X) # one value for each patient
    epsilon1 <- epsilon1 + sweep(sweep(X, MARGIN = 2, STATS = Etempo, FUN = "-"), 
                                 MARGIN = 1, STATS = d_martinRes[,iterT], FUN = "*")
  }
  
  return(epsilon1)
}  

calcEpsilon3 <- function(x, data, eXb, X, status, time){
  
  epsilon3 <- matrix(0, nrow = NROW(X), ncol = 1)
  d_martinRes <- calcMartinRes(x = x, data = data, status = status, time = time, type = "hazard") # one value for each unique time and each patient
  
  Utime <- unique(time)
  for(iterT in 1:length(Utime)){
    S0tempo <- calcS0(t = Utime[iterT], eventtime = time, eXb = eXb) # one value for each patient
    epsilon3 <- epsilon3 + d_martinRes[,iterT]/S0tempo
  }
  
  return(epsilon3)
}

calcH <- function(t, eventtime, eXb, X, status){
  
  H <- rep(0, times = NCOL(X))
  
  Utime <- unique(eventtime[eventtime<=t])
  
  for(iterT in 1:length(Utime)){
    S1tempo <- calcS1(Utime[iterT], eventtime, eXb, X)
    S0tempo <- calcS0(Utime[iterT], eventtime, eXb)
    
    # handle ties
    indexPat <- which(eventtime == eventtime[iterT])
    H <- H - sum(status[indexPat])*S1tempo/S0tempo^2
  }
  
  return(H)
} 





# iid2.coxph <- function(x, newdata = eval(x$call$data), time = NULL, ...){
# 
#   #### extract information from the Cox model
#   
#   ## training
#   resInfo <- getCoxInfo(x, design = FALSE)
#   
#   olddata <- eval(x$call$data)
#   oldtime <- x$y[,"time"]
#   oldstatus <- x$y[,"status"]
#   oldeXb <- exp(lpCox(x, data = olddata, xvars = resInfo$xvars, stratavars = resInfo$stratavars))
#   oldX <- centerData(x, data = eval(x$call$data), stratavars = resInfo$stratavar)
#    
#   # reorder data
#   oldorder <- order(oldtime)
#   olddata <- olddata[oldorder,,drop = FALSE]
#   oldtime <- oldtime[oldorder]
#   oldstatus <- oldstatus[oldorder]
#   oldeXb <- oldeXb[oldorder]
#   oldX <- oldX[oldorder,,drop = FALSE]
#   
#   ## test
#   f.surv <- x$call$formula[[2]]
#   if("time" %in% names(f.surv) == FALSE && "event" %in% names(f.surv) == FALSE ){
#     names(f.surv)[2:3] <- c("time","event")
#   }else if("time" %in% names(f.surv) == FALSE){
#     names(f.surv)[setdiff(2:3,which(names(f.surv)=="event"))] <- "time"
#   }else if("event" %in% names(f.surv) == FALSE){
#     names(f.surv)[setdiff(2:3,which(names(f.surv)=="time"))] <- "event"
#   }
#   timeVar <- all.vars(f.surv[[which(names(f.surv)=="time")]])
#   statusVar <- all.vars(f.surv[[which(names(f.surv)=="event")]])
#   
#   newtime <- newdata[[timeVar]]
#   newstatus <- newdata[[statusVar]]
#   neweXb <- exp(lpCox(x, data = newdata, xvars = resInfo$xvars, stratavars = resInfo$stratavars))
#   newX <- centerData(x, data = newdata, stratavars = resInfo$stratavar)
#   
#   # reorder data
#   neworder <- order(newtime)
#   newdata <- newdata[neworder,,drop = FALSE]
#   newtime <- newtime[neworder]
#   newstatus <- newstatus[neworder]
#   neweXb <- neweXb[neworder]
#   newX <- newX[neworder,,drop = FALSE]
#   
#   if(is.null(time)){
#     time <- tail(oldtime,1)
#   }
# 

