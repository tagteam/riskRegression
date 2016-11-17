`iid` <- function(x,...) UseMethod("idd")

#' @examples
#' library(survival)
#' set.seed(10)
#' d <- sampleData(2e1, outcome = "survival")[,.(eventtime,event,X1,X6)]
#' setkey(d, eventtime)
#' d[,event := 1]
#' 
#' library(timereg)
#' mGS.cox <- cox.aalen(Surv(eventtime, event) ~ prop(X1)+prop(X6), data = d, resample.iid = TRUE)
#'  
#' m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE)
#' IC.cox <- iid.coxph(m.cox, center = FALSE)
#' var(IC.cox$IC)
#' m.cox$var
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
#' 
iid.coxph <- function(x, ...){
  
  #### extract information from the Cox model
  resInfo <- getCoxInfo(x, design = FALSE)
  data <- eval(x$call$data)
  eventtime <- x$y[,"time"]
  status <- x$y[,"status"]
  eXb <- exp(lpCox(x, data = data, xvars = resInfo$xvars, stratavars = resInfo$stratavars))
  X <- centerData(x, data = eval(x$call$data), stratavars = resInfo$stratavar, center = FALSE)
  
  tau <- max(eventtime)[1]
  iInfo <- x$var
  
  # reorder data
  order <- order(eventtime)
  data <- data[order,,drop = FALSE]
  eventtime <- eventtime[order]
  status <- status[order]
  eXb <- eXb[order]
  X <- X[order,,drop = FALSE]
  
  #### Compute S0, S1, S2, U, E
  calcS0(t = tau, eventtime, eXb)
  calcS1(t = tau, eventtime, eXb, X)
  calcE(t = tau, eventtime, eXb, X)
  ScoreLv <- calcU(t = tau, eventtime = eventtime, eXb = eXb, X = X, status = status, aggregate = TRUE)
  cat("Sample average of the individual scores \n")
  print(ScoreLv)
  InfoLv <- calcI(t = tau, eventtime = eventtime, eXb = eXb, X = X, status = status, aggregate = TRUE)
  cat("Difference between the second derivative of the loglikelihood and the inverse of the Variance Covariance matrix \n")
  print(InfoLv-solve(x$var))
  
  martinRes <- calcMartinRes(x = x, data = data, status = status, time = eventtime, type = "cumHazard")
  cat("Difference between the Martingal residuals from RR and survival \n")
  print(range(residuals(x, type = "martingale")[order]-diag(martinRes)))
  d_martinRes <- calcMartinRes(x = x, data = data, status = status, time = eventtime, type = "hazard")
  
  # fields:::image.plot(martinRes)
  # fields:::image.plot(d_martinRes)
  
  #### Influence function
  ICobs <- calcIC1(t = tau, eventtime = eventtime, eXb = eXb, X = X, status = status, iInfo = iInfo)
  # Sigma_IC1 <- var(ICobs)
  
  # compute the martingale residuals and apply the formula of torben  
  epsilon1 <- calcEpsilon1(x = x, data = data, eXb = eXb, X = X, status = status, time = eventtime)
  
  #Var_epsilon1 <- matrix(apply(apply(epsilon1, 1, tcrossprod), 1, sum),2,2)
  #Sigma_Epsilon1 <- NROW(data) * iInfo %*% Var_epsilon1 %*% iInfo
  
  return(list(IC = ICobs, 
              epsilon1 = epsilon1,
              epsilon11 = t(apply(epsilon1, 1, function(row){iInfo %*% row}))
              ))
  
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
calcU <- function(t, eventtime, eXb, X, status, aggregate){
  status[eventtime>t] <- 0 # artificial censoring of the observations after t
  indexEvent <- which(status==1)
  evalX <- X[indexEvent,,drop = FALSE]
  evaltime <- unique(eventtime[indexEvent])
  n.evaltime <- length(evaltime)
  
  score <- matrix(0, nrow = NROW(X), ncol = NCOL(X))
  
  for(iterT in 1:n.evaltime){
    Etempo <- calcE(t = evaltime[iterT], eventtime, eXb, X)
    # handle ties
    indexPat <- which(eventtime[indexEvent] == evaltime[iterT])
    score[indexEvent[indexPat],] <- sweep(evalX[indexPat, ,drop = FALSE], MARGIN = 2, STATS = Etempo, FUN = "-") 
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

fctU0 <- function(t, X, eventtime =  oldtime, status = oldstatus){
  # here t should be min(t, t_obs) where t is the time at which the IC is computed and t_obs is the event time of the observation
  status[eventtime>t] <- 0 # artificial censoring of the observations after t
  evaltime <- eventtime[status == 1]
  n.evaltime <- sum(status)
  
  score0 <- rep(0,NCOL(X))
  
  for(iterT in 1:n.evaltime){
    Etempo <- fctE(t = evaltime[iterT])
    S0tempo <- fctS0(t = evaltime[iterT])
    
    score0 <- score0 + (X - Etempo)/S0tempo
  }
  
  return(score0)
}

calcMartinRes <- function(x, data, status, time, type){ # dNi-dLambdai
  Utime <- unique(time)
  martinRes <- -predictCox(x, newdata = data, times = Utime, type = type, keep.times = FALSE)[[1]]
  
  for(iterI in 1:NROW(data)){
    
    martinRes[iterI,time[iterI]<unique(time)] <- 0 # no more at risk => Lambda = 0
    
    if(type == "cumHazard"){
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
  for(iterT in 1:NCOL(d_martinRes)){
    Etempo <- calcE(t = Utime[iterT], eventtime = time, eXb = eXb, X = X) # one value for each patient
    epsilon1 <- epsilon1 + sweep(sweep(X, MARGIN = 2, STATS = Etempo, FUN = "-"), 
                                 MARGIN = 1, STATS = d_martinRes[,iterT], FUN = "*")
  }

  return(epsilon1)
}  

calcIC1 <- function(t, eventtime, eXb, X, status, iInfo, version = "1"){
  Score <- calcU(t = t, X = X, eXb = eXb, eventtime = eventtime, status = status, aggregate = FALSE)
  E <- calcE(t = t, eventtime = eventtime, eXb = eXb, X = X)
  E_alltimes <- t(sapply(eventtime, FUN = function(tt){calcE(tt, eventtime = eventtime, eXb = eXb, X = X)}))
  
  S0 <- calcS0(t = t, eventtime = eventtime, eXb = eXb)
  S0_alltimes <- sapply(eventtime, FUN = calcS0, eventtime = eventtime, eXb = eXb)
  
  IC <- matrix(NA, nrow = NROW(X), ncol = NCOL(X))
  
  for(iterI in 1:NROW(X)){
    Xtempo <- X[iterI,,drop=FALSE]
    
    term1 <- Score[iterI,,drop=FALSE]
    term2 <- sweep(sweep(E_alltimes, MARGIN = 2, FUN = "-", STATS = Xtempo),
                   MARGIN = 1, FUN = "*", STATS = (eventtime <= eventtime[iterI]) * status / S0_alltimes)
    
    
    IC[iterI,] <- iInfo %*% t(term1 - eXb[iterI] * colSums(term2))
  }
  
  return(IC)
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
#   #### Compute S0, S1, S2, U, E
#   fctS0 <- function(t, eXb = oldeXb, eventtime = oldtime){
#     res <- (t<=eventtime) * eXb 
#     return(sum(res))
#   }
#   fctS1 <- function(t, eXb = oldeXb, X = oldX, eventtime = oldtime){
#     res <- sweep(X, MARGIN = 1, STATS = (t<=eventtime) * eXb, FUN = "*")
#     return(colSums(res))
#   }
#   fctE <- function(t, eXb = oldeXb, X = oldX, eventtime = oldtime){
#     fctS1(t)/fctS0(t)
#   }
#   fctU <- function(t, X, eventtime, status, aggregate){
#     status[eventtime>t] <- 0 # artificial censoring of the observations after t
#     indexEvent <- which(status==1)
#     evalX <- X[indexEvent,,drop = FALSE]
#     evaltime <- unique(eventtime[indexEvent])
#     n.evaltime <- length(evaltime)
#     
#     score <- matrix(0, nrow = NROW(X), ncol = NCOL(X))
#     
#     for(iterT in 1:n.evaltime){
#       Etempo <- fctE(t = evaltime[iterT])
#       # handle ties
#       indexPat <- which(eventtime[indexEvent] == evaltime[iterT])
#       score[indexEvent[indexPat],] <- sweep(evalX[indexPat, ,drop = FALSE], MARGIN = 2, STATS = Etempo, FUN = "-") 
#     }
#     
#     if(aggregate){
#       score <- colSums(score)
#     }
#     return(score)
#   }
#   fctU0 <- function(t, X, eventtime =  oldtime, status = oldstatus){
#     # here t should be min(t, t_obs) where t is the time at which the IC is computed and t_obs is the event time of the observation
#     status[eventtime>t] <- 0 # artificial censoring of the observations after t
#     evaltime <- eventtime[status == 1]
#     n.evaltime <- sum(status)
#     
#     score0 <- rep(0,NCOL(X))
#     
#     for(iterT in 1:n.evaltime){
#       Etempo <- fctE(t = evaltime[iterT])
#       S0tempo <- fctS0(t = evaltime[iterT])
#       
#       score0 <- score0 + (X - Etempo)/S0tempo
#     }
#     
#     return(score0)
#   }
#   
#   fctS2 <- function(t, eXb, X, eventtime){
#     vecH <- apply(X, 1, function(row){ tcrossprod(row)})
#     res <- sweep(vecH, MARGIN = 2, STATS = (t<=eventtime) * eXb, FUN = "*")
#     res <- rowSums(res)
#     
#     return(matrix(res, nrow = NCOL(X), ncol = NCOL(X)))
#   }
#   fctI <- function(t, eXb, X, eventtime, status, aggregate){
#     indexEvent <- which(status==1)
#     evaltime <- eventtime[indexEvent]
#     n.evaltime <- sum(status)
#     info <- array(0, dim = c(NCOL(X), NCOL(X), NROW(X)))
#     
#     for(iterT in 1:n.evaltime){
#       Etempo <- fctE(t = evaltime[iterT], eXb = eXb, X = X, eventtime = eventtime)
#       S0tempo <- fctS0(t = evaltime[iterT], eXb = eXb, eventtime = eventtime)
#       S2tempo <- fctS2(t = evaltime[iterT], eXb = eXb, X = X, eventtime = eventtime)
#       
#       # handle ties
#       indexPat <- which(eventtime[indexEvent] == evaltime[iterT])
#       info[,,indexEvent[indexPat]] <- S2tempo/S0tempo-tcrossprod(Etempo)
#     }
#     
#     if(aggregate){
#        info <- apply(info,1:2,sum)
#     }
#     
#     return(info)
#   }
#   
#   fctS0(t = time)
#   fctS1(t = time)
#   fctE(t = time)
#   ScoreLv <- fctU(t = time, X = oldX, eventtime = oldtime, status = oldstatus, aggregate = TRUE)
#   print(ScoreLv)
#   InfoLv <- fctI(t = time, eXb = oldeXb, X = oldX, eventtime = oldtime, status = oldstatus, aggregate = TRUE)
#   print(InfoLv-solve(x$var))
#   
#   
#   martingalRes <- -predictCox(x, newdata = olddata, times = unique(oldtime), type = "cumHazard", keep.times = FALSE)[[1]]
#   dmartingalRes <- -predictCox(x, newdata = olddata, times = unique(oldtime), type = "hazard", keep.times = FALSE)[[1]]
#   for(iterI in 1:NROW(olddata)){
#     martingalRes[iterI,oldtime[iterI]<unique(oldtime)] <- 0
#     martingalRes[iterI,oldtime[iterI]<=unique(oldtime)] <- oldstatus[iterI] + martingalRes[iterI,oldtime[iterI]<=unique(oldtime)]
#     dmartingalRes[iterI,oldtime[iterI]<unique(oldtime)] <- 0
#     dmartingalRes[iterI,oldtime[iterI]==unique(oldtime)] <- oldstatus[iterI] + dmartingalRes[iterI,oldtime[iterI]==unique(oldtime)]
#   }
#   print(range(residuals(x, type = "martingale")[neworder]-diag(martingalRes)))
#   
#   # fields:::image.plot(martingalRes)
#   # fields:::image.plot(dmartingalRes[1:15,1:10], x = 1:15, y = 1:10)
#   
#   epsilon1 <- matrix(0, nrow = NROW(newdata), ncol = NCOL(newX))
#   for(iterT in 1:NCOL(dmartingalRes)){
#     epsilon1 <- epsilon1 + sweep(sweep(oldX, MARGIN = 2, STATS = fctE(t = unique(oldtime)[iterT]), FUN = "-"), MARGIN = 1, STATS = dmartingalRes[,iterT], FUN = "*")
#   }
#   browser()
#   
#   #### Influence function
#   Info <- x$var
#   iInfo <- solve(Info)
#   Score <- fctU(t = time, X = newX, eventtime = newtime, status = newstatus, aggregate = FALSE)
#   fctI(t = time, eXb = oldeXb, X = oldX, eventtime = oldtime, status = oldstatus, aggregate = FALSE)[,,1]
#   
#   fctIC <- function(index){
#     iInfo %*% t(Score[index,,drop=TRUE] - neweXb[index] * fctU0(t = min(time, newtime[index]), X = newX[index,,drop = FALSE]))
#   }
#   
#   IC <- sapply(1:NROW(newX),fctIC)
#     var(t(IC))
#                
#   # compute the martingale residuals and apply the formula of torben  
#     Veps <- matrix(apply(apply(epsilon1, 1, tcrossprod), 1, sum),2,2)
#     
#       
#     iInfo %*% Veps %*% iInfo
#     
#   
# 
# }
# 
