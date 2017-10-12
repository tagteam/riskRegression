### influence2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 19 2017 (16:01) 
## Version: 
## last-updated: Oct 12 2017 (17:02) 
##           By: Thomas Alexander Gerds
##     Update #: 119
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ doc
#' @title Influence test [Experimental!!]
#' @description Compare two estimates using their influence function
#' @name influenceCoxTest
#' 
#' @param object either a list of models or an object of class predictCox or predictCSC.
#' @param object2 same as predict1 but for another model.
#' @param type the type of predicted value. 
#' @param nsim.band the number of simulations used to compute the quantiles
#' for the confidence bands.
#' @param newdata A \code{data.frame} or \code{data.table} containing
#'     the values of the predictor variables defining subject specific
#'     predictions.
#' @param times Time points at which to evaluate the predictions.
#' @param cause Identifies the cause of interest among the competing
#'     events.
#' @param conf.level Level of confidence.
#' @param tanh.transform Should the comparison be made on the arc-tangent scale or on the original scale ?
#' @param ... additional arguments to be passed to lower level functions.
#' 
#' @examples
#' library(lava)
#' library(survival)
#' n <- 100
#'
#' #### Under H1
#' set.seed(1)
#' newdata <- data.frame(X1=1)
#'
#' ## simulate non proportional hazard using lava
#' m <- lvm()
#' regression(m) <- y ~ 1
#' regression(m) <- s ~ exp(-2*X1)
#' distribution(m,~X1) <- binomial.lvm()
#' distribution(m,~cens) <- coxWeibull.lvm(scale=1)
#' distribution(m,~y) <- coxWeibull.lvm(scale=1,shape=~s)
#' eventTime(m) <- eventtime ~ min(y=1,cens=0)
#' d <- as.data.table(sim(m,n))
#' setkey(d, eventtime)
#'
#' ## fit cox models
#' m.cox <- coxph(Surv(eventtime, status) ~ X1, 
#'                data = d, y = TRUE, x = TRUE)
#'
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), 
#'                      data = d, y = TRUE, x = TRUE)
#'
#' ## compare models
#' # one time point
#' influenceCoxTest(list(m.cox, mStrata.cox), tanh.transform = FALSE,
#'                  type = "survival", newdata = newdata, times = 0.5)
#' influenceCoxTest(list(m.cox, mStrata.cox), tanh.transform = TRUE,
#'                  type = "survival", newdata = newdata, times = 0.5)
#'                                  
#' # several timepoints
#' influenceCoxTest(list(m.cox, mStrata.cox), tanh.transform = TRUE,
#'                  type = "survival", newdata = newdata, times = seq(0.1,1,by=0.1))
#' influenceCoxTest(list(m.cox, mStrata.cox), tanh.transform = FALSE,
#'                  type = "survival", newdata = newdata, times = seq(0.1,1,by=0.1))
#'
#' #### Under H0 (Cox) ####
#' set.seed(1)
#' ## simulate proportional hazard using lava
#' m <- lvm()
#' regression(m) <- y ~ 1
#' distribution(m,~X1) <- binomial.lvm()
#' distribution(m,~cens) <- coxWeibull.lvm()
#' distribution(m,~y) <- coxWeibull.lvm()
#' eventTime(m) <- eventtime ~ min(y=1,cens=0)
#' d <- as.data.table(sim(m,n))
#' setkey(d, eventtime)
#' 
#' ## fit cox models
#' Utime <- sort(unique(d$eventtime))
#' m.cox <- coxph(Surv(eventtime, status) ~ X1, 
#'                data = d, y = TRUE, x = TRUE)
#'
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), 
#'                      data = d, y = TRUE, x = TRUE)
#'
#' p.cox <- predictCox(m.cox, newdata = newdata, time = Utime, type = "survival")
#' p.coxStrata <- predictCox(mStrata.cox, newdata = newdata, time = Utime, type = "survival")
#' autoplot(p.cox)
#' autoplot(p.coxStrata)
#'  
#' ## compare models
#' influenceCoxTest(list(m.cox, mStrata.cox), tanh.transform = TRUE,
#'                  type = "survival", newdata = newdata, times = Utime[1:(n-10)])
#' influenceCoxTest(list(m.cox, mStrata.cox), tanh.transform = FALSE,
#'                  type = "survival", newdata = newdata, times = Utime[1:(n-10)])
#'
#' #### Under H0 (CSC) ####
#' set.seed(1)
#' ff <- ~ f(X1,2) + f(X2,-0.033)
#' ff <- update(ff, ~ .+ f(X3,0) + f(X4,0) + f(X5,0))
#' ff <- update(ff, ~ .+ f(X6,0) + f(X7,0) + f(X8,0) + f(X9,0))
#' d <- sampleData(n, outcome = "competing.risk", formula = ff)
#' setkey(d, time)
#'
#' Utime <- sort(unique(d$time))
#' 
#' ## fit cox models
#' m.CSC <- CSC(Hist(time, event) ~ X1 + X2, data = d)
#' mStrata.CSC <- CSC(Hist(time, event) ~ strata(X1) + X2 + X3, data = d)
#'
#' ## compare models
#' influenceCoxTest(list(m.CSC, mStrata.CSC), tanh.transform = TRUE,
#'                  cause = 1, newdata = data.frame(X1=1,X2=1,X3=1), times = Utime[1:20])
#' influenceCoxTest(list(m.CSC, mStrata.CSC), tanh.transform = FALSE,
#'                  cause = 1, newdata = data.frame(X1=1,X2=1,X3=1), times = Utime[1:20]) 
#' 
#' @export
`influenceCoxTest` <-
  function(object,...) UseMethod("influenceCoxTest")
# }}}

# {{{ influenceCoxTest.list
#' @rdname influenceCoxTest
#' @export
influenceCoxTest.list <- function(object, newdata, times, type, cause, ...){
  
  validCox <- c("coxph","cph","phreg")
  ls.pred <- lapply(object, function(x){
    
    if(any(validCox %in% class(x))){ # Cox
      pred <- predictCox(x,
                         type = type,
                         newdata = newdata,
                         times = times,
                         iid = TRUE,
                         band = FALSE,
                         se = TRUE,
                         log.transform = FALSE)
    }else if("CauseSpecificCox" %in% class(x)){ # CSC
      pred <- predict(x,
                      cause = cause,
                      newdata = newdata,
                      times = times,
                      iid = TRUE,
                      band = FALSE,
                      se = TRUE,
                      log.transform = FALSE)
    }else{
      stop("can only handle Cox and Cause specific Cox models \n")
    }
    return(pred)
  })
  out <- influenceCoxTest(ls.pred[[1]], ls.pred[[2]], ...)
  return(out)
  
}
# }}}

# {{{ influenceCoxTest.default
#' @rdname influenceCoxTest
#' @export
influenceCoxTest.default <- function(object, object2, tanh.transform = FALSE, conf.level = 0.95, nsim.band = 1e4, ...){
  
  ## type
  if(class(object)=="predictCox"){
    if(object$type!=object2$type){
      stop("Cannot compare different types of prediction \n")
    }
    type <- object$type
    if(length(type)!=1){
      stop("Cannot analyse simulatenously several types of predictions \n",
           "here: \"",paste0(type,collapse ="\" \""),"\"\n")
    }
  }else if(class(object)=="predictCSC"){
    type <- "absRisk"
  }else{
    stop("Can only deal with object of class \'predictCox\' or \'predictCSC\'")
  }
  estimator1 <- object[[type]]
  estimator2 <- object2[[type]]
  if(NROW(estimator1)!=1 || NROW(estimator2)!=1){
    stop("Cannot analyse simulatenously predictions conditional on several set of covariates \n")
  }
  
  # transformation
  if(tanh.transform){
    if(type=="cumhazard"){
      stop("tanh.transform cannot be used when the argument \'type\' is set to \"cumhazard\"\n")
    }
  }
  
  ## influence function
  type.iid <- paste0(type,".iid")
  if(type.iid %in% names(object) == FALSE){
    stop("No iid decomposition for type ",type," in object \n")
  }else{
    estimator1.iid <- object[[type.iid]]
  }
  if(type.iid %in% names(object2) == FALSE){
    stop("No iid decomposition for type ",type," in object2 \n")
  }else{
    estimator2.iid <- object2[[type.iid]]
  }
  
  ## check sample size
  n <- length(estimator1.iid[1,1,])
  if(n!=length(estimator2.iid[1,1,])){
    stop("Cannot compare estimates on different sample size \n")
  }
  
  ## time
  time <- object$time
  n.time <- length(time)
  
  ## compute the difference
  delta <- as.numeric(estimator1-estimator2)
  delta.iid <- matrix((estimator1.iid-estimator2.iid)[1,,], nrow = n, ncol = n.time, byrow = TRUE)
  
  if(tanh.transform){
    delta.iid <- rowMultiply_cpp(delta.iid, scale = 1/(1+delta^2))
    delta <- atanh(delta)
  }
  
  delta.se <- sqrt(colSums(delta.iid^2))
  
  ## test (punctual)
  zval <- qnorm(1-(1-conf.level)/2, 0,1)
  tableRes <- data.frame(cbind(time = time,
                               delta = delta,
                               se = delta.se,
                               inf = delta - zval*delta.se,
                               sup = delta + zval*delta.se,
                               z = delta/delta.se))
  if(tanh.transform){
    tableRes$delta <- tanh(tableRes$delta)
    tableRes$inf <- tanh(tableRes$inf)
    tableRes$sup <- tanh(tableRes$sup)
  }else{
    tableRes$inf <- pmax(tableRes$inf,-1)
    tableRes$sup <- pmin(tableRes$sup,1)
  }
  tableRes$p <- 2*(1-pnorm(abs(tableRes$z)))  
  tableRes$p[delta==0] <- 1 # otherwise 0/0 returns NA
  
  index.NA <- unique(c(which(is.na(delta)),which(is.infinite(delta)),which(is.na(delta.se))))
  
  ## test (band)
  if(length(time)>0){
    
    if(length(index.NA)>0){
      n.timeNNA <- n.time - length(index.NA)
      A.iid <- array(t(delta.iid[,-index.NA,drop=FALSE]), dim = c(n.timeNNA,n,1))
      M.se <- matrix(delta.se[-index.NA], ncol = 1)
    }else{
      A.iid <- array(t(delta.iid), dim = c(n.time,n,1))
      M.se <- matrix(delta.se, ncol = 1)
    }
    
    sample.max <- sampleMaxProcess_cpp(nObject = n,
                                       nNew = 1,
                                       nSim = nsim.band,
                                       iid = A.iid,
                                       se = M.se)
    
    quantile.band <- quantile(sample.max, 0.95)
    tableRes$infBand <- delta - quantile.band*delta.se
    tableRes$supBand <- delta + quantile.band*delta.se
    
    if(tanh.transform){
      tableRes$infBand <- tanh(tableRes$infBand)
      tableRes$supBand <- tanh(tableRes$supBand)
    }else{
      tableRes$infBand <- pmax(tableRes$infBand,-1)
      tableRes$supBand <- pmin(tableRes$supBand,1)
    }
    
    if(length(index.NA)>0){
      pBand <- mean(sample.max>max(abs(tableRes$z[-index.NA]), na.rm=TRUE))
    }else{
      pBand <- mean(sample.max>max(abs(tableRes$z), na.rm=TRUE))
    }
    
  }else{
    quantile.band <- NA
    pBand <- NA
  }
  
  ## export
  tableRes[index.NA,-(1:2)] <- NA # otherwise 0/0 returns NA
  
  out <- list(table = tableRes, 
              quantile.band = quantile.band,
              pBand = pBand)
  class(out) <- "influenceCoxTest"
  
  return(out)
}
# }}}

#----------------------------------------------------------------------
### influence2.R ends here
