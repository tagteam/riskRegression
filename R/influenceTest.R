### influence2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 19 2017 (16:01) 
## Version: 
## last-updated: jul 18 2017 (15:56) 
##           By: Brice Ozenne
##     Update #: 97
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
#' @name influenceTest
#' 
#' @param object either a list of models or an object of class predictCox or predictCSC.
#' @param object2 same as predict1 but for another model.
#' @param type the type of predicted value. 
#' @param nSim.band the number of simulations used to compute the quantiles
#' for the confidence bands.
#' @param newdata A \code{data.frame} or \code{data.table} containing
#'     the values of the predictor variables defining subject specific
#'     predictions.
#' @param times Time points at which to evaluate the predictions.
#' @param cause Identifies the cause of interest among the competing
#'     events.
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
#' influenceCoxTest(list(m.cox, mStrata.cox),
#'                  type = "survival", newdata = newdata, times = seq(0.1,1,by=0.1))
#'
#' #### Under H0 (Cox) ####
#' set.seed(1)
#' d <- sampleData(n, outcome = "survival")
#' setkey(d, eventtime)
#' 
#' ## fit cox models
#' m.cox <- coxph(Surv(eventtime, event) ~ X1, 
#'                data = d, y = TRUE, x = TRUE)
#'
#' mStrata.cox <- coxph(Surv(eventtime, event) ~ strata(X1), 
#'                      data = d, y = TRUE, x = TRUE)
#'
#' ## compare models
#' resH0 <- influenceCoxTest(list(m.cox, mStrata.cox),
#'                  type = "survival", newdata = newdata, times = unique(d$eventtime))
#' resH0[1:10,]
#'
#' #### Under H0 (CSC) ####
#' set.seed(1)
#' d <- sampleData(n, outcome = "competing.risk")
#' setkey(d, time)
#' 
#' ## fit cox models
#' m.CSC <- CSC(Hist(time, event) ~ X1 + X2, data = d)
#' mStrata.CSC <- CSC(Hist(time, event) ~ strata(X1) + X2, data = d)
#'
#' ## compare models
#' resH0 <- influenceCoxTest(list(m.CSC, mStrata.CSC),
#'                  cause = 1, newdata = data.frame(X1=1,X2=1), times = unique(d$time))
#' resH0[1:10,]
#' 
#' @export
`influenceCoxTest` <-
    function(object,...) UseMethod("influenceCoxTest")
# }}}

#' @rdname influenceTest
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
                               logTransform = FALSE)
        }else if("CauseSpecificCox" %in% class(x)){ # CSC
            pred <- predict(x,
                            cause = cause,
                            newdata = newdata,
                            times = times,
                            iid = TRUE,
                            band = FALSE,
                            se = TRUE,
                            logTransform = FALSE)
        }else{
            stop("can only handle Cox and Cause specific Cox models \n")
        }
        return(pred)
    })
    out <- influenceCoxTest(ls.pred[[1]], ls.pred[[2]])
    return(out)
   
}

#' @rdname influenceTest
#' @export
influenceCoxTest.default <- function(object, object2, nSim.band = 1e4, ...){

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
        stop("can only deal with object of class \'predictCox\' or \'predictCSC\'")
    }
    estimator1 <- object[[type]]
    estimator2 <- object2[[type]]
    if(NROW(estimator1)!=1 || NROW(estimator2)!=1){
        stop("Cannot analyse simulatenously predictions conditional on several set of covariates \n")
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
        stop("cannot compare estimates on different sample size \n")
    }

    ## time
    time <- object$time
    n.time <- length(time)
    
    ## log transform
    if(object$logTransform || object2$logTransform){
        stop("No compatible (yet!) with log-transformed influence function \n")
    }

    ## prepare
    delta <- as.numeric(estimator1-estimator2)
    delta.iid <- t((estimator1.iid-estimator2.iid)[1,,])
    delta.se <- sqrt(colSums(delta.iid^2))
    
    ## test (punctual)
    out <- data.frame(cbind(time = time,
                            delta = delta,
                            se.delta = delta.se,
                            z.delta = delta/delta.se))
    out$p <- 2*(1-pnorm(abs(out$z.delta)))  

    ## test (band)
    if(length(time)>0){

        sample.max <- sampleMaxProcess_cpp(nObject = n,
                                           nNew = 1,
                                           nSim = nSim.band,
                                           iid = array(t(delta.iid), dim = c(n.time,n,1)),
                                           se = matrix(delta.se, ncol = 1))

        out$zBand.delta <- quantile(sample.max, 0.95)
        out$pBand.delta <- mean(sample.max>max(abs(out$z.delta), na.rm = TRUE))
                
        
    }
    
  return(out)
}

#----------------------------------------------------------------------
### influence2.R ends here
