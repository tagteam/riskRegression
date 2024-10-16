## * predictCox (documentation)
#' @title Survival probabilities, hazards and cumulative hazards from Cox regression models 
#' @name predictCox
#' 
#' @description Routine to get baseline hazards and subject specific hazards
#' as well as survival probabilities from a \code{survival::coxph} or \code{rms::cph} object.
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata [data.frame or data.table]  Contain the values of the predictor variables
#' defining subject specific predictions.
#' Should have the same structure as the data set used to fit the \code{object}.
#' @param times [numeric vector] Time points at which to return
#' the estimated hazard/cumulative hazard/survival.
#' @param centered [logical] If \code{TRUE} return prediction at the
#'     mean values of the covariates \code{fit$mean}, if \code{FALSE}
#'     return a prediction for all covariates equal to zero.  in the
#'     linear predictor. Will be ignored if argument \code{newdata} is
#'     used. For internal use.
#' @param type [character vector] the type of predicted value. Choices are \itemize{
#'     \item \code{"hazard"} the baseline hazard function when
#'     argument \code{newdata} is not used and the hazard function
#'     when argument \code{newdata} is used.  \item \code{"cumhazard"}
#'     the cumulative baseline hazard function when argument
#'     \code{newdata} is not used and the cumulative hazard function
#'     when argument \code{newdata} is used.  \item \code{"survival"}
#'     the survival baseline hazard function when argument
#'     \code{newdata} is not used and the cumulative hazard function
#'     when argument \code{newdata} is used.  } Several choices can be
#'     combined in a vector of strings that match (no matter the case)
#'     strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param keep.strata [logical] If \code{TRUE} add the (newdata) strata
#'     to the output. Only if there any.
#' @param keep.times [logical] If \code{TRUE} add the evaluation times
#'     to the output.
#' @param keep.newdata [logical] If \code{TRUE} add the value of the
#'     covariates used to make the prediction in the output list.
#' @param keep.infoVar [logical] For internal use.
#' @param se [logical] If \code{TRUE} compute and add the standard errors to the output.
#' @param band [logical] If \code{TRUE} compute and add the quantiles for the confidence bands to the output.
#' @param iid [logical] If \code{TRUE} compute and add the influence function to the output.
#' @param confint [logical] If \code{TRUE} compute and add the confidence intervals/bands to the output.
#' They are computed applying the \code{confint} function to the output.
#' @param diag [logical] when \code{FALSE} the hazard/cumlative hazard/survival for all observations at all times is computed,
#' otherwise it is only computed for the i-th observation at the i-th time.
#' @param average.iid [logical] If \code{TRUE} add the average of the influence function over \code{newdata} to the output.
#' @param product.limit [logical]. If \code{TRUE} the survival is computed using the product limit estimator.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
#' @param store [vector of length 2] Whether prediction should only be computed for unique covariate sets and mapped back to the original dataset (\code{data="minimal"}) and whether the influence function should be stored in a memory efficient way (\code{iid="minimal"}). Otherwise use \code{data="full"} and/or \code{iid="full"}.
#' 
#' @details
#' When the argument \code{newdata} is not specified, the function computes the baseline hazard estimate.
#' See (Ozenne et al., 2017) section "Handling of tied event times".
#'
#' Otherwise the function computes survival probabilities with confidence intervals/bands.
#' See (Ozenne et al., 2017) section "Confidence intervals and confidence bands for survival probabilities".
#' The survival is computed using the exponential approximation (equation 3).
#'
#' A detailed explanation about the meaning of the argument \code{store = c(iid="full")} can be found
#' in (Ozenne et al., 2017) Appendix B "Saving the influence functions".
#' 
#' The function is not compatible with time varying predictor variables.
#' 
#' The centered argument enables us to reproduce the results obtained with the \code{basehaz}
#' function from the survival package but should not be modified by the user.
#'
#' @return The iid decomposition is output in an attribute called \code{"iid"},
#' using an array containing the value of the influence of each subject used to fit the object (dim 1), for each subject in newdata (dim 3), and each time (dim 2).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#'
#' @references
#' Brice Ozenne, Anne Lyngholm Sorensen, Thomas Scheike, Christian Torp-Pedersen and Thomas Alexander Gerds.
#' riskRegression: Predicting the Risk of an Event using Cox Regression Models.
#' The R Journal (2017) 9:2, pages 440-460.
#' 
#' @seealso
#' \code{\link{confint.predictCox}} to compute confidence intervals/bands.
#' \code{\link{autoplot.predictCox}} to display the predictions.

## * predictCox (examples)
#' @examples
##' 
#' library(survival)
#' library(data.table)
#'
#' #######################
#' #### generate data ####
#' #######################
#' 
#' set.seed(10)
#' d <- sampleData(50,outcome="survival") ## training dataset
#' nd <- sampleData(5,outcome="survival") ## validation dataset
#'
#' ## add ties
#' d$time.round <- round(d$time,1)
#' any(duplicated(d$time.round))
#'
#' ## add categorical variables
#' d$U <- sample(letters[1:3],replace=TRUE,size=NROW(d))
#' d$V <- sample(letters[5:8],replace=TRUE,size=NROW(d))
#' nd <- cbind(nd,d[sample(NROW(d),replace=TRUE,size=NROW(nd)),c("U","V")])
#' 
#' ######################
#' #### Kaplan Meier ####
#' ######################
#' 
#' #### no ties
#' fit.KM <- coxph(Surv(time,event)~ 1, data=d, x = TRUE, y = TRUE)
#' predictCox(fit.KM) ## exponential approximation
#' ePL.KM <- predictCox(fit.KM, product.limit = TRUE) ## product-limit
#' as.data.table(ePL.KM)
#' 
#' range(survfit(Surv(time,event)~1, data = d)$surv - ePL.KM$survival) ## same
#'
#' #### with ties (exponential approximation, Efron estimator for the baseline hazard)
#' fit.KM.ties <- coxph(Surv(time.round,event)~ 1, data=d, x = TRUE, y = TRUE)
#' predictCox(fit.KM)
#' 
#' ## retrieving survfit results with ties
#' fit.KM.ties <- coxph(Surv(time.round,event)~ 1, data=d,
#'                      x = TRUE, y = TRUE, ties = "breslow")
#' ePL.KM.ties <- predictCox(fit.KM, product.limit = TRUE)
#' range(survfit(Surv(time,event)~1, data = d)$surv - ePL.KM.ties$survival) ## same
#'
#' #################################
#' #### Stratified Kaplan Meier ####
#' #################################
#' 
#' fit.SKM <- coxph(Surv(time,event)~strata(X2), data=d, x = TRUE, y = TRUE)
#' ePL.SKM <- predictCox(fit.SKM, product.limit = TRUE)
#' ePL.SKM
#' 
#' range(survfit(Surv(time,event)~X2, data = d)$surv - ePL.SKM$survival) ## same
#'
#' ###################
#' #### Cox model ####
#' ###################
#'
#' fit.Cox <- coxph(Surv(time,event)~X1 + X2 + X6, data=d, x = TRUE, y = TRUE)
#' 
#' #### compute the baseline cumulative hazard
#' Cox.haz0 <- predictCox(fit.Cox)
#' Cox.haz0
#' 
#' range(survival::basehaz(fit.Cox)$hazard-Cox.haz0$cumhazard) ## same
#'
#' #### compute individual specific cumulative hazard and survival probabilities
#' ## exponential approximation
#' fit.predCox <- predictCox(fit.Cox, newdata=nd, times=c(3,8),
#'                           se = TRUE, band = TRUE)
#' fit.predCox
#' 
#' ## product-limit
#' fitPL.predCox <- predictCox(fit.Cox, newdata=nd, times=c(3,8),
#'                             product.limit = TRUE, se = TRUE)
#' fitPL.predCox
#'
#' ## Note: product limit vs. exponential does not affect uncertainty quantification
#' range(fitPL.predCox$cumhazard.se - fit.predCox$cumhazard.se) ## same
#' range(fitPL.predCox$survival.se - fit.predCox$survival.se) ## different
#' ## except through a multiplicative factor (multiply by survival in the delta method)
#' fitPL.predCox$survival.se - fitPL.predCox$cumhazard.se * fitPL.predCox$survival
#'
#'
#' #### left truncation
#' test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
#'               stop=c(2,3,6,7,8,9,9,9,14,17), 
#'               event=c(1,1,1,1,1,1,1,0,0,0), 
#'               x=c(1,0,0,1,0,1,1,1,0,0)) 
#' m.cph <- coxph(Surv(start, stop, event) ~ 1, test2, x = TRUE)
#' as.data.table(predictCox(m.cph))
#'
#' basehaz(m.cph)
#'
#' ##############################
#' #### Stratified Cox model ####
#' ##############################
#' 
#' #### one strata variable
#' fit.SCox <- coxph(Surv(time,event)~X1+strata(X2)+X6, data=d, x = TRUE, y = TRUE)
#' predictCox(fit.SCox, newdata=nd, times = c(3,12))
#' predictCox(fit.SCox, newdata=nd, times = c(3,12), product.limit = TRUE)
#' ## NA: timepoint is beyond the last observation in the strata and censoring
#' tapply(d$time,d$X2,max)
#' 
#' #### two strata variables
#' fit.S2Cox <- coxph(Surv(time,event)~X1+strata(U)+strata(V)+X2,
#'                    data=d, x = TRUE, y = TRUE)
#'
#' base.S2Cox <- predictCox(fit.S2Cox)
#' range(survival::basehaz(fit.S2Cox)$hazard - base.S2Cox$cumhazard) ## same
#' predictCox(fit.S2Cox, newdata=nd, times = 3)

## * predictCox (code)
#' @rdname predictCox
#' @export
predictCox <- function(object,
                       times,
                       newdata = NULL,
                       centered = TRUE,
                       type=c("cumhazard","survival"),
                       keep.strata = TRUE,
                       keep.times = TRUE,
                       keep.newdata = FALSE,
                       keep.infoVar = FALSE,
                       se = FALSE,
                       band = FALSE,
                       iid = FALSE,
                       confint = (se+band)>0,
                       diag = FALSE,
                       average.iid = FALSE,
                       product.limit = FALSE,
                       store = NULL){

    call <- match.call()
    newdata <- data.table::copy(newdata)
    ## centering
    if(!is.null(newdata)){
        if(inherits(centered,"data.frame")){
            df.reference <- centered
            centered2 <- TRUE ## for the linear predictor of the hazard
        }else{
            df.reference <- NULL
            centered2 <- centered ## for the linear predictor of the hazard
        }
        centered <- TRUE ## for the linear predictor of the baseline hazard
    }else{
        centered2 <- FALSE
    }
    
    ## ** Extract elements from object
    if (missing(times)) {
        nTimes <- 0
        times <- numeric(0)
    }else{
        nTimes <- length(times)
    }
    needOrder <- (nTimes[1]>0 && is.unsorted(times))
    if (all(!is.na(times)) && needOrder) {
        order.times <- order(times)
        oorder.times <- order(order.times)
        times.sorted <- sort(times)
    }else{
        if (nTimes==0){
            times.sorted <- numeric(0)
            order.times <- integer(0)
            oorder.times <- integer(0)
        }else{
            times.sorted <- times
            order.times <- 1:nTimes
            oorder.times <- 1:nTimes
        }
    }
    object.n <- coxN(object)  
    object.modelFrame <- coxModelFrame(object)
    infoVar <- coxVariableName(object, model.frame = object.modelFrame)
    object.baseEstimator <- coxBaseEstimator(object)

    ## ease access
    is.strata <- infoVar$is.strata
    object.levelStrata <- levels(object.modelFrame$strata) ## levels of the strata variable
    nStrata <- length(object.levelStrata) ## number of strata
    nVar.lp <- length(infoVar$lpvars) ## number of variables in the linear predictor

    ## ** normalize model frame
    ## convert strata to numeric
    object.modelFrame[,c("strata.num") := as.numeric(.SD$strata) - 1]

    ## linear predictor
    ## if we predict the hazard for newdata then there is no need to center the covariates
    object.modelFrame[,c("eXb") := exp(coxLP(object, data = NULL, center = if(is.null(newdata)){centered}else{FALSE}))]
    ## add linear predictor and remove useless columns
    rm.name <- setdiff(names(object.modelFrame),c("start","stop","status","eXb","strata","strata.num"))
    if(length(rm.name)>0){
        object.modelFrame[,c(rm.name) := NULL]
    }
  
    ## sort the data
    object.modelFrame[, c("statusM1") := 1-.SD$status] ## sort by statusM1 such that deaths appear first and then censored events
    object.modelFrame[, c("XXXindexXXX") := 1:.N] ## keep track of the initial positions (useful when calling calcSeCox)

    if(inherits(object,"prodlim") && object$reverse){
        reverse <- TRUE
        data.table::setkeyv(object.modelFrame, c("strata.num","stop","start","status"))
    }else{
        reverse <- FALSE
        data.table::setkeyv(object.modelFrame, c("strata.num","stop","start","statusM1"))
    }
    
    ## last event time in each strata
    if(!is.null(attr(times,"etimes.max"))){ ## specified by the user
        etimes.max <- attr(times,"etimes.max")
        attr(times,"etimes.max") <- NULL
        attr(times.sorted,"etimes.max") <- etimes.max
    }else if(is.strata){ ## based on the data
        if(nVar.lp==0){
            iDTtempo <- object.modelFrame[, .SD[which.max(.SD$stop)], by = "strata.num"]
            etimes.max <- iDTtempo[,if(.SD$status==1){1e12}else{.SD$stop}, by = "strata.num"][[2]]
        }else{
            etimes.max <- object.modelFrame[, max(.SD$stop), by = "strata.num"][[2]]
        }
    }else{
        if(nVar.lp==0 && (utils::tail(object.modelFrame$status,1)==1)){ ## no covariates and ends by a death
            etimes.max <- 1e12
        }else{
            etimes.max <- max(object.modelFrame[["stop"]])
        }
    }

    ## ** checks
    ## *** object
    if((reverse == FALSE) && ("risk" %in% object$type)){
        stop("predictCox can only deal with prodlim objects applied to survival data or targeting the censoring distribution. \n")
    }
    ## predictCox is not compatible with all coxph/cph object (i.e. only handle only simple cox models)
    if(!is.null(object$weights) && !all(object$weights==1)){
        stop("predictCox does not know how to handle Cox models fitted with weights \n")
    }
    if(!is.null(object$naive.var)){
        stop("predictCox does not know how to handle frailty.") 
    }
    if(object.baseEstimator == "exact"){
        stop("Prediction with exact handling of ties is not implemented.\n")
    }
    if(!is.null(object$call$tt)){
        stop("predictCox does not know how to handle time varying effects.\n") 
    }
    ## convergence issue
    if(!is.null(coef(object)) && any(is.na(coef(object)))){
        print(coef(object))
        stop("One or several parameters of the regression model have no value, i.e., a value 'NA'.\n")
    }

    ## *** times
    if(nTimes[1]>0 && any(is.na(times))){
        stop("Missing (NA) values in argument \'times\' are not allowed.\n")
    }

    ## *** type
    type <- tolower(type)
    if(any(type %in% c("lp","hazard","cumhazard","survival") == FALSE)){
        stop("type can only be \"lp\", \"hazard\", \"cumhazard\" or/and \"survival\" \n") 
    }
    if(is.null(newdata) && "lp" %in% type){
        stop("Cannot evaluate the linear predictor when argument \'newdata\' is missing. \n")
    }
    if(length(times)>1 && "lp" %in% type){
        stop("Cannot evaluate the linear predictor when there are multiple timepoints. \n")
    }

    ## *** newdata 
    if (missing(newdata) && (se || iid || average.iid)){
        stop("Argument 'newdata' is missing. Cannot compute standard errors in this case.")
    }
    if(!is.null(newdata)){
        if(nTimes[1]==0 && !identical(as.character(type),"lp")){
            stop("Time points at which to evaluate the predictions are missing \n")
        }
        if(!is.vector(times)){
            stop("Argument \'times\' must be a vector \n")
        }
        name.regressor <- c(infoVar$lpvars.original, infoVar$stratavars.original)
        if(length(name.regressor) > 0 && any(name.regressor %in% names(newdata) == FALSE)){
            stop("Missing variables in argument \'newdata\': \"",
                 paste0(setdiff(name.regressor,names(newdata)), collapse = "\" \""),
                 "\"\n")
        }
        if(se[1] && ("hazard" %in% type)){
            stop("Standard error cannot be computed for the hazard \n")
        }
        if(band[1] && ("hazard" %in% type)){
            stop("Confidence bands cannot be computed for the hazard \n")
        }
    }

    ## *** diag
    if(!is.logical(diag)){
        stop("Argument \'diag\' must be of type logical \n")
    }
    if(diag){
        if(NROW(newdata)!=length(times)){
            stop("When argument \'diag\' is TRUE, the number of rows in \'newdata\' must equal the length of \'times\' \n")
        }
    }

    ## *** store
    store.data <- NULL
    store.iid <- "full"
    if(!is.null(store)){
        if(length(store) > 2){
            stop("Argument \'store\' should contain at most two elements. \n",
                 "For instance store = c(data = \"full\", iid = \"full\") or store = c(data = \"minimal\", iid = \"minimal\").\n")
        }
        if(is.null(names(store)) || any(names(store) %in% c("data","iid") == FALSE)){
            stop("Incorrect names for argument \'store\': should \"data\" and \"iid\". \n",
                 "For instance store = c(data = \"full\", iid = \"full\") or store = c(data = \"minimal\", iid = \"minimal\").\n")
        }
        if("data" %in% names(store)){
            if(store[["data"]] %in% c("minimal","full") == FALSE){
                stop("Element in argument \'store\' should take value \'minimal\' or \'full\'.\n",
                     "For instance store = c(data = \"full\") or store = c(data = \"minimal\").\n")
            }
            store.data <- store[["data"]]
        }
        if("iid" %in% names(store)){
            if(store[["iid"]] %in% c("minimal","full") == FALSE){
                stop("Element in argument \'store\' should take value \'minimal\' or \'full\'.\n",
                     "For instance store = c(iid = \"full\") or store = c(iid = \"minimal\").\n")
            }
            store.iid <- store[["iid"]]
        }
    }

    ## *** average.iid
    if(average.iid==TRUE && !is.null(attr(average.iid,"factor"))){
        if(is.null(store.iid) && !is.null(object$iid$store.iid)){
            store.iid <- object$iid$store.iid
        }
        test.list <- !is.list(attr(average.iid,"factor"))
        if(test.list){
            stop("Attribute \"factor\" of argument \'average.iid\' must be a list \n")
        }
        test.matrix <- any(unlist(lapply(attr(average.iid,"factor"), is.matrix))==FALSE)
        if(test.matrix){
            stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices \n")
        }
        for(iFactor in 1:length(attr(average.iid,"factor"))){ ## iFactor <- 1
            ## when only one column and diag = FALSE, use the same weights at all times
            if((diag == FALSE) && (NCOL(attr(average.iid,"factor")[[iFactor]])==1) && (nTimes > 1)){
                attr(average.iid,"factor")[[iFactor]] <- matrix(attr(average.iid,"factor")[[iFactor]][,1],
                                                                nrow = NROW(attr(average.iid,"factor")[[iFactor]]),
                                                                ncol = nTimes, byrow = FALSE)
            }
            ## check dimensions
            if(any(dim(attr(average.iid,"factor")[[iFactor]])!=c(NROW(newdata), diag + (1-diag)*nTimes))){
                stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices of size ",new.n,",",diag + (1-diag)*nTimes," \n")
            }
        }
    }

    ## ** compress newdata into unique patient profile
    outCompress <- compressData(object, newdata = newdata, times = times, diag = diag, average.iid = average.iid,
                                oorder.times = oorder.times, times.sorted = times.sorted,
                                nStrata = nStrata, infoVar = infoVar, level = store.data)
    if(!is.null(outCompress)){
        newdata <- outCompress$newdata
        times.sorted <- outCompress$times.sorted
        diag <- outCompress$diag
        order.times <- outCompress$order.times
        oorder.times <- outCompress$oorder.times
        nTimes <- outCompress$nTimes
        average.iid <- outCompress$average.iid
    }

    ## ** baseline hazard
    ## *** evaluate at requested times
    ## if(inherits(object,"prodlim")){
    ##     Lambda0 <- baseHaz_prodlim(object,
    ##                                times = times.sorted,
    ##                                etimes.max = etimes.max)
    ## }else{
    Lambda0 <- baseHaz_cpp(starttimes = object.modelFrame$start,
                           stoptimes = object.modelFrame$stop,
                           status = object.modelFrame$status,
                           eXb = object.modelFrame$eXb,
                           strata = object.modelFrame$strata.num,
                           nPatients = object.n,
                           nStrata = nStrata,
                           emaxtimes = etimes.max,
                           predtimes = times.sorted,
                           cause = 1,
                           Efron = (object.baseEstimator == "efron"),
                           reverse = reverse)
    ## }

    ## *** evaluate at all jump times    
    if(product.limit){        
        if(length(times.sorted)==0){
            AllLambda0 <- Lambda0
            ## }else if(inherits(object,"prodlim")){
            ##     AllLambda0 <- baseHaz_prodlim(object,
            ##                                   times = numeric(0),
            ##                                   etimes.max = pmin(etimes.max,max(c(times.sorted,Inf))) ## times.sorted may be numeric(0)
            ##                                   )
        }else{
            AllLambda0 <- baseHaz_cpp(starttimes = object.modelFrame$start,
                                      stoptimes = object.modelFrame$stop,
                                      status = object.modelFrame$status,
                                      eXb = object.modelFrame$eXb,
                                      strata = object.modelFrame$strata.num,
                                      nPatients = object.n,
                                      nStrata = nStrata,
                                      emaxtimes = pmin(etimes.max,max(c(times.sorted,Inf))), ## times.sorted may be numeric(0),
                                      predtimes = numeric(0),
                                      cause = 1,
                                      Efron = (object.baseEstimator == "efron"),
                                      reverse = reverse)
        }
    }
    
    ## *** post-process
    ## restaure strata levels (numeric to character)
    Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
    if(product.limit && "survival" %in% type){
        AllLambda0$strata <- factor(AllLambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
    }
    ## update AllLambda with time 0, timepoints at which the hazard can no more be computed (NA), and evaluate baseline survival
    if(product.limit && "survival" %in% type){
        ls.AllLambda0 <- lapply(levels(AllLambda0$strata), function(iS){ ## iS <- levels(AllLambda0$strata)[1]
            iIndex <- which(AllLambda0$strata==iS)
            iOut <- list(times = AllLambda0$times[iIndex],
                         hazard = AllLambda0$hazard[iIndex],
                         cumhazard = AllLambda0$cumhazard[iIndex],
                         strata = AllLambda0$strata[iIndex])
            if(iOut$times[1]>0 && (!is.null(newdata) || length(times.sorted)>0)){ ## add time 0 when it is not there
                iOut$times <- c(0,iOut$times)
                iOut$hazard <- c(0,iOut$hazard)
                iOut$cumhazard <- c(0,iOut$cumhazard)
                iOut$strata <- c(iOut$strata[1],iOut$strata)
            }
            if(any(is.na(Lambda0$hazard[Lambda0$strata==iS])) && (!is.null(newdata) || length(times.sorted)>0)){ ## add time where the first NA occured
                iOut$times <- c(iOut$times, etimes.max[levels(Lambda0$strata)==iS]+1e-12)
                iOut$hazard <- c(iOut$hazard,NA)
                iOut$cumhazard <- c(iOut$cumhazard,NA)
                iOut$strata <- c(iOut$strata,iOut$strata[1])
            }
            if(is.null(newdata)){ ## compute baseline survival, i.e., survival for covariate values of 0
                iOut$survival <- cumprod(1-iOut$hazard)
            }
            return(iOut)        
        })
        AllLambda0 <- list(times = unlist(lapply(ls.AllLambda0,"[[","times"), use.names = FALSE),
                           hazard = unlist(lapply(ls.AllLambda0,"[[","hazard"), use.names = FALSE),
                           cumhazard = unlist(lapply(ls.AllLambda0,"[[","cumhazard"), use.names = FALSE),
                           survival = unlist(lapply(ls.AllLambda0,"[[","survival"), use.names = FALSE),
                           strata = unlist(lapply(ls.AllLambda0,"[[","strata"), use.names = FALSE))
    }

    ## *** special case: return baseline hazard/cumulative hazard/survival
    if (is.null(newdata)){
        if (!("hazard" %in% type)){
            Lambda0$hazard <- NULL
        }
        if ("survival" %in% type){  ## must be before cumhazard
            if(product.limit){
                if(length(times.sorted)==0){
                    Lambda0$survival <- AllLambda0$survival
                }else{
                    Lambda0$survival <- unlist(lapply(unique(AllLambda0$strata), function(iS){ ## iS <- 1
                        AllLambda0$survival[AllLambda0$strata==iS][prodlim::sindex(jump.times = AllLambda0$times[AllLambda0$strata==iS], eval.times = times.sorted[oorder.times])]
                    }))
                }
            }else{
                Lambda0$survival <- exp(-Lambda0$cumhazard)
            }
        }else{
            Lambda0$survival <- NULL
        }
        if (!("cumhazard" %in% type)){
            Lambda0$cumhazard <- NULL
        } 
        if (keep.times==FALSE){
            Lambda0$times <- NULL
        }
        if (keep.strata[[1]]==FALSE ||(is.null(call$keep.strata) && !is.strata)){
            Lambda0$strata <- NULL
        }
        if( keep.newdata[1]==TRUE){
            Lambda0$newdata <- object.modelFrame
        }
        add.list <- list(lastEventTime = etimes.max,
                         se = FALSE,
                         band = FALSE,                         
                         type = type,
                         nTimes = nTimes,
                         baseline = TRUE,
                         var.lp = infoVar$lpvars.original,
                         var.strata = infoVar$stratavars.original)
        if(keep.infoVar){
            add.list$infoVar <- infoVar
        }
        
        Lambda0[names(add.list)] <- add.list
        class(Lambda0) <- "predictCox"
        return(Lambda0)

    } 

    ## ** compute subject specific hazard / cumlative hazard / survival
    out <- list()

    ## *** reformat newdata (compute linear predictor and strata)
    new.n <- NROW(newdata)
    newdata <- data.table::as.data.table(newdata)

    Xb <- coxLP(object, data = newdata, center = FALSE)
    if ("lp" %in% type){
        out$lp <- cbind(Xb)
        lp.iid <- centered2 ## when ask for centered then we need the iid for exporting the correct se
    }else{
        lp.iid <- FALSE
    }
    new.eXb <- exp(Xb)
    new.strata <- coxStrata(object, data = newdata, 
                            sterms = infoVar$strata.sterms, 
                            strata.vars = infoVar$stratavars, 
                            strata.levels = infoVar$strata.levels)
    new.levelStrata <- levels(droplevels(new.strata))

    ## *** estimates
    if (is.strata==FALSE && any(c("hazard","cumhazard","survival") %in% type)){ ## no strata
        if(diag && (!product.limit || "hazard" %in% type || "cumhazard" %in% type)){
            iTimes <- prodlim::sindex(jump.times = Lambda0$times, eval.times = times.sorted[oorder.times])
        }
        if("survival" %in% type && product.limit){
            iAllTimes <- prodlim::sindex(jump.times = AllLambda0$times, eval.times = times.sorted[oorder.times])
        }                
            
        if ("hazard" %in% type){
            if(diag){
                out$hazard <- cbind(new.eXb * Lambda0$hazard[iTimes])
            }else{
                out$hazard <- (new.eXb %o% Lambda0$hazard[oorder.times])
            }
        }

        if ("cumhazard" %in% type || ("survival" %in% type && product.limit == FALSE)){
            if(diag){
                cumhazard <- cbind(new.eXb * Lambda0$cumhazard[iTimes])
            }else{
                cumhazard <- new.eXb %o% Lambda0$cumhazard[oorder.times]
            }
            if ("cumhazard" %in% type){
                out$cumhazard <- cumhazard
            }
            if ("survival" %in% type && product.limit == FALSE){
                out$survival <- exp(-cumhazard)
            }
        }

        if ("survival" %in% type && product.limit){
            if(diag){
                out$survival <- cbind(mapply(x = new.eXb, t = iAllTimes, FUN = function(x,t){
                    prod(1 - x * AllLambda0$hazard[1:t])
                }))
            }else{
                out$survival <- do.call(rbind,apply(1 - new.eXb %o% AllLambda0$hazard, MARGIN = 1, cumprod, simplify = FALSE))[,iAllTimes,drop=0L]
                    
            }
        }

    }else if(any(c("hazard","cumhazard","survival") %in% type)){ ## strata

        ## initialization
        if ("hazard" %in% type){
            out$hazard <- matrix(0, nrow = new.n, ncol = nTimes*(1-diag)+diag)
        }
        if ("cumhazard" %in% type){
            out$cumhazard <- matrix(NA, nrow = new.n, ncol = nTimes*(1-diag)+diag)                
        }
        if ("survival" %in% type){
            out$survival <- matrix(NA, nrow = new.n, ncol = nTimes*(1-diag)+diag)
        }

        ## loop across strata
        for(S in new.levelStrata){ ## S <- 1
            id.S <- which(Lambda0$strata==S)
            newid.S <- which(new.strata==S)
            if(product.limit && "survival" %in% type){
                jump.times.S <- AllLambda0$times[AllLambda0$strata==S]
                hazard0.S <- AllLambda0$hazard[AllLambda0$strata==S]
            }

            if(diag && (!product.limit || "hazard" %in% type || "cumhazard" %in% type)){
                iSTimes <- prodlim::sindex(jump.times = Lambda0$times[id.S], eval.times = times.sorted[oorder.times[newid.S]])                  
            }
            if("survival" %in% type && product.limit){
                if(diag){
                    iSAllTimes <- prodlim::sindex(jump.times = jump.times.S, eval.times = times.sorted[oorder.times[newid.S]])
                }else{
                    iSAllTimes <- prodlim::sindex(jump.times = jump.times.S, eval.times = times.sorted[oorder.times])
                }
            }
        
            if ("hazard" %in% type){
                if(diag){
                    out$hazard[newid.S] <- new.eXb[newid.S] * Lambda0$hazard[id.S][iSTimes]
                }else{
                    out$hazard[newid.S,] <- new.eXb[newid.S] %o% Lambda0$hazard[id.S][oorder.times]
                }
            }
            if ("cumhazard" %in% type || ("survival" %in% type && product.limit == FALSE)){
                if(diag){
                    cumhazard.S <- cbind(new.eXb[newid.S] * Lambda0$cumhazard[id.S][iSTimes])
                }else{
                    cumhazard.S <- new.eXb[newid.S] %o% Lambda0$cumhazard[id.S][oorder.times]
                }
                if ("cumhazard" %in% type){
                    out$cumhazard[newid.S,] <- cumhazard.S
                }
                if ("survival" %in% type && product.limit == FALSE){
                    out$survival[newid.S,] <- exp(-cumhazard.S)
                }
            }
            if("survival" %in% type && product.limit){
                if(diag){
                    out$survival[newid.S,] <- mapply(x = new.eXb[newid.S], t = iSAllTimes, FUN = function(x,t){
                        prod(1 - x * hazard0.S[1:t])
                    })
                }else{
                    out$survival[newid.S,] <- do.call(rbind,apply(1 - new.eXb[newid.S] %o% hazard0.S, MARGIN = 1, cumprod, simplify = FALSE))[,iSAllTimes,drop=0L]
                }
            }
        }
    }
    
    ## *** uncertainty
    if(se[[1]] || band[[1]] || iid[[1]] || average.iid[[1]]){
        if(nVar.lp > 0){
            ## get the (new) design matrix
            new.LPdata <- model.matrix(object, data = newdata)
            if(NROW(new.LPdata)!=NROW(newdata)){
                stop("NROW of the design matrix and newdata differ. \n",
                     "Maybe because newdata contains NA values \n")
            }
            if(any(sort(colnames(new.LPdata))!=sort(names(coef(object))))){
                stop("Names of the design matrix and model parameters differ. \n",
                     "Possible error in model.matrix due to special operator in the formula. \n")
            }
        }else{
            new.LPdata <- matrix(0, ncol = 1, nrow = new.n)
        }

        ## restaure original ordering
        data.table::setkeyv(object.modelFrame,"XXXindexXXX")
        if(diag){
            Lambda0$oorder.times <- oorder.times
        }else{
            Lambda0$oorder.times <- 1:nTimes
        }

        ## Computation of the influence function and/or the standard error
        export <- c("iid"[(iid+band+lp.iid)>0],"se"[(se+band)>0],"average.iid"[average.iid==TRUE])
        if(!is.null(attr(average.iid,"factor"))){
            if(diag){
                attr(export,"factor") <- attr(average.iid,"factor")
            }else{
                ## re-order columns according to times
                attr(export,"factor") <- lapply(attr(average.iid,"factor"), function(iF){ ## iF <- attr(average.iid,"factor")[[1]]
                    iF[,order.times,drop=FALSE]
                })
            }
        }
        if(diag){
            times2 <- times.sorted[oorder.times]
        }else{
            times2 <- times.sorted
        }
        attr(times2,"etimes.max") <- attr(times.sorted,"etimes.max")
        outSE <- calcSeCox(object,
                           times = times2,
                           nTimes = nTimes,
                           type = type,
                           diag = diag,
                           Lambda0 = Lambda0,
                           object.n = object.n,
                           object.time = object.modelFrame$stop,
                           object.eXb = object.modelFrame$eXb,
                           object.strata =  object.modelFrame$strata, 
                           nStrata = nStrata,
                           new.n = new.n,
                           new.eXb = new.eXb,
                           new.LPdata = new.LPdata,
                           new.strata = new.strata,
                           new.survival = if(diag){out$survival}else{out$survival[,order.times,drop=FALSE]},
                           nVar.lp = nVar.lp, 
                           export = export,
                           store.iid = store.iid)

        ## restaure orginal time ordering
        if((iid+band)>0){
            if ("lp" %in% type){
                out$lp.iid <- outSE$lp.iid
            }
            if ("hazard" %in% type){
                if (needOrder[1] && (diag[1] == FALSE))
                    out$hazard.iid <- outSE$hazard.iid[,oorder.times,,drop=0L]
                else
                    out$hazard.iid <- outSE$hazard.iid
            }
            if ("cumhazard" %in% type){
                if (needOrder[1] && (diag[1] == FALSE))
                    out$cumhazard.iid <- outSE$cumhazard.iid[,oorder.times,,drop=0L]
                else
                    out$cumhazard.iid <- outSE$cumhazard.iid
            }
            if ("survival" %in% type){
                if (needOrder[1] && (diag[1] == FALSE))
                    out$survival.iid <- outSE$survival.iid[,oorder.times,,drop=0L]
                else
                    out$survival.iid <- outSE$survival.iid
            }
        }
        if(average.iid == TRUE){
            if("lp" %in% type){
                out$lp.average.iid <- outSE$lp.average.iid
            }
            if ("hazard" %in% type){
                if (needOrder && (diag[1] == FALSE)){
                    if(is.list(outSE$hazard.average.iid)){
                        out$hazard.average.iid <- lapply(outSE$hazard.average.iid, function(iIID){iIID[,oorder.times,drop=0L]})
                    }else{
                        out$hazard.average.iid <- outSE$hazard.average.iid[,oorder.times,drop=0L]
                    }
                }else{
                    out$hazard.average.iid <- outSE$hazard.average.iid
                }
            }
            if ("cumhazard" %in% type){
                if (needOrder && (diag[1] == FALSE)){
                    if(is.list(outSE$cumhazard.average.iid)){
                        out$cumhazard.average.iid <- lapply(outSE$cumhazard.average.iid, function(iIID){iIID[,oorder.times,drop=0L]})
                    }else{
                        out$cumhazard.average.iid <- outSE$cumhazard.average.iid[,oorder.times,drop=0L]
                    }
                }else{
                    out$cumhazard.average.iid <- outSE$cumhazard.average.iid
                }
            }
            if ("survival" %in% type){
                if (needOrder && (diag[1] == FALSE)){
                    if(is.list(outSE$survival.average.iid)){
                        out$survival.average.iid <- lapply(outSE$survival.average.iid, function(iIID){iIID[,oorder.times,drop=0L]})
                    }else{
                        out$survival.average.iid <- outSE$survival.average.iid[,oorder.times,drop=0L]
                    }
                }else {
                    out$survival.average.iid <- outSE$survival.average.iid
                }
            }

            if(is.list(attr(average.iid,"factor"))){
                if("hazard" %in% type){
                    names(out$hazard.average.iid) <- names(attr(average.iid,"factor"))
                }
                if("cumhazard" %in% type){
                    names(out$cumhazard.average.iid) <- names(attr(average.iid,"factor"))
                }
                if("survival" %in% type){
                    names(out$survival.average.iid) <- names(attr(average.iid,"factor"))
                }
            }

        }
        if((se+band)>0){
            if("lp" %in% type){
                out$lp.se <- outSE$lp.se
            }
            if ("cumhazard" %in% type){
                if (needOrder && (diag[1] == FALSE)){
                    out$cumhazard.se <- outSE$cumhazard.se[,oorder.times,drop=0L]
                }else{
                    out$cumhazard.se <- outSE$cumhazard.se
                }          
            }
            if ("survival" %in% type){
                if (needOrder && (diag[1] == FALSE)){
                    out$survival.se <- outSE$survival.se[,oorder.times,drop=0L]
                } else{
                    out$survival.se <- outSE$survival.se
                }          
            }
        }      
    }

    ## ** substract reference
    if("lp" %in% type && centered2){
        if(is.null(df.reference)){
            data <- try(eval(object$call$data), silent = TRUE)
            if(inherits(x=data,what="try-error")){
                stop("Could not evaluate the dataset used to fit the model to define a reference level. \n",
                     "Set argument \'centered\' to FALSE or to a data.frame definining the reference level. \n")
            }
            var.original <- infoVar$lpvars.original
                
                ls.ref <- lapply(var.original, function(iVar){
                    if(is.numeric(data[[iVar]])){
                        return(unname(mean(data[[iVar]])))
                    }else if(is.factor(data[[iVar]])){
                        return(unname(factor(levels(data[[iVar]])[1],levels(data[[iVar]]))))
                    }else if(is.character(data[[iVar]])){
                        return(unname(sort(unique(data[[iVar]])[1])))
                    }
                })
                df.reference <- as.data.frame(setNames(ls.ref,var.original))
            }

            ls.args <- as.list(call)[-1]
            ls.args$newdata <- df.reference
            ls.args$centered <- FALSE
            ls.args$type <- "lp"
            ls.args$se <- (se+band>0)
            ls.args$iid <- (iid+se+band>0)
            ls.args$band <- FALSE
            ls.args$confint <- FALSE
            outRef <- do.call(predictCox, args = ls.args)

        out$lp <- out$lp - as.double(outRef$lp)
        if(band[1] || se[1]){
            out$lp.se <- cbind(sqrt(colSums(colCenter_cpp(outSE$lp.iid, outRef$lp.iid)^2))) ## use outSe$lp.iid instead of out$lp.iid for when not exporting the iid
        }
        if(band[1] || iid[1]){
            out$lp.iid <- colCenter_cpp(out$lp.iid, outRef$lp.iid)
        }
            
    }

    ## ** add information to the predictions
    add.list <- list(lastEventTime = etimes.max,
                     se = se,
                     band = band,
                     type = type,
                     diag = diag,
                     nTimes = nTimes,
                     baseline = FALSE,
                     var.lp = infoVar$lpvars.original,
                     var.strata = infoVar$stratavars.original)
    if (keep.times==TRUE){
        add.list$times <- times
    }
    if( keep.infoVar){
        add.list$infoVar <- infoVar
    }
    out[names(add.list)] <- add.list
    class(out) <- "predictCox"

    ## ** confidence intervals/bands
    if(confint){
        out <- stats::confint(out)
    }
    if(band[1] && se[1]==FALSE){
        out[paste0(type,".se")] <- NULL
    }
    if(band[1] && iid[1]==FALSE){
        out[paste0(type,".iid")] <- NULL
    }

    ## ** retrieve original patient profile from unique patient profiles
    if(!is.null(outCompress)){
        out <- decompressData(out, newdata = newdata, type = type, diag = outCompress$diag.save, times = times, se = se, confint = confint, band = band, iid = iid, average.iid = average.iid,
                              newdata.index = outCompress$newdata.index, times.sorted = times.sorted, needOrder = needOrder)
        newdata <- outCompress$newdata.save
    }
    all.covars <- c(infoVar$stratavars.original,infoVar$lpvars.original)
    if(keep.newdata[1]==TRUE && length(all.covars)>0){
        if(data.table::is.data.table(newdata)){
            out$newdata <- newdata[, all.covars, with = FALSE]
        }else{
            out$newdata <- newdata[, all.covars, drop = FALSE]
        }
    }
    if(keep.strata[1]==TRUE && length(infoVar$stratavars.original)>0){
        if(!is.null(outCompress)){
            out$strata <- new.strata[attr(outCompress$newdata.index,"vectorwise")]
        }else{
            out$strata <- new.strata
        }
    }
    
    ## ** export
    return(out)
}



## * baseHaz_prodlim
baseHaz_prodlim <- function(object, times, etimes.max){

    ## ** check input
    if(!inherits(object,"prodlim")){
        stop("baseHaz_prodlim can only handle prodlim objects. \n")
    }

    ## ** gather baseline hazard, cumulative hazard and corresponding strata and jump times
    n.strata <- length(object$first.strata)

    Lambda0 <- list(times = object$time,
                    hazard = object$hazard,
                    cumhazard = NA,
                    strata = unlist(mapply(s = 0:(n.strata-1), size = object$size.strata, function(s,size){rep(s,size)}, SIMPLIFY = FALSE)))
    Lambda0$cumhazard <- unname(unlist(tapply(Lambda0$hazard, Lambda0$strata, cumsum, simplify = FALSE)))

    ## ** subset to match user specific times
    if(length(times)>0){
        ls.Lambda0 <- lapply(0:(n.strata-1), function(iS){ ## iS <- 0
            iIndex <- which(Lambda0$strata==iS)
            iSindex <- prodlim::sindex(jump.times = Lambda0$times[iIndex], eval.times = times)
            iOut <- list(times = times,
                         hazard = c(0,Lambda0$hazard[iIndex])[iSindex+1],
                         cumhazard = c(0,Lambda0$cumhazard[iIndex])[iSindex+1],
                         strata =  rep(iS, length(times)))
            return(iOut)
        })

        Lambda0 <- list(times = unlist(lapply(ls.Lambda0,"[[","times")),
                        hazard = unlist(lapply(ls.Lambda0,"[[","hazard")),
                        cumhazard = unlist(lapply(ls.Lambda0,"[[","cumhazard")),
                        strata = unlist(lapply(ls.Lambda0,"[[","strata")))

        ## set to NA after the last event in each strata
        index.NA <- which(Lambda0$times>etimes.max[Lambda0$strata+1])
        if(length(index.NA)>0){
            Lambda0$hazard[index.NA] <- NA
            Lambda0$cumhazard[index.NA] <- NA
        }
    }

    ## ** export
    return(Lambda0)
}

