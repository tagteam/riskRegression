## * iidCox - documentation
#' @title Extract iid decomposition from a Cox model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iidCox
#' 
#' @param object object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata [data.frame] Optional new data at which to do iid decomposition
#' @param baseline.iid [logical] Should the influence function for the baseline hazard be computed.
#' @param tau.hazard [numeric vector] the vector of times at which the i.i.d decomposition of the baseline hazard will be computed
#' @param tau.max [numeric] latest time at which the i.i.d decomposition of the baseline hazard will be computed. Alternative to \code{tau.hazard}.
#' @param store.iid [character] the method used to compute the influence function and the standard error.
#' Can be \code{"full"} or \code{"minimal"}. See the details section.
#' @param keep.times [logical] If \code{TRUE} add the evaluation times to the output.
#' @param return.object [logical] If \code{TRUE} return the object where the iid decomposition has been added.
#' Otherwise return a list (see the return section)
#' 
#' @details
#' 
#' This function implements the first three formula (no number,10,11) of the subsection
#' "Empirical estimates" in Ozenne et al. (2017).
#' If there is no event in a strata, the influence function for the baseline hazard is set to 0. 
#'
#' \bold{Argument store.iid}:
#' If \code{n} denotes the sample size, \code{J} the number of jump times, and \code{p} the number of coefficients:
#' \itemize{
#' \item \code{store.iid="full"} exports the influence function for the coefficients and the baseline hazard at each event time.
#' \item \code{store.iid="minimal"} exports the influence function for the coefficients. For the
#' baseline hazard it only computes the quantities necessary to compute the influence function in order to save memory. }
#'
#' More details can be found in appendix B of Ozenne et al. (2017).
#'
#' @references
#' Brice Ozenne, Anne Lyngholm Sorensen, Thomas Scheike, Christian Torp-Pedersen and Thomas Alexander Gerds.
#' riskRegression: Predicting the Risk of an Event using Cox Regression Models.
#' The R Journal (2017) 9:2, pages 440-460.
#' 
#' @return For Cox models, it returns the object with an additional iid slot (i.e. \code{object$iid}).
#' It is a list containing:
#' \itemize{
#'  \item IFbeta: Influence function for the regression coefficient.
#'  \item IFhazard: Time differential of the influence function of the hazard.
#'  \item IFcumhazard: Influence function of the cumulative hazard.
#'  \item calcIFhazard: Elements used to compute the influence function at a given time.
#'  \item time: Times at which the influence function has been evaluated.
#'  \item etime1.min: Time of first event (i.e. jump) in each strata.
#'  \item etime.max: Last observation time (i.e. jump or censoring) in each strata.
#'  \item indexObs: Index of the observation in the original dataset.}
#' 
#' For Cause-Specific Cox models,
#' it returns the object with an additional iid slot for each model (e.g. \code{object$models[[1]]iid}).
#'                 
#' @examples
#' library(survival)
#' library(data.table)
#' library(prodlim)
#' set.seed(10)
#' d <- sampleData(100, outcome = "survival")[,.(eventtime,event,X1,X6)]
#' setkey(d, eventtime)
#' 
#' m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
#' system.time(IF.cox <- iidCox(m.cox))
#' 
#' IF.cox.all <- iidCox(m.cox, tau.hazard = sort(unique(c(7,d$eventtime))))
#' IF.cox.beta <- iidCox(m.cox, baseline.iid = FALSE)
#'
#' @export
`iidCox` <-
    function(object, newdata,
             baseline.iid, tau.hazard, tau.max, store.iid, 
             keep.times, return.object) UseMethod("iidCox")



## * iidCox.coxph
#' @rdname iidCox
#' @export
iidCox.coxph <- function(object, newdata = NULL,
                         baseline.iid = TRUE, tau.hazard = NULL, tau.max = NULL, store.iid = "full", 
                         keep.times = TRUE, return.object = TRUE){

    ## ** extract elements from object
    iInfo <- coxVarCov(object)
    object.modelFrame <- coxModelFrame(object)
    infoVar <- coxVariableName(object, model.frame = object.modelFrame)
    nVar.lp <- length(infoVar$lpvars)
    nObs <- NROW(object.modelFrame)
    
    object.strata <- coxStrata(object, data = NULL, strata.vars = infoVar$stratavars)
    object.levelStrata <- levels(object.strata)
    object.eXb <- exp(coxLP(object, data = NULL, center = FALSE))
    if(nVar.lp>0){
        object.LPdata <- as.matrix(object.modelFrame[,.SD,.SDcols = infoVar$lpvars])
    }else{
        object.LPdata <- NULL
    }
    nStrata <- length(levels(object.strata))

    if(!is.null(object$reverse) && object$reverse){
        reverse <- TRUE
    }else{
        reverse <- FALSE
    }
  
    ## ** Extract new observations
    if(!is.null(newdata)){
    
        if(!inherits(x=newdata,what="data.frame")){
            stop("class of \'newdata\' must inherit from data.frame \n")
        }
    
                                        # if(infoVar$status %in% names(newdata)){ # call Cox model with with event==1
        tempo <- with(newdata, eval(coxFormula(object)[[2]]))
        if(!is.null(object$reverse) && object$reverse){ ## only for reverse product-limit (prodlim object)
            new.status <- 1-tempo[,2]
        }else{
            new.status <- tempo[,2]
        }
        new.time <- tempo[,1]
                                        # }else{ # Cox model from CSC 
                                        #    new.status <- newdata[[infoVar$status]]
                                        #    new.time <- newdata[[infoVar$time]]
                                        # }
    
        new.strata <- coxStrata(object,
                                data = newdata, 
                                sterms = infoVar$strata.sterms,
                                strata.vars = infoVar$stratavars,
                                strata.levels = infoVar$strata.levels)
        new.eXb <- exp(coxLP(object, data = newdata, center = FALSE))
        new.LPdata <- model.matrix(object, data = newdata)
    
    }else{
        new.status <-  object.modelFrame[["status"]]
        new.time <-  object.modelFrame[["stop"]]
        new.strata <-  object.strata
        new.eXb <- object.eXb
        new.LPdata <- object.LPdata
    
    }

    ## ** tests 
    store.iid <- match.arg(store.iid, c("minimal","full"))

    ## time at which the influence function is evaluated
    if(is.list(tau.hazard)){
        if(length(tau.hazard)!=nStrata){
            stop("If a list, argument \"tau.hazard\" must be a list with ",nStrata," elements \n",
                 "each element being the vector of times for each strata \n")
        }

        need.order <- vector(mode = "logical", length = nStrata)
        tau.oorder <- vector(mode = "list", length = nStrata)
        etimes.max <- vector(mode = "numeric", length = nStrata)

        for(iStrata in 1:nStrata){ ## iStrata <- 1
            etimes.max[iStrata] <- attr(tau.hazard[[iStrata]],"etimes.max")
            need.order[iStrata] <- is.unsorted(tau.hazard[[is.unsorted]])  || any(duplicated(tau.hazard[[is.unsorted]]))
            if(need.order){
                Utau.hazard[[iStrata]] <- sort(unique(tau.hazard[[iStrata]]))
            }else{
                Utau.hazard[[iStrata]] <- tau.hazard[[iStrata]]
            }
            tau.oorder[[iStrata]] <- match(tau.hazard[[iStrata]],Utau.hazard[[iStrata]])
        }
    }else if(!is.null(tau.hazard)){
        etimes.max <- attr(tau.hazard,"etimes.max")        

        need.order <- is.unsorted(tau.hazard) || any(duplicated(tau.hazard))
        if(need.order){
            Utau.hazard <- sort(unique(tau.hazard))
        }else{
            Utau.hazard <- tau.hazard
        }
        tau.oorder <- match(tau.hazard,Utau.hazard)
    }else{
        etimes.max <- NULL        
        need.order <- FALSE
    }

    ## ** Compute quantities of interest
  
    ## baseline hazard
    lambda0 <- predictCox(object,
                          type = "hazard",
                          centered = FALSE,
                          keep.strata = TRUE)
    if(is.null(etimes.max)){
        etimes.max <- lambda0$lastEventTime
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
        if(reverse){
            object.order_strata[[iStrata]] <- order(object.modelFrame[object.index_strata[[iStrata]], .SD$stop],
                                                    object.modelFrame[object.index_strata[[iStrata]], .SD$status])
        }else{
            object.order_strata[[iStrata]] <- order(object.modelFrame[object.index_strata[[iStrata]], .SD$stop])
        }

        indexTempo <- object.index_strata[[iStrata]][object.order_strata[[iStrata]]]
        object.eXb_strata[[iStrata]] <- object.eXb[indexTempo]
        if(nVar.lp>0){
            object.LPdata_strata[[iStrata]] <- object.LPdata[indexTempo,,drop = FALSE]
        }else{
            object.LPdata_strata[[iStrata]] <- matrix(nrow = 0, ncol = 0)
        }
        object.status_strata[[iStrata]] <- object.modelFrame[indexTempo, .SD$status]
        object.time_strata[[iStrata]] <- object.modelFrame[indexTempo, .SD$stop]
            
        ## reorder new data
        if(!is.null(newdata)){
            new.index_strata[[iStrata]] <- which(new.strata == object.levelStrata[iStrata])
            if(reverse){
                new.order_strata[[iStrata]] <- order(new.time[new.index_strata[[iStrata]]],
                                                     new.status[new.index_strata[[iStrata]]])
            }else{
                new.order_strata[[iStrata]] <- order(new.time[new.index_strata[[iStrata]]])
            }

            new.eXb_strata[[iStrata]] <- new.eXb[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
            if(nVar.lp>0){
                new.LPdata_strata[[iStrata]] <- new.LPdata[new.index_strata[[iStrata]][new.order_strata[[iStrata]]],,drop = FALSE]
            }else{
                new.LPdata_strata[[iStrata]] <- matrix(nrow = 0, ncol = 0)
            }
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
                                      p = nVar.lp, add0 = TRUE, reverse = reverse)
        
        if(reverse){
            ## when ties, the censored observations correspond to the previous timepoint
            new.indexJump[[iStrata]] <- prodlim::sindex(Ecpp[[iStrata]]$Utime1, new.time) - 1 - reverse * (1-new.status) * (new.time %in% Ecpp[[iStrata]]$Utime1)
        }else{
            new.indexJump[[iStrata]] <- prodlim::sindex(Ecpp[[iStrata]]$Utime1, new.time) - 1
        }
        ## print(head(cbind(time = new.time, status = new.status, jump = new.indexJump[[iStrata]])))
        ## print(new.indexJump[[iStrata]])

        ## if event/censoring is before the first event in the training dataset 
        ## then sindex return 0 thus indexJump is -1
        ## the following 3 lines convert -1 to 0
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

    ## ** Prepare output
    out <- list(IFbeta = matrix(NA, nrow = nObs, ncol = nVar.lp,
                                dimnames = list(NULL, infoVar$lpvars)),  
                IFhazard = NULL,
                IFcumhazard = NULL,
                calcIFhazard = list(delta_iS0 = NULL,
                                    Elambda0 = NULL,
                                    cumElambda0 = NULL,
                                    lambda0_iS0= NULL,
                                    cumLambda0_iS0= NULL,
                                    time1 = NULL),
                obstime = new.time,
                time = vector(mode = "list", length = nStrata),  # time at which the IF is assessed
                etime1.min = rep(NA, nStrata),
                etime.max = etimes.max,
                indexObs = new.order,
                store.iid = store.iid
                )

    if(store.iid=="minimal"){
        out$calcIFhazard$Elambda0 <- vector(mode = "list", length = nStrata)
        out$calcIFhazard$cumElambda0 <- vector(mode = "list", length = nStrata)
        out$calcIFhazard$eXb <- matrix(NA, nrow = nObs, ncol = nStrata)
        out$calcIFhazard$lambda0_iS0 <- vector(mode = "list", length = nStrata)
        out$calcIFhazard$cumLambda0_iS0 <- vector(mode = "list", length = nStrata)
        out$calcIFhazard$delta_iS0 <- matrix(NA, nrow = nObs, ncol = nStrata)
        out$calcIFhazard$time1 <- vector(mode = "list", length = nStrata)
    }else{
        out$IFhazard <- vector(mode = "list", length = nStrata)
        out$IFcumhazard <- vector(mode = "list", length = nStrata)
    }

    ## ** Computation of the influence function (coefficients)
    for(iStrata in 1:nStrata){

        iSubset <- new.index_strata[[iStrata]]
        iOrder <- new.order_strata[[iStrata]]
        iOrderBack <- order(new.order_strata[[iStrata]])
        new.indexJump_strata <- new.indexJump[[iStrata]][iSubset[iOrder]]
    
        ## IF
        if(nVar.lp > 0){
            out$IFbeta[iSubset,] <- IFbeta_cpp(newT = new.time_strata[[iStrata]],
                                               neweXb = new.eXb_strata[[iStrata]],
                                               newX = new.LPdata_strata[[iStrata]],
                                               newStatus = new.status_strata[[iStrata]], 
                                               newIndexJump = new.indexJump_strata, 
                                               S01 = Ecpp[[iStrata]]$S0,
                                               E1 = Ecpp[[iStrata]]$E,
                                               time1 = Ecpp[[iStrata]]$Utime1,
                                               iInfo = iInfo,
                                               p = nVar.lp)[iOrderBack,,drop=FALSE]
        }else{
            out$IFbeta[iSubset,] <- matrix(NA, ncol = 1, nrow = length(new.index_strata[[iStrata]]))
        }
        
    }
    
    # }}}
    
    ## ** Computation of the influence function (baseline hazard)
    if(baseline.iid){
        for(iStrata in 1:nStrata){

            ## hazard
            if(nStrata==1){ # select only the time,lambda corresponding to the events and not censored observations
                timeStrata <- lambda0$time[lambda0$time %in% Ecpp[[1]]$Utime1]
                lambda0Strata <- lambda0$hazard[lambda0$time %in% Ecpp[[1]]$Utime1]
            }else{ # same within the strata
                index.strata <- which(lambda0$strata == object.levelStrata[iStrata])
                index.keep <- index.strata[lambda0$time[index.strata] %in% Ecpp[[iStrata]]$Utime1]
                
                timeStrata <- lambda0$time[index.keep]
                lambda0Strata <- lambda0$hazard[index.keep]
            }
        
            if(length(timeStrata)>0){
                out$etime1.min[iStrata] <- timeStrata[1]
            }else{ ## case of no event in strata
                out$etime1.min[iStrata] <- min(Ecpp[[iStrata]]$Utime1)
            }
            ## tau.hazard
            if(is.list(tau.hazard)){
                iTau.oorder <- tau.oorder[[iStrata]]
                tau.hazard_strata <- Utau.hazard[[nStrata]]
            }else if(!is.null(tau.hazard)){
                iTau.oorder <- tau.oorder
                tau.hazard_strata <- Utau.hazard
            }else if(any(object.status_strata[[iStrata]]==1)){
                tau.hazard_strata <- unique(object.time_strata[[iStrata]][object.status_strata[[iStrata]] == 1])
                if(!is.null(tau.max)){
                    if(any(tau.hazard_strata<=tau.max)){
                        tau.hazard_strata <- tau.hazard_strata[tau.hazard_strata<=tau.max]
                    }else{
                        tau.hazard_strata <- 0
                    }
                }
            }else{
                tau.hazard_strata <- 0
            }
            

            ## E
            nUtime1_strata <- length(Ecpp[[iStrata]]$Utime1)
            if(nVar.lp > 0){
                Etempo <- Ecpp[[iStrata]]$E[-NROW(Ecpp[[iStrata]]$E),,drop = FALSE]
            }else{
                Etempo <- matrix(0, ncol = 1, nrow = nUtime1_strata-1)
            }

            ## IF
            if(length(timeStrata)==0){
                ## no event in the strata
                if(store.iid=="full"){
                    IFlambda_res <- list(hazard = matrix(c(0,NA)[(tau.hazard_strata > etimes.max[iStrata])+1], byrow = TRUE, nrow = nObs, ncol = length(tau.hazard_strata)),
                                         cumhazard = matrix(c(0,NA)[(tau.hazard_strata > etimes.max[iStrata])+1], byrow = TRUE, nrow = nObs, ncol = length(tau.hazard_strata)))
                }else if(store.iid=="minimal"){
                    IFlambda_res <- list(delta_iS0 = rep(0, nObs),
                                         Elambda0 = matrix(0, nrow = nVar.lp, ncol = length(tau.hazard_strata)),
                                         cumElambda0 = matrix(0, nrow = nVar.lp, ncol = length(tau.hazard_strata)),
                                         eXb = rep(0, nObs),
                                         lambda0_iS0 = rep(0, length(timeStrata)),
                                         cumLambda0_iS0 = rep(0, length(timeStrata)),
                                         time1 = timeStrata)
                }
            }else{
                IFlambda_res <- IFlambda0_cpp(tau = tau.hazard_strata,
                                              IFbeta = out$IFbeta,
                                              newT = new.time, neweXb = new.eXb, newStatus = new.status, newIndexJump = new.indexJump[[iStrata]], newStrata = as.numeric(new.strata),
                                              S01 = Ecpp[[iStrata]]$S0,
                                              E1 = Etempo,
                                              time1 = timeStrata, lastTime1 = etimes.max[iStrata],
                                              lambda0 = lambda0Strata,
                                              p = nVar.lp, strata = iStrata,
                                              minimalExport = (store.iid=="minimal"),
                                              reverse = reverse)
            }
            ## output
            if(need.order){
                out$time[[iStrata]] <- tau.hazard_strata[iTau.oorder]
            }else{
                out$time[[iStrata]] <- tau.hazard_strata
            }
            if(length(out$time[[iStrata]])==0) browser()
            if(store.iid=="minimal"){
                if(need.order && nVar.lp>0){
                    out$calcIFhazard$Elambda0[[iStrata]] <- IFlambda_res$Elambda0[,iTau.oorder,drop=FALSE]
                    out$calcIFhazard$cumElambda0[[iStrata]] <- IFlambda_res$cumElambda0[,iTau.oorder,drop=FALSE]
                }else{
                    out$calcIFhazard$Elambda0[[iStrata]] <- IFlambda_res$Elambda0
                    out$calcIFhazard$cumElambda0[[iStrata]] <- IFlambda_res$cumElambda0
                }
                out$calcIFhazard$eXb[,iStrata] <- IFlambda_res$eXb
                out$calcIFhazard$lambda0_iS0[[iStrata]] <- c(0,IFlambda_res$lambda0_iS0)
                out$calcIFhazard$cumLambda0_iS0[[iStrata]] <- c(0,IFlambda_res$cumLambda0_iS0)
                out$calcIFhazard$delta_iS0[,iStrata] <- IFlambda_res$delta_iS0
                out$calcIFhazard$time1[[iStrata]] <- c(0,IFlambda_res$time1) # event time by strata
            }else{
                if(keep.times){
                    colnames(IFlambda_res$hazard) <- tau.hazard_strata
                    colnames(IFlambda_res$cumhazard) <- tau.hazard_strata
                } 
                if(need.order){
                    out$IFhazard[[iStrata]] <- IFlambda_res$hazard[,iTau.oorder,drop=FALSE]
                    out$IFcumhazard[[iStrata]] <- IFlambda_res$cumhazard[,iTau.oorder,drop=FALSE]
                }else{
                    out$IFhazard[[iStrata]] <- IFlambda_res$hazard
                    out$IFcumhazard[[iStrata]] <- IFlambda_res$cumhazard
                }
            }
        }
    }

    ## ** export
    if(return.object){
        object$iid <- out
        return(object)
    }else{
        return(out)
    }
}

## * iidCox.cph
#' @rdname iidCox
#' @export
iidCox.cph <- iidCox.coxph

## * iidCox.phreg
#' @rdname iidCox
#' @export
iidCox.phreg <- iidCox.coxph

## * iidCox.prodlim
#' @rdname iidCox
#' @export
iidCox.prodlim <- iidCox.coxph 

## * iidCox.CauseSpecificCox
#' @rdname iidCox
#' @export
iidCox.CauseSpecificCox <- function(object, newdata = NULL,
                                    baseline.iid = TRUE, tau.hazard = NULL, tau.max = NULL, store.iid = "full", 
                                    keep.times = TRUE, return.object = TRUE){

    nCause <- length(object$causes)
    
    for(iCause in 1:nCause){
        object$models[[iCause]]$iid <- iidCox(object$models[[iCause]],
                                              newdata = newdata,
                                              baseline.iid = baseline.iid,
                                              tau.hazard = tau.hazard,
                                              tau.max = tau.max,
                                              store.iid = store.iid,
                                              keep.times = keep.times,
                                              return.object = FALSE)
    }

    if(return.object){
        return(object)
    }else{
        return(lapply(object$models, function(iM){iM$iid}))
    }
}


## * is.iidCox
##' @title Check Computation of the Influence Function in a Cox Model
##' @description Check whether the influence function of the Cox model or cause specific Cox models has been stored in the object.
##' @param object fitted Cox regression model object either obtained with \code{coxph} (survival package), \code{cph} (rms package), \code{CSC} (riskRegression package).
##' @export
`is.iidCox` <- function(object) UseMethod("is.iidCox")

##' @export
is.iidCox.default <- function(object){
    return(NA)
}

##' @export
is.iidCox.coxph <- function(object){
    return(!is.null(object$iid))
}
##' @export
is.iidCox.cph <- is.iidCox.coxph
##' @export
is.iidCox.phreg <- is.iidCox.coxph
##' @export
is.iidCox.prodlim <- is.iidCox.coxph
##' @export
is.iidCox.CauseSpecificCox <- function(object){
    out <- all(unlist(lapply(object$models, function(iM){!is.null(iM$iid)})))
    return(out)
}
