#' Cause-specific Cox proportional hazard regression
#' 
#' Interface for fitting cause-specific Cox proportional hazard regression
#' models in competing risk.
#' 
#' @param formula A list of formulae, one for each cause, each specifying a
#' cause-specific Cox regression model.
#' @param data A data in which to fit the models.
#' @param cause The cause of interest. Defaults to the first cause.
#' @param survtype Either \code{"hazard"} (the default) or \code{"survival"}.
#' If \code{"hazard"} fit cause-specific Cox regression models for all causes.
#' If \code{"survival"} fit one cause-specific Cox regression model for the
#' cause of interest and also a Cox regression model for event-free survival.
#' @param \dots Arguments given to \code{coxph}.
#' @return %% ~Describe the value returned % If it is a LIST, use \item{call
#' }{the call} \item{models }{a list with the fitted (cause-specific) Cox
#' regression objects} \item{response }{the event history response }
#' \item{eventTimes }{the sorted (unique) event times } \item{survtype }{the
#' value of \code{survtype}} \item{theCause }{the cause of interest. see
#' \code{cause}} \item{causes }{the other causes} %% ...
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' 
#' Ulla B. Mogensen \email{ulmo@@biostat.ku.dk}
#' @seealso \code{\link{coxph}}
#' @keywords survival
##' @examples
##' 
##' library(riskRegression)
##' library(prodlim)
##' library(pec)
##' library(survival)
##' data(Melanoma)
##' ## fit two cause-specific Cox models
##' ## different formula for the two causes
##' fit1 <- CSC(list(Hist(time,status)~sex,Hist(time,status)~invasion+epicel+age),
##'             data=Melanoma)
##' print(fit1)
##' fit1a <- CSC(list(Hist(time,status)~sex,Hist(time,status)~invasion+epicel+age),
##'              data=Melanoma,
##'              survtype="surv")
##' print(fit1a)
##' 
##' ## same formula for both causes
##' fit2 <- CSC(Hist(time,status)~invasion+epicel+age,
##'             data=Melanoma)
##' print(fit2)
##' 
##' ## combine a cause-specific Cox regression model for cause 2 
##' ## and a Cox regression model for the event-free survival:
##' ## different formula for cause 2 and event-free survival
##' fit3 <- CSC(list(Hist(time,status)~sex+invasion+epicel+age,
##'                  Hist(time,status)~invasion+epicel+age),
##'             data=Melanoma)
##' print(fit3)
##' 
##' ## same formula for both causes
##' fit4 <- CSC(Hist(time,status)~invasion+epicel+age,
##'             data=Melanoma,
##'             survtype="surv")
##' print(fit4)
##' #'
#' @export
CSC <- function (formula,data,cause,survtype="hazard",...){
    # {{{ type
    survtype <- match.arg(survtype,c("hazard","survival"))
    # }}}
    # {{{ formulae
    if (class(formula)=="formula") formula <- list(formula)
    else if (!length(formula)<=4) stop("Formula can either be a single formula (to be used for both cause specific hazards) or a list of formulae, one for each cause.")
    responseFormula <- reformulate("1", formula[[1]][[2]])
    # }}}
    # {{{ response
    ## require(survival)
    call <- match.call()
    # get information from formula
    mf <- model.frame(responseFormula, data = data, na.action = na.omit)
    response <- model.response(mf)
    time <- response[, "time"]
    status <- response[, "status"]
    event <- getEvent(response)
    # }}}
    # {{{ event times
    eventTimes <- unique(sort(as.numeric(time[status != 0])))
    # }}}
    # {{{ causes
    causes <- getStates(response)
    if (survtype=="hazard")
        NC <- length(causes)
    else
        NC <- 2 # cause of interest and overall survival
    if (length(formula)!=NC && length(formula)>1) stop("Wrong number of formulae. Should be ",NC,".")
    if (length(formula)==1) {
        ## warning("The same formula used for all causes")
        formula <- lapply(1:NC,function(x)formula[[1]])
    }
    # }}}
    # {{{ find the cause of interest
    if (missing(cause)){
        theCause <- causes[1]
    }
    else{
        if ((foundCause <- match(as.character(cause),causes,nomatch=0))==0)
            stop(paste("Requested cause: ",cause," Available causes: ", causes))
        else
            theCause <- foundCause
    }
    otherCauses <- causes[-match(theCause,causes)]
    # }}}
    # {{{ fit Cox models
    if (survtype=="hazard"){
        CoxModels <- lapply(1:NC,function(x){
            if (x==1)
                causeX <- theCause
            else
                causeX <- otherCauses[x-1]
            ## EHF <- prodlim::EventHistory.frame(formula=formula[[x]],
            EHF <- EventHistory.frame(formula=formula[[x]],
                                      data=data,
                                      unspecialsDesign=FALSE,
                                      specials="strata",
                                      specialsDesign=FALSE)
            formulaX <- formula[[x]]
            covData <- cbind(EHF$design,EHF$strata)
            ## response <- model.response(covData)
            time <- as.numeric(EHF$event.history[, "time",drop=TRUE])
            event <- getEvent(EHF$event.history)
            statusX <- as.numeric(event==causeX)
            workData <- data.frame(time=time,status=statusX)
            workData <- cbind(workData,covData)
            formulaXX <- as.formula(paste("survival::Surv(time,status)",as.character(delete.response(terms.formula(formulaX)))[[2]],sep="~"))
            fit <- survival::coxph(formulaXX, data = workData,...)
            ## fit <- rms::cph(formulaXX, data = workData,...)
            ## fit$call$formula <- fit$formula
            fit$call$data <- workData
            ## fit$call$data <- NULL
            fit$call$formula <- formulaXX
            fit
        })
        names(CoxModels) <- paste("Cause",c(theCause,otherCauses))
    }
    else{
        CoxModels <- lapply(1:2,function(x){
            causeX <- theCause
            formulaX <- formula[[x]]
            covData <- model.frame(formulaX,data=data)
            response <- model.response(covData)
            time <- response[, "time"]
            event <- getEvent(response)
            if (x==1){
                statusX <- event==causeX
            }
            else{
                statusX <- response[,"status"]
            }
            workData <- cbind(time=time,status=statusX,covData[,-1,drop=FALSE])
            formulaXX <- as.formula(paste("Surv(time,status)",as.character(delete.response(terms.formula(formulaX)))[[2]],sep="~"))
            fit <- survival::coxph(formulaXX, data = workData,...)
            ## fit <- rms::cph(formulaXX, data = workData,...)
            ## fit$call$formula <- fit$formula
            fit$call$formula <- formulaXX
            fit
        })
        names(CoxModels) <- c(paste("Cause",theCause),"OverallSurvival")
    }
    # }}}
    out <- list(call=call,
                models=CoxModels,
                response=response,
                eventTimes=eventTimes,
                survtype=survtype,
                theCause=theCause,
                causes=c(theCause,otherCauses))
    class(out) <- "CauseSpecificCox"
    out
}


