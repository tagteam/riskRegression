#' Cause-specific Cox proportional hazard regression
#'
#' Interface for fitting cause-specific Cox proportional hazard regression
#' models in competing risk.
#'
#' The causes and their order
#' are determined by \code{prodlim::getStates()} applied to the Hist object.
#' @param formula Either a single \code{Hist} formula or a list of formulas.
#' If it is a list it must contain as many \code{Hist} formulas as there are
#' causes when \code{surv.type="hazard"} and exactly two formulas when \code{surv.type="survival"}.
#' If it is a list the first formula is used for the cause of interest specific Cox regression
#' and the other formula(s) either for the other cause specific Cox regression(s) or for the
#' Cox regression of the combined event where each cause counts as event. Note that when only one
#' formula is given the covariates enter in exactly the same way into all Cox regression analyses.
#' @param data A data in which to fit the models.
#' @param cause The cause of interest. Defaults to the first cause (see Details).
#' @param surv.type Either \code{"hazard"} (the default) or
#' \code{"survival"}.  If \code{"hazard"} fit cause-specific Cox
#' regression models for all causes to predict the event free survival probabilities. If \code{"survival"} fit 
#' a Cox regression model for the hazard function of the combined endpoint (all causes combined with 'or')
#' to predict the event-free survival probabilities.
#' @param fitter Either a single character string specifying the routine to fit all Cox regression models,
#' or a list of character strings specifying (potentially different) routines for each model fitted in the
#' internal loop.
#'
#' Available routines are \code{"coxph"} for \link[survival]{coxph}, \code{"cph"} for \link[rms]{cph},
#' \code{"phreg"} for \link[mets]{phreg}, and \code{"glmnet"} for \link[glmnet]{glmnet}.
#'
#' If \code{fitter} is a list, its length must match the number of models fitted:
#' \itemize{
#' \item \code{surv.type="hazard"}: one model per cause
#' \item \code{surv.type="survival"}: two models (cause of interest + overall survival)
#' }
#' @param fitter_arguments Optional list of per-model argument lists passed to the corresponding fitter.
#' Must be either \code{NULL} (default) or a list whose length matches the number of fitted models (see
#' \code{fitter}). Each element must itself be a list (named or unnamed) of arguments for the fitter used
#' for that model.
#'
#' @param ... Additional arguments passed to the fitter *for all fitted models*. These arguments are
#' appended to each element of \code{fitter_arguments}, i.e., they are used to *augment* the per-model
#' arguments. If the same argument name is provided both in \code{fitter_arguments[[k]]} and in \code{...},
#' the value in \code{...} takes precedence.
#' @return \item{models }{a list with the fitted (cause-specific) Cox
#' regression objects} \item{response }{the event history response }
#' \item{eventTimes }{the sorted (unique) event times } \item{surv.type }{the
#' value of \code{surv.type}} \item{theCause }{the cause of interest. see
#' \code{cause}} \item{causes }{the other causes} %% ...
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk} and Ulla B. Mogensen
#' @references
#' B. Ozenne, A. L. Soerensen, T.H. Scheike, C.T. Torp-Pedersen,
#' and T.A. Gerds. riskregression: Predicting the risk
#' of an event using Cox regression models. R Journal, 9(2):440--460, 2017.
#'
#' J Benichou and Mitchell H Gail. Estimates of absolute cause-specific risk
#' in cohort studies. Biometrics, pages 813--826, 1990.
#'
#' T.A. Gerds, T.H. Scheike, and P.K. Andersen. Absolute risk regression for
#' competing risks: Interpretation, link functions, and prediction. Statistics
#' in Medicine, 31(29):3921--3930, 2012.
#'
#' @seealso \code{\link[survival]{coxph}}
#' @keywords survival
#' @examples
#'
#' library(prodlim)
#' library(survival)
#' data(Melanoma)
#' ## fit two cause-specific Cox models
#' ## different formula for the two causes
#' fit1 <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
#'             data=Melanoma)
#' print(fit1)
#'
#' \dontrun{
#' library(rms)
#' fit1a <- CSC(list(Hist(time,status)~sex+rcs(age,3),Hist(time,status)~invasion+epicel+log(thick)),
#'             data=Melanoma,fitter="cph")
#' print(fit1a)
#' }
#' \dontrun{
#' library(glmnet)
#' # lasso
#' fit1b <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
#'             data=Melanoma,fitter="glmnet")
#' # rigde regression
#' fit1c <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
#'             data=Melanoma,fitter="glmnet")
#' print(fit1b)
#' }
##' \dontrun{
##' library(Publish)
##' publish(fit1)
##' }
##' 
##' ## model hazard of all cause mortality instead of hazard of type 2
##' fit1a <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~invasion+epicel+log(thick)),
##'              data=Melanoma,
##'              surv.type="surv")
##'
##' ## the predicted probabilities are similar
##' plot(predictRisk(fit1,times=500,cause=1,newdata=Melanoma),
##'      predictRisk(fit1a,times=500,cause=1,newdata=Melanoma))
##'
##' ## special case where cause 2 has no covariates
##' fit1b <- CSC(list(Hist(time,status)~sex+age,Hist(time,status)~1),
##'              data=Melanoma)
##' print(fit1b)
##' predict(fit1b,cause=1,times=100,newdata=Melanoma)
##'
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
##'             surv.type="surv",
##'             data=Melanoma)
##' print(fit3)
##'
##' ## same formula for both causes
##' fit4 <- CSC(Hist(time,status)~invasion+epicel+age,
##'             data=Melanoma,
##'             surv.type="surv")
##' print(fit4)
##'
##' ## strata
##' fit5 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
##'             data=Melanoma,
##'             surv.type="surv")
##' print(fit5)
##'
##' ## sanity checks
##'
##' cox1 <- coxph(Surv(time,status==1)~invasion+epicel+age+strata(sex),data=Melanoma)
##' cox2 <- coxph(Surv(time,status!=0)~invasion+epicel+age+strata(sex),data=Melanoma)
##' all.equal(coef(cox1),coef(fit5$models[[1]]))
##' all.equal(coef(cox2),coef(fit5$models[[2]]))
##'
##' ## predictions
##' ##
##' ## surv.type = "hazard": predictions for both causes can be extracted
##' ## from the same fit
##' fit2 <- CSC(Hist(time,status)~invasion+epicel+age, data=Melanoma)
##' predict(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
##' predictRisk(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
##' predictRisk(fit2,cause=2,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
##' predict(fit2,cause=1,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
##' predict(fit2,cause=2,newdata=Melanoma[c(17,99,108),],times=c(100,1000,10000))
##'
##' ## surv.type = "surv" we need to change the cause of interest
##' library(survival)
##' fit5.2 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
##'             data=Melanoma,
##'             surv.type="surv",cause=2)
##' ## now this does not work because the object was fitted with surv.type='surv'
##' try(predictRisk(fit5.2,cause=1,newdata=Melanoma,times=4000))
##'
##' ## but this does
##' predictRisk(fit5.2,cause=2,newdata=Melanoma,times=100)
##' predict(fit5.2,cause=2,newdata=Melanoma,times=100)
##' predict(fit5.2,cause=2,newdata=Melanoma[4,],times=100)
##'
##' fit5.2 <- CSC(Hist(time,status)~invasion+epicel+age+strata(sex),
##'             data=Melanoma,
##'             surv.type="hazard",cause=2)
##'
#' @export
CSC <- function(formula,
                data,
                cause,
                surv.type="hazard",
                fitter="coxph",
                fitter_arguments = NULL,
                ...){
    allowed_fitters <- c("coxph","cph","phreg","glmnet")
    dots_args <- list(...)
    fitter_is_list <- is.list(fitter)
    if (!fitter_is_list){
        fitter <- match.arg(fitter, allowed_fitters)
    } else {
        ## validate basic structure early; full length check happens after nmodels is known
        if (!all(vapply(fitter, function(z) is.character(z) && length(z)==1, logical(1)))){
            stop("If 'fitter' is a list, each element must be a single character string.")
        }
    }
    # {{{ type
    surv.type <- match.arg(surv.type,c("hazard","survival"))
    # }}}
    # {{{ formulae & response
    if (inherits(x=formula,what="formula")) formula <- list(formula)
    call <- match.call()
    # get outcome information from formula
    Rform <- update(formula[[1]],".~1")
    ## mf <- stats::model.frame(Rform, data = data, na.action = na.omit)
    ## response <- stats::model.response(mf)
    response <- eval(Rform[[2]],envir=data)
    if (any(is.na(response)))
        stop("Event history response may not contain missing values")
    time <- response[, "time"]
    status <- response[, "status"]
    event <- prodlim::getEvent(response)
    if ("entry" %in% colnames(response))
        entry <- response[, "entry"]
    else{
        ## phreg requires entry (start time) even if not supplied; keep old behavior
        ## but generalize to the case where 'fitter' is a list
        any_phreg <- if (fitter_is_list) any(vapply(fitter, identical, logical(1), "phreg")) else identical(fitter, "phreg")
        entry <- if (any_phreg) 0 else NULL
    }
    if (any(entry>time)) stop("entry > time detected. Entry time into the study must be strictly greater than outcome time.")
    ## remove event history variables from data
    if(any((this <- match(all.vars(Rform),names(data),nomatch=0))>0)){
        if (data.table::is.data.table(data))
            data <- data[,-this,with=FALSE]
        else
            data <- data[,-this]
    }
    # }}}
    # {{{ sorted unique event times
    eventTimes <- unique(sort(as.numeric(time[status != 0])))
    # }}}
    # {{{ causes
    causes <- prodlim::getStates(response)
    if (surv.type=="hazard")
        NC <- length(causes)
    else
        NC <- 2 # cause of interest and overall survival
    if (length(formula)!=NC[1] && length(formula)>1) stop("Wrong number of formulae. Should be one for each cause ",NC,".")
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
            stop(paste0("Cannot find all requested cause(s) ...\n\n",
                        "Requested cause(s): ", paste0(cause, collapse = ", "),
                        "\n Available causes: ", paste(causes, collapse = ", "),
                        "\n"))
        else{
            theCause <- causes[foundCause]
        }
    }
    otherCauses <- causes[-match(theCause,causes)]
    # }}}
    # {{{ fit Cox models
    if (surv.type=="hazard"){
        nmodels <- NC
    }else {
       nmodels <- 2}
    ## --- normalize fitter to a list of length nmodels ---
    if (fitter_is_list){
        if (length(fitter) != nmodels){
            stop("If 'fitter' is a list, its length must match the number of fitted models (", nmodels, ").")
        }
        fitter_list <- lapply(fitter, function(z) match.arg(z, allowed_fitters))
    } else {
        fitter_list <- rep(list(fitter), nmodels)
    }

    ## --- normalize fitter_arguments to a list of length nmodels ---
    if (is.null(fitter_arguments)){
        fitter_arguments <- vector("list", nmodels)
    } else {
        if (!is.list(fitter_arguments)){
            stop("'fitter_arguments' must be NULL or a list.")
        }
        if (length(fitter_arguments) == 1L && nmodels > 1L){
            fitter_arguments <- rep(fitter_arguments, nmodels)
        }
        if (length(fitter_arguments) != nmodels){
            stop("'fitter_arguments' must have length 1 or match the number of fitted models (", nmodels, ").")
        }
        fitter_arguments <- lapply(fitter_arguments, function(a){
            if (is.null(a)) return(list())
            if (!is.list(a)) stop("Each element of 'fitter_arguments' must be a list (or NULL).")
            if (any("" %in% names(a))) stop("fitter_arguments must be named")
            a
        })
    }

    CoxModels <- lapply(1:nmodels,function(x){
        fitterX <- fitter_list[[x]]
        if (surv.type=="hazard"){
            if (x==1)
                causeX <- theCause
            else
                causeX <- otherCauses[x-1]
        } else{
            causeX <- theCause
        }
        if (surv.type=="hazard"){
            statusX <- as.numeric(event==causeX)
        }else{
            if (x==1){
                statusX <- 1*(event==causeX)
            }
            else{ ## event-free status
                statusX <- response[,"status"]
            }
        }
        if (is.null(entry))
            workData <- data.frame(time=time,status=statusX)
        else
            workData <- data.frame(time=time,status=statusX,entry=entry)
        if(any(this <- match(names(data),names(workData),nomatch=0)>0)){
            warning(paste("Variables named",paste(names(data)[this],collapse=", "),"in data will be ignored."))
            if (is.data.table(data))
                data <- data[,-this,with=FALSE]
            else
                data <- data[,-this,drop=FALSE]
        }
        workData <- cbind(workData,data)
        if (is.null(entry))
            survresponse <- "survival::Surv(time, status)"
        else
            survresponse <- "survival::Surv(entry, time, status)"
        ## check whether right hand side of formula includes ~.
        allvars <- all.vars(formula[[x]])
        if (any(grepl("^\\.$",allvars))){
            formulaXX <- as.formula(paste0(survresponse,"~."))
        }
        else {
            formulaXX <- update(formula[[x]],paste0(survresponse,"~."))
        }
        # previous
        # formulaXX <- update(formula[[x]],paste0(survresponse,"~."))
        ## as.formula(paste(survresponse,
        ## as.character(delete.response(terms.formula(formulaX)))[[2]],
        ## sep="~"))
        args <- list(formulaXX, data = workData)
        ## augment per-model fitter_arguments[[x]] with dots_args (dots take precedence)
        extra.args <- c(dots_args,fitter_arguments[[x]])
        extra.args <- extra.args[!duplicated(names(extra.args))]
        if (fitterX=="coxph"){
            fit <- do.call("coxph",c(args,list(x=TRUE,y=TRUE),extra.args))
        } else if (fitterX=="cph") {
            fit <- do.call("cph",c(args,list(surv=TRUE,x=TRUE,y=TRUE),extra.args))
        } else if (fitterX=="phreg") {
            fit <- do.call("phreg",c(args,extra.args))
        } else if (fitterX=="glmnet"){
            fit <- do.call("GLMnet", c(args, extra.args))
        }
        fit$call$formula <- formulaXX
        fit$call$data <- workData
        fit
    })
    if (surv.type=="hazard"){
        names(CoxModels) <- paste("Cause",c(theCause,otherCauses))
    }else{
        names(CoxModels) <- c(paste("Cause",theCause),"OverallSurvival")
    }
    # }}}
    out <- list(call=call,
                models=CoxModels,
                response=response,
                eventTimes=eventTimes,
                surv.type=surv.type,
                fitter=fitter_list,
                theCause=theCause,
                causes=c(theCause,otherCauses))
    class(out) <- "CauseSpecificCox"
    out
}


