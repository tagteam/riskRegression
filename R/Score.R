# {{{ roxy header
##' Methods to score the predictive performance of risk markers and risk prediction models
##'
##' The function implements a toolbox for the risk prediction modeller:
##' all tools work for the three outcomes: (1) binary (uncensored),
##' (2) right censored time to event without competing risks,
##' (3) right censored time to event with competing risks
##' 
##' Computed are the (time-dependent) Brier score and the (time-dependent)
##' area under the ROC curve for a list of risk prediction models either in
##' external validation data or in the learning data using bootstrap
##' cross-validation. The function optionally provides results for plotting (time-point specific)
##' ROC curves, for (time-point specific) calibration curves and for (time-point specific) retrospective boxplots.
##'
##' For uncensored binary outcome the Delong-Delong test is used to contrast AUC of rival models in all other cases
##' the p-values correspond to Wald tests based on standard errors obtained with an estimate of the influence function
##' as described in detain in the appendix of Blanche et al. (2015).
##' @title Score risk predictions
##' @aliases Score Score.list
##' @param object List of risk predictions (see details and examples).
##' @param formula A formula which identifies the outcome (left hand
##'     side). E.g., \code{Y ~ 1} for binary and \code{Hist(time,status) ~ 1} for time-to-event outcome.
##' In right censored data, the right hand side of the
##'     formula is used to estimate the inverse probability of censoring weights (IPCW) model.
##' @param data \code{data.frame} or \code{data.table} in which the formula can be
##'     interpreted. 
##' @param metrics Character vector specifying which metrics to
##'     apply. Case does not matter. Choices are \code{"AUC"} and \code{"Brier"}.
##' @param summary Character vector specifying which summary
##'     statistics to apply to the predicted risks. The only choice is 
##'     \code{"riskQuantile"} which enables (time-point specific) boxplots of
##'     predicted risks conditional on the outcome (at the time-point).
##'     Set to \code{NULL} to avoid estimation of retrospective risk quantiles.
##' @param plots Character vector specifying which plots to prepare.
##'     Currently implemented are \code{"ROC"}, \code{"Calibration"} and \code{"boxplot"}
##' @param cause Event of interest. Used for binary outcome \code{Y}
##'     to specify that risks are risks of the event \code{Y=event}
##'     and for competing risks outcome to specify the cause of
##'     interest.
##' @param times For survival and competing risks outcome: list of
##'     prediction horizons. All times which are greater than the
##'     maximal observed time in the data set are automatically removed.
##'     Note that the object returned by the function may become huge when
##'     the prediction performance is estimated at many prediction horizons.
##' @param landmarks Not yet implemented.
##' @param use.event.times If \code{TRUE} merge all unique event times with 
##'     the vector given by argument \code{times}.
##' @param null.model If \code{TRUE} fit a risk prediction model which ignores
##'     the covariates and predicts the same value for all subjects. The model is fitted using \code{data}
##' and the left hand side of \code{formula}. For binary outcome this is just the empirical prevalence. For (right censored) time to event outcome, the null models are
##' equal to the Kaplan-Meier estimator (no competing risks) and the Aalen-Johansen estimator (with competing risks).
##' @param conf.int Either logical or a numeric value between 0 and 1. In right censored data,
##'     confidence intervals are based on Blanche et al (see references). Setting \code{FALSE} prevents the 
##'     computation confidence intervals. \code{TRUE} means compute 95 percent confidence 
##'     intervals and corresponding p-values for AUC and Brier score. If set to 0.87, the 
##'     level of significance is 13 percent. So, do not set it to 0.87.
##' @param contrasts Either logical or a list of contrasts. A list of contrasts defines which risk prediction models (markers)
##'    should be contrasted with respect to their prediction performance.
##'   If \code{TRUE} do all possible comparisons. For
##'     example, when \code{object} is a list with two risk prediction models and
##'     \code{null.model=TRUE} setting \code{TRUE} is equivalent to
##'     \code{list(c(0,1,2),c(1,2))} where \code{c(0,1,2)} codes for the
##'     two comparisons: 1 vs 0 and 2 vs 0 (positive integers refer to
##'     elements of \code{object}, 0 refers to the benchmark null
##'     model which ignores the covariates).  This again is equivalent
##'     to explicitly setting \code{list(c(0,1),c(0,2),c(1,2))}. A
##'     more complex example: Suppose \code{object} has 7 elements and you
##'     want to do the following 3 comparisons: 6 vs 3, 2 vs 5 and 2
##'     vs 3, you should set \code{contrasts=c(6,3),c(2,5,3)}.
##' @param probs Quantiles for retrospective summary statistics of the
##'     predicted risks. This affects the result of the function \code{boxplot.Score}.
##' @param cens.method Method for dealing with right censored
##'     data. Either \code{"outside.ipcw"} or \code{"inside.ipcw"}. \code{"outside.ipcw"} 
##'     Here IPCW refers to inverse probability of censoring weights. The attributes inside and outside are
##'      only relevant for cross-validation (\code{split.method}) where the inside method computes the IPCW model
##' in the the loop in each of the current test data sets whereas the outside method computes the IPCW model only once in the full data.
##'     Also, jackknife pseudo-values may become another option for all statistics. Right now
##'     they are only used for calibration curves.
##' @param cens.model Model for estimating inverse probability of
##'     censored weights. Implemented are the Kaplan-Meier method (\code{"km"}) and
##' Cox regression (\code{"cox"}) both applied to the censored times. If the right hand side of \code{formula} does not specify covariates,
##' the Kaplan-Meier method is used even if this argument is set to \code{"cox"}.
##' @param split.method Method for cross-validation. Right now the only choice is \code{bootcv} in which case bootstrap learning sets
##' are drawn with our without replacement (argument \code{M}) from \code{data}. The data not included in the current bootstrap learning
##' set are used as validation set to compute the prediction performance.
##' @param B Number of bootstrap sets for cross-validation.
##' @param M Size of subsamples for bootstrap cross-validation. If specified it
##'     has to be an integer smaller than the size of \code{data}.
##' @param seed Super seed for setting training data seeds when
##'     randomly splitting (bootstrapping) the data during cross-validation.
##' @param trainseeds Seeds for training models during cross-validation.
##' @param keep list of characters (not case sensitive) which determines additional output.
##' \code{"residuals"} provides Brier score residuals and
##' \code{"splitindex"} provides sampling index used to split the data into training and validation sets.  
##' @param predictRisk.args
##'  A list of argument-lists to control how risks are predicted. 
##'  The names of the lists should be the S3-classes of the \code{object}.
##'  The argument-lists are then passed on to the S3-class specific predictRisk method. 
##'  For example, if your object contains one or several random forest model fitted with the function randomForestSRC::rfsrc then you can
##'  specify additional arguments for the function riskRegression::predictRisk.rfsrc which will pass
##'  these on to the function randomForestSRC::predict.rfsrc. A specific example in this case would be
##'  \code{list(rfsrc=list(na.action="na.impute"))}.
##'
##'  A more flexible approach is to write a new predictRisk S3-method. See Details.
##' @param ... Named list containging additional arguments that are passed on to the \code{predictRisk} methods corresponding to object. See examples.
##' @return List with scores and assessments of contrasts, i.e.,
##'     tests and confidence limits for performance and difference in performance (AUC and Brier),
##'     summaries and plots. Most elements are in\code{data.table} format.
##' @details
##' The already existing predictRisk methods (see methods(predictRisk)) may not cover all models and methods
##' for predicting risks. But users can quickly extend the package as explained in detail in Mogensen et al. (2012) for
##' the predecessors \code{pec::predictSurvProb} and \code{pec::predictEventProb} which have been unified as
##' \code{riskRegression::predictRisk}.
##' 
##' @examples
##' # binary outcome
##' library(lava)
##' set.seed(18)
##' learndat <- sampleData(48,outcome="binary")
##' testdat <- sampleData(40,outcome="binary")
##'
##' # score logistic regression models
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
##' lr2 = glm(Y~X3+X5,data=learndat,family=binomial)
##' Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5)"=lr2),formula=Y~1,data=testdat)
##'
##' xb=Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,
##'         data=testdat,plots=c("calibration","ROC"))
##' plotROC(xb)
##' plotCalibration(xb)
##' 
##' # compute AUC for a list of continuous markers
##' markers = as.list(testdat[,.(X6,X7,X8,X9,X10)])
##' u=Score(markers,formula=Y~1,data=testdat,metrics=c("auc"))
##'
##' # cross-validation
##' \dontrun{
##' lr1a = glm(Y~X6,data=learndat,family=binomial)
##' lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
##' Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=3)
##'}
##' # survival outcome
##' 
##' # Score Cox regression models
##' library(survival)
##' library(rms)
##' library(prodlim)
##' set.seed(18)
##' trainSurv <- sampleData(100,outcome="survival")
##' testSurv <- sampleData(40,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
##' xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,times=c(5,8))
##' xs
##' 
##' # time-dependent AUC for list of markers
##' survmarkers = as.list(testSurv[,.(X6,X7,X8,X9,X10)])
##' Score(survmarkers,
##'       formula=Surv(time,event)~1,metrics="auc",data=testSurv,
##'       conf.int=TRUE,times=c(5,8))
##' 
##' # compare models on test data
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=testSurv,conf.int=TRUE,times=c(5,8))
##'
##' # crossvalidation models in traindata
##' \dontrun{
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
##'       split.method="bootcv",B=3)
##' }
##'
##' # restrict number of comparisons
##'\dontrun{
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=trainSurv,contrasts=TRUE,
##' null.model=FALSE,conf.int=TRUE,times=c(5,8),split.method="bootcv",B=3)
##'}
##' # competing risks outcome
##' set.seed(18)
##' trainCR <- sampleData(40,outcome="competing.risks")
##' testCR <- sampleData(40,outcome="competing.risks")
##' library(riskRegression)
##' library(cmprsk)
##' # Cause-specific Cox regression
##' csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR)
##' csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR)
##' # Fine-Gray regression
##' fgr1 = FGR(Hist(time,event)~X1+X2+X7+X9,data=trainCR,cause=1)
##' fgr2 = FGR(Hist(time,event)~X3+X5+X6,data=trainCR,cause=1)
##' Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2,
##'            "FGR(X1+X2+X7+X9)"=fgr1,"FGR(X3+X5+X6)"=fgr2),
##'       formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(5,8))
##' 
##' @author Thomas A Gerds \email{tag@@biostat.ku.dk} and Paul Blanche \email{paul.blanche@@univ-ubs.fr}
##' @references
##'
##' Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
##' Evaluating Random Forests for Survival Analysis Using Prediction Error
##' Curves. Journal of Statistical Software, 50(11), 1-23. URL
##' http://www.jstatsoft.org/v50/i11/.
##'
##' Paul Blanche, Cecile Proust-Lima, Lucie Loubere, Claudine Berr, Jean- Francois Dartigues, and
##' Helene Jacqmin-Gadda. Quantifying and comparing
##' dynamic predictive accuracy of joint models for longitudinal marker and
##' time-to-event in presence of censoring and competing risks. Biometrics, 71
##' (1):102--113, 2015.
##' P. Blanche, J-F Dartigues, and H. Jacqmin-Gadda. Estimating and comparing
##'     time-dependent areas under receiver operating characteristic curves for
##'     censored event times with competing risks. Statistics in Medicine, 32(30):
##'     5381--5397, 2013.
##' 
#' E. Graf et al.  (1999), Assessment and comparison of prognostic
#' classification schemes for survival data. Statistics in Medicine, vol 18,
#' pp= 2529--2545.
#' 
#' Efron, Tibshirani (1997) Journal of the American Statistical Association 92,
#' 548--560 Improvement On Cross-Validation: The .632+ Bootstrap Method.
#' 
#' Gerds, Schumacher (2006), Consistent estimation of the expected Brier score
#' in general survival models with right-censored event times. Biometrical
#' Journal, vol 48, 1029--1040.
#' 
#' Thomas A. Gerds, Martin Schumacher (2007) Efron-Type Measures of Prediction
#' Error for Survival Analysis Biometrics, 63(4), 1283--1287
#' doi:10.1111/j.1541-0420.2007.00832.x
#' 
#' Martin Schumacher, Harald Binder, and Thomas Gerds. Assessment of survival
#' prediction models based on microarray data. Bioinformatics, 23(14):1768-74,
#' 2007.
#' 
#' Mark A. van de Wiel, Johannes Berkhof, and Wessel N. van Wieringen Testing
#' the prediction error difference between 2 predictors Biostatistics (2009)
#' 10(3): 550-560 doi:10.1093/biostatistics/kxp011
##'
##' @export
##'
# }}}
# {{{
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
Score.list <- function(object,
                       formula,
                       data,
                       metrics=c("auc","brier"),
                       summary=NULL, # riskQuantile, risk
                       plots= NULL, # c("roc","calibration","boxplot","p-values"),
                       cause,
                       times,
                       landmarks,
                       use.event.times=FALSE,
                       null.model=TRUE,
                       conf.int=.95,
                       contrasts=TRUE,
                       probs=c(0,0.25,0.5,0.75,1),
                       cens.method="outside.ipcw",
                       cens.model="cox",
                       split.method,
                       B,
                       M,
                       seed,
                       trainseeds,
                       keep,
                       predictRisk.args,
                       ...){
    id=time=status=id=WTi=b=time=status=model=reference=p=model=pseudovalue=ReSpOnSe=event=NULL
    # }}}
    theCall <- match.call()
    # ----------------------------find metrics and plots ----------------------
    # {{{
    ## Metrics <- lapply(metrics,grep,c("AUC","Brier"),ignore.case=TRUE,value=TRUE)
    plots[grep("^box",plots,ignore.case=TRUE)] <- "boxplot"
    summary[grep("^riskQuantile",summary,ignore.case=TRUE)] <- "riskQuantile"
    # force riskQuantile to make boxplots
    if (("boxplot"%in% plots) && !("riskQuantile" %in% summary))
        summary <- c(summary,"riskQuantile")
    summary[grep("^risk$|^risks$",summary,ignore.case=TRUE)] <- "risks"
    metrics[grep("^auc$",metrics,ignore.case=TRUE)] <- "AUC"
    metrics[grep("^brier$",metrics,ignore.case=TRUE)] <- "Brier"
    plots[grep("^roc$",plots,ignore.case=TRUE)] <- "ROC"
    plots[grep("^cal",plots,ignore.case=TRUE)] <- "Calibration"
    if (length(posRR <- grep("^rr$|^r2|rsq",summary,ignore.case=TRUE))>0){
        if (!null.model) stop("Need the null model to compute R^2 but argument 'null.model' is FALSE.")
        summary <- summary[-posRR]
        if (!("Brier" %in% metrics)) metrics <- c(metrics,"Brier")
        if (!null.model) {
            null.model <- TRUE 
            warning("Value of argument 'null.model' ignored as the null model is needed to compute R^2.")
        }
        rsquared <- TRUE
    }else{
        rsquared <- FALSE
    }
    ## Plots <- lapply(plots,grep,c("Roc","Cal"),ignore.case=TRUE,value=TRUE)
    if ("ROC" %in% plots) {
        ## add AUC if needed
        if (!("AUC" %in% metrics)) metrics <- c(metrics,"AUC")
    }
    if ("Calibration" %in% plots) {
        ## add pseduo if needed
        if (!("pseudo" %in% cens.method)) cens.method <- c(cens.method,"pseudo")
    }
    # }}}
    # -----------------parse other arguments and prepare data---------
    # {{{ censoring model arguments
    if (length(grep("^km|^kaplan|^marg",cens.model,ignore.case=TRUE))>0)
        cens.model <- "KaplanMeier"
    else cens.model <- "cox"
    # }}}
    # {{{ Response
    if (missing(formula)){stop("Argument formula is missing.")}    
    formula.names <- try(all.names(formula),silent=TRUE)
    if (!(formula.names[1]=="~")
        ||
        (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
        stop("Invalid specification of formula.\n Could be that you forgot the right hand side:\n ~covariate1 + covariate2 + ...?\nNote that any subsetting, ie data$var or data[,\"var\"], is not supported by this function.")
    }
    if (missing(data)){stop("Argument data is missing.")}
    data <- data.table(data)
    responseFormula <- stats::update(formula,~1)
    ## if (missing(event)) event <- 1
    responsevars <- all.vars(responseFormula)
    if (!all(responsevars%in%names(data)))stop("Response variable(s) ",paste(responsevars,collapse=", ")," not found in data.")
    response <- getResponse(formula=responseFormula,
                            data=data,
                            cause=cause,
                            vars=responsevars)
    response.dim <- NCOL(response)
    responseType <- attr(response,"model")
    states <- attr(response,"states")
    if (missing(cause)) 
        if (responseType=="binary")
            cause <- attr(response,"event")
        else
            cause <- states[[1]]
    # add null model and find names for the object
    if (null.model==TRUE){
        nullobject <- getNullModel(formula=formula,data=data,responseType=responseType)
    } else{
        nullobject <- NULL
    }
    ## put ReSpOnSe for binary and (time, event, status) in the first column(s) 
    ## data[,eval(responsevars):=NULL]
    data <- cbind(response,data)
    if (responseType=="binary")
        data[,event:=factor(ReSpOnSe,levels=1:(length(states)+1),labels=c(states,"censored"))]
    if (responseType=="survival")
        formula <- stats::update(formula,"Hist(time,status)~.")
    if (responseType=="competing.risks")
        formula <- stats::update(formula,"Hist(time,event)~.")
    N <- NROW(response)
    ## stop("Dont know how to predict response of type ",responseType))
    censType <- attr(response,"cens.type")
    if (is.null(censType)) censType <- "uncensoredData"
    # }}}
    # {{{ SplitMethod
    if (!missing(seed)) set.seed(seed)
    split.method <- getSplitMethod(split.method=split.method,B=B,N=N,M=M)
    B <- split.method$B
    splitIndex <- split.method$index
    do.resample <- !(is.null(splitIndex))
    # }}}
    # {{{ Checking the ability of the elements of object to predict risks
    # {{{ number of models and their labels
    NF <- length(object)
    # }}}
    allmethods <- utils::methods(predictRisk)
    if (is.null(names(object))){
        names(object) <- sapply(object,function(o)class(o)[1])}
    else {
            names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
        }
    names.object <- names(object) <- make.unique(names(object))
    ## sanity checks
    object.classes <- lapply(object,function(x)class(x))
    lapply(1:NF,function(f){
        name <- names.object[[f]]
        if (any(c("integer","factor","numeric","matrix") %in% object.classes[[f]])){
            if (split.method$internal.name!="noPlan")
                stop(paste0("Cannot crossvalidate performance of deterministic risk predictions:",name))
        }
        if (c("factor") %in% object.classes[[f]]){
            if (length(grep("^brier$",metrics,ignore.case=TRUE)>0) || length(grep("^cal",plots,ignore.case=TRUE)>0)){
                stop(paste0("Cannot compute Brier score or calibration plots for predictions that are coded as factor: ",name))
            }
        }
        candidateMethods <- paste("predictRisk",unlist(object.classes[[f]]),sep=".")
        if (all(match(candidateMethods,allmethods,nomatch=0)==0)){
            stop(paste("Cannot find function (S3-method) called ",
                       candidateMethods,
                       sep=""))
        }
    })
    # }}}
    # {{{ additional arguments for predictRisk methods
    if (!missing(predictRisk.args)){
        if (!(all(names(predictRisk.args) %in% unlist(object.classes))))
            stop(paste0("Argument predictRisk.args should be a list whose names match the S3-classes of the argument object.
  For example, if your object contains a random forest model fitted with the function randomForestSRC::rfsrc then you can
  specify additional arguments for the function riskRegression::predictRisk.rfsrc which will pass
  these on to the function randomForestSRC::predict.rfsrc. A specific example in this case would be
  predictRisk.args=list(\"rfsrc\"=list(na.action = \"na.impute\"). \n\nThe classes of your current object are: ",
  paste(unlist(object.classes),collapse=", ")))
    }else{
        predictRisk.args <- NULL
    }
    # }}}
    # {{{ add null model and check resampling ability
    if (!is.null(nullobject)) {
        mlevs <- 0:NF
        mlabels <- c(names(nullobject),names(object))
    } else{
        mlevs <- 1:NF
        mlabels <- names(object)
    }
    if (do.resample){
        nix <- lapply(1:length(object),function(f){
            fit <- object[[f]]
            if(is.null(fit$call))
                stop(paste("model",names(object)[f],"does not have a call argument."))
        })
    }
    # }}}
    # {{{ resolve keep statements
    if (!missing(keep) && is.character(keep)){
        if("residuals" %in% tolower(keep)) keep.residuals=TRUE else keep.residuals = FALSE
        if ("splitindex" %in% tolower(keep)) keep.splitindex=TRUE else keep.splitindex = FALSE
        if ("cv" %in% tolower(keep)) keep.cv=TRUE else keep.cv = FALSE
    }else{
        keep.residuals=FALSE
        keep.cv=FALSE
        keep.splitindex=FALSE
    }
    # }}}
    # {{{ resolve se.fit and contrasts
    if (is.logical(conf.int) && conf.int==FALSE
        || conf.int<=0
        || conf.int>1) 
        se.fit <- 0L
    else
        se.fit <- 1L
    if (se.fit==1L) {
        if (is.numeric(conf.int) && conf.int<1 && conf.int>0)
            alpha <- 1-conf.int
        else
            alpha <- .05
    }else{
        alpha <- NA
    }
    if ((NF+length(nullobject))<=1) dolist <- NULL 
    else{
        if (se.fit==0L || (is.logical(contrasts) && contrasts==FALSE)){
            dolist <- NULL
        } else{
            if (is.logical(contrasts) && contrasts==TRUE){
                if (is.null(nullobject)){
                    dolist <- lapply(1:(NF-1),function(i){c(i,(i+1):NF)})
                } else{
                    dolist <- lapply(0:(NF-1),function(i){c(i,(i+1):NF)})
                }
            }else{
                dolist <- contrasts
                if (!is.list(contrasts)) contrasts <- list(contrasts)
                if (!(all(sapply(dolist,function(x){x<=NF && x>=0}))))
                    stop(paste("Argument contrasts should be a list of positive integers possibly mixed with 0 that refer to elements of object.\nThe object has ",NF,"elements but "))
            }
        }
    }
    # }}}
    # {{{ add id *before* ordering the data. data have to be ordered when ipcw is called
    ID <- 1:N
    data[,ID:=ID]
    # }}}
    # {{{ Evaluation landmarks and horizons (times)
    if (responseType %in% c("survival","competing.risks")){
        ## in case of a tie, events are earlier than right censored
        data.table::setorder(data,time,-status)
        eventTimes <- unique(data[,time])
        maxtime <- eventTimes[length(eventTimes)]
        ## print(maxtime)
        include.times <- NULL
        if (missing(landmarks)){
            start <- 0
            if (missing(times)){
                if (use.event.times==TRUE)
                    times <- unique(c(start,eventTimes))
                else{
                    ## times <- seq(start,maxtime,(maxtime - start)/100)
                    times <- median(eventTimes)
                }
            } else{
                if (use.event.times==TRUE) 
                    times <- sort(c(start,unique(times),eventTimes))
                else
                    ## times <- sort(unique(c(start,times)))
                    times <- sort(unique(times))
            }
            (if (any(times>maxtime))
                 message(paste0("Upper limit of followup is ",
                                maxtime,"\nResults at evaluation time(s) beyond this time point are not computed.")))
            ## need to save indices to modify matrix input
            include.times <- times<=maxtime
            times <- times[include.times]
            NT <-  length(times)
            if (NT==0)
                stop("No evaluation time before end of followup.")
        }
        else{
            stop("Landmark updating not yet implemented.")
        }
    } else{
        if (!missing(times) && (!is.null(times)) && (times!=FALSE)) warning("Function 'Score': Response type is not time-to-event: argument 'times' will be ignored.",call.=FALSE)
        times <- NULL
        NT <- 1
    }
    # }}}

    # -----------------Dealing with censored data outside the loop -----------------------
    # {{{
    if (responseType %in% c("survival","competing.risks")){
        if (censType=="rightCensored"){
            if ("outside.ipcw" %in% cens.method){
                Weights <- getCensoringWeights(formula=formula,
                                               data=data,
                                               response=data[,1:response.dim,with=FALSE],
                                               times=times,
                                               cens.model=cens.model,
                                               responseType=responseType)
            }
        } else{ if (censType=="uncensored"){
                    cens.method <- c("none","inside")
                    Weights <- NULL
                }
                else{
                    stop("Cannot handle this type of censoring.")
                }
        }
        if ("pseudo" %in% cens.method){
            if (censType=="rightCensored"
                && (responseType %in% c("survival","competing.risks"))){        
                margForm <- update(formula,paste(".~1"))
                margFit <- prodlim::prodlim(margForm,data=data)
                jack <- data.table(ID=data[["ID"]],
                                   times=rep(times,rep(N,NT)),
                                   pseudovalue=c(prodlim::jackknife(margFit,cause=cause,times=times)))
                if (responseType=="survival") jack[,pseudovalue:=1-pseudovalue]
                                    
            }
        }
    }
    # }}}
    # -----------------performance program----------------------
    # {{{ 
    trainModel <- function(model,data){
        model$call$data <- data
        try(eval(model$call),silent=FALSE)
    }
    computePerformance <- function(object,
                                   nullobject,
                                   testdata,
                                   testweights,
                                   traindata=NULL,
                                   trainseed=NULL,
                                   metrics,
                                   plots,
                                   times,
                                   cause,
                                   se.fit,
                                   keep.residuals,
                                   alpha,
                                   probs,
                                   dolist,
                                   DT.residuals=NULL,
                                   item,
                                   NF,
                                   NT){
        Brier=Rsquared=NULL
        N <- NROW(testdata)
        # split into response and predictors
        response <- testdata[,1:response.dim,with=FALSE]
        response[,ID:=testdata[["ID"]]]
        X <- testdata[,-c(1:response.dim),with=FALSE]
        # {{{ -----------------IPCW inner loop-----------------------
        if (responseType %in% c("survival","competing.risks")){
            if (censType=="rightCensored"){
                if ("inside.ipcw" %in% cens.method){
                    Weights <- getCensoringWeights(formula=formula,
                                                   data=testdata,
                                                   response=response,
                                                   times=times,
                                                   cens.model=cens.model,
                                                   responseType=responseType)
                } else{
                    if ("outside.ipcw" %in% cens.method){
                        Weights <- testweights
                    }else{
                        stop("don't know how to refit pseudo values in crossvalidation")
                        cens.method <- "jackknife.pseudo.values"
                    }
                }
                ## for both inside and outside IPCW we
                ## add subject specific weights here
                ## and time specific weights later
                response[,WTi:=Weights$IPCW.subject.times]
            } else {
                if (censType=="uncensored"){
                    Weights <- list(IPCW.times=rep(1,NT),
                                    IPCW.subject.times=matrix(1,ncol=NT,nrow=N))
                    Weights$method <- "marginal"
                    response[,WTi:=1]
                } else{
                    stop("Cannot handle this type of censoring.")
                }
            }
        }
        if (responseType=="binary") Weights <- NULL
        # }}}
        # extract predictions as data.table
        args <- switch(responseType,"binary"={list(newdata=X)},
                       "survival"={list(newdata=X,times=times)},
                       "competing.risks"={list(newdata=X,times=times,cause=cause)},
                       stop("Unknown responseType."))
        pred <- data.table::rbindlist(lapply(mlevs, function(f){
            if (f>0 && (length(extra.args <- unlist(lapply(object.classes[[f]],function(cc){predictRisk.args[[cc]]})))>0)){
                args <- c(args,extra.args)
            }
            if (f!=0 && any(c("integer","factor","numeric","matrix") %in% object.classes[[f]])){
                ## sort predictions by ID
                if (!is.null(dim(object[[f]]))) {## input matrix
                    if(!is.null(include.times)){ ## remove columns at times beyond max time
                        p <- c(do.call("predictRisk",c(list(object=object[[f]][,include.times,drop=FALSE]),args))[testdata[["ID"]],])
                    } else{ 
                        p <- c(do.call("predictRisk",c(list(object=object[[f]]),args))[testdata[["ID"]],])
                    }
                }
                else{ ## either binary or only one time point
                    p <- do.call("predictRisk", c(list(object=object[[f]]),args))[testdata[["ID"]]]
                }
            }else{
                if (!is.null(traindata)){
                    set.seed(trainseed)
                    ## remove response from traindata to avoid clash when model uses Hist(time,status)
                    ## where status has 0,1,2 but now event history response has status=0,1
                    nResp <- switch(responseType,"binary"=1,"survival"=2,"competing.risks"=3)
                    if (f==0)
                        trained.model <- trainModel(model=nullobject[[1]],data=traindata[,-c(1:nResp),with=FALSE])
                    else
                        trained.model <- trainModel(model=object[[f]],data=traindata[,-c(1:nResp),with=FALSE])
                    if ("try-error" %in% class(trained.model)){
                        stop(paste0("Failed to fit model ",f,ifelse(try(b>0,silent=TRUE),paste0(" in cross-validation step ",b,"."))))
                    }
                }
                else{
                    if (f==0)
                        trained.model <- nullobject[[1]]
                    else
                        trained.model <- object[[f]]
                }
                p <- c(do.call("predictRisk", c(list(object=trained.model),args)))
                if (f==0 && responseType!="binary") {## glm predicts the same value for all subjects
                    p <- rep(p,rep(N,NT))
                }
                ## ## predict risks not survival
                ## if (predictHandlerFun=="predictRisk") p <- 1-p
            }
            if (!is.null(times)){
                data.table(ID=testdata[["ID"]],model=f,risk=p,times=rep(times,rep(N,NT)))
            } else {
                data.table(ID=testdata[["ID"]],model=f,risk=p)
            }
        }))
        out <- vector(mode="list",length=length(c(summary,metrics,plots)))
        names(out) <- c(summary,metrics,plots)
        if (item ==0 & "Calibration" %in% plots){
            if (responseType=="binary" || censType=="uncensored")
                out[["Calibration"]]$plotframe <- merge(response,pred[model!=0],by="ID")
            else{
                out[["Calibration"]]$plotframe <- merge(jack,pred[model!=0],by=c("ID","times"))
            }
            out[["Calibration"]]$plotframe[,model:=factor(model,levels=mlevs,mlabels)]
        }
        if (!is.null(Weights)){
            if (Weights$method=="marginal"){
                Wt <- data.table(times=times,Wt=Weights$IPCW.times)
                ## FIXME: many digits in times may cause merge problems
                ## Wt$times <- factor(Wt$times)
                ## pred$times <- factor(pred$times)
                pred <- merge(pred,Wt,by=c("times"))
            }else{
                Wt <- data.table(times=rep(times,rep(N,NT)),
                                 Wt=c(Weights$IPCW.times),
                                 ID=testdata$ID)
                pred <- merge(pred,Wt,by=c("ID","times"))
            }
        }
        # compute and test performance
        ## input <- list("prediction"=pred,response=response)
        input <- list(DT=merge(response,pred,by="ID"),
                      N=N,
                      NT=NT,
                      NF=NF,
                      alpha=alpha,
                      se.fit=se.fit,
                      keep.residuals=keep.residuals,
                      dolist=dolist,Q=probs,ROC=FALSE)
        if (responseType=="competing.risks")
            input <- c(input,list(cause=cause,states=states))
        if (responseType %in% c("survival","competing.risks") && se.fit>0L){
            if (censType=="rightCensored"){
                input <- c(input,list(MC=response[,getInfluenceCurve.KM(time=time,status=status)]))
            }
            else{
                input <- c(input,list(MC=matrix(0,ncol=length(response$time),nrow=length(unique(response$time)))))
            }
        }
        ## make sure that brier comes first, so that we can remove the null.model afterwards
        for (m in sort(metrics,decreasing=TRUE)){
            if (m=="AUC" && ("ROC" %in% plots)){
                input <- replace(input, "ROC",TRUE)
                ## call AUC method
                out[[m]] <- do.call(paste(m,responseType,sep="."),input)
                out[["ROC"]]$plotframe <- out[[m]]$ROC
                out[["ROC"]]$plotframe[,model:=factor(model,levels=mlevs,mlabels)]
                out[[m]]$ROC <- NULL
            }else{                
                input <- replace(input, "ROC",FALSE)
                ## call Brier or AUC method
                out[[m]] <- do.call(paste(m,responseType,sep="."),input)
            }
            out[[m]]$score[,model:=factor(model,levels=mlevs,mlabels)]
            ## set model and reference in model comparison results
            if (NROW(out[[m]]$contrasts)>0){
                out[[m]]$contrasts[,model:=factor(model,levels=mlevs,mlabels)]
                out[[m]]$contrasts[,reference:=factor(reference,levels=mlevs,mlabels)]
            }
        }
        ## summary should be after metrics because R^2 depends on Brier score
        if (rsquared){
            if (responseType=="binary")
                out[["Brier"]][["score"]][,Rsquared:=1-Brier/Brier[model=="Null model"]]
            else
                out[["Brier"]][["score"]][,Rsquared:=1-Brier/Brier[model=="Null model"],by=times]
        }
        
        for (s in summary){
            if (s=="risks") {
                out[[s]] <- list(score=copy(input$DT),contrasts=NULL)
            } else{
                out[[s]] <- do.call(paste(s,responseType,sep="."),input)
            }
            out[[s]]$score[,model:=factor(model,levels=mlevs,mlabels)]
            if (NROW(out[[s]]$contrasts)>0){
                out[[s]]$contrasts[,model:=factor(model,levels=mlevs,mlabels)]
                out[[s]]$contrasts[,reference:=factor(reference,levels=mlevs,mlabels)]
            }
        }
        
        out
    }
    # }}}
    # -----------------apparent nosplit performance---------------------
    # {{{
    noSplit <- computePerformance(object=object,
                                  nullobject=nullobject,
                                  testdata=data,
                                  testweights=Weights,
                                  metrics=metrics,
                                  plots=plots,
                                  times=times,
                                  cause=cause,
                                  se.fit=se.fit,
                                  keep.residuals=keep.residuals,
                                  alpha=alpha,
                                  probs=probs,
                                  dolist=dolist,
                                  DT.residuals=NULL,
                                  item=0,
                                  ## DT.bootcount=NULL
                                  NF=NF,
                                  NT=NT)
    # }}}
    # -----------------crossvalidation performance---------------------
    # {{{ 
    crossval <- NULL
    if (split.method$name%in%c("BootCv","LeaveOneOutBoot")){
        if (missing(trainseeds)||is.null(trainseeds)){
            if (!missing(seed)) set.seed(seed)
            trainseeds <- sample(1:1000000,size=B,replace=FALSE)
        }
        crossval <- foreach (b=1:B) %dopar%{
            traindata=data[split.method$index[,b],,drop=FALSE]
            ## subset.data.table preserves order
            if (split.method$name=="BootCv"){
                DT.residuals=NULL
                ## DT.bootcount=NULL
            }else{
                DT.residuals=data.table(residuals=rep(0,N))
                DT.influence=data.table(influence=rep(0,N*N))
                ## DT.bootcount=NULL
            }
            testids <- (match(1:N,unique(split.method$index[,b]),nomatch=0)==0)
            testdata <- subset(data,testids,drop=FALSE)
            if ((censType=="rightCensored")
                &&
                (responseType %in% c("survival","competing.risks"))
                &&
                ("outside.ipcw" %in% cens.method)){
                testweights <- Weights
                ## browser(skipCalls=1L)
                testweights$IPCW.subject.times <- subset(testweights$IPCW.subject.times,testids,drop=FALSE)
                if (Weights$dim>0){
                    testweights$IPCW.times <- subset(testweights$IPCW.times,testids,drop=FALSE)
                }
            } else {
                testweights <- NULL
            }
            ## case se.fit=we 1 need only p-values for multi-split tests
            se.fit.cv <- se.fit*2
            cb <- computePerformance(object=object,
                                     nullobject=nullobject,
                                     testdata=testdata,
                                     testweights=testweights,
                                     traindata=traindata,
                                     trainseed=trainseeds[b],
                                     metrics=metrics,
                                     plots=plots,
                                     times=times,
                                     cause=cause,
                                     se.fit=se.fit.cv,
                                     keep.residuals=keep.residuals,
                                     alpha=alpha,
                                     probs=probs,
                                     dolist=dolist,
                                     item=b,
                                     NF=NF,
                                     NT=NT)
            cb
        }
        ##  # ------ Leave-one-out bootstrap cens.method
        if (split.method$name=="LeaveOneOutBoot"){  
            crossvalPerf <- lapply(metrics,function(m){
                if (m=="Brier"){
                    ## if (test==TRUE){crossval[[1]][[m]]$score)>0
                    ## Brier <- data.table::rbindlist(lapply(crossval,function(x)x[["Brier"]]))
                    ## Brier[,list(looboot=mean(ipcwResiduals)),by=list(model,times,ID)]
                    ## bootcv=Brier[,list(bootcv=mean(ipcwResiduals)),by=list(model,times,b)]
                    if (length(crossval[[1]][[m]]$score)>0){
                        Ibi=perf=.SD=ipcwResiduals=NULL
                        if (null.model==TRUE) init <- data.table(ID=rep(1:N, each=(NF+1)*B), b=rep(1:B, each=(NF+1), times=N), model=rep(0:NF, times=B))
                        else init <- data.table(ID=rep(1:N, each=NF*B), b=rep(1:B, each=NF, times=N), model=rep(1:NF, times=B))
                        
                        looBoot <- data.table::rbindlist(lapply(seq_along(crossval), function(x) crossval[[x]][[m]]$residual[,b:=x]))
                        looBoot <- merge(looBoot, init, by=c("ID", "b", "model"), all=T)
                        looBoot[,Ibi:=1*(!is.na(ipcwResiduals))]
                        data.table::set(looBoot, i=which(is.na(looBoot[["ipcwResiduals"]])), j="ipcwResiduals", value=0)
                        looBoot[,perf:=do.call("/", c(lapply(.SD, sum))), by=c("ID", "model"), .SDcols=c("ipcwResiduals", "Ibi")]
                        data.table::set(looBoot, j=c("b", "ipcwResiduals","Ibi"), value=NULL)
                        looBoot <- unique(looBoot)
                        looBoot <- looBoot[, sum(perf)/N, by="model"]
                        out <- list(looBoot)
                    }
                }
            })
        }else{ ## split.method bootcv
            crossvalPerf <- lapply(metrics,function(m){
                if (responseType %in% c("survival","competing.risks")){
                    byvars <- c("model","times")
                } else{
                    byvars <- c("model")
                }
                ## score
                if (length(crossval[[1]][[m]]$score)>0){
                    cv.score <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$score}))
                    bootcv.score <- cv.score[,data.table::data.table(mean(.SD[[m]],na.rm=TRUE),
                                                                     lower=quantile(.SD[[m]],alpha/2,na.rm=TRUE),
                                                                     upper=quantile(.SD[[m]],(1-alpha/2),na.rm=TRUE)),by=byvars,.SDcols=m]
                    data.table::setnames(bootcv.score,c(byvars,m,"lower","upper"))
                }else{
                    cv.score <- NULL
                    bootcv.score <- NULL
                }
                ## contrasts and multi-split test
                if (length(crossval[[1]][[m]]$contrasts)>0){
                    cv.contrasts <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$contrasts}))
                    delta.m <- paste0("delta.",m)
                    bootcv.contrasts <- cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE),
                                                                             lower=quantile(.SD[[delta.m]],alpha/2,na.rm=TRUE),
                                                                             upper=quantile(.SD[[delta.m]],(1-alpha/2),na.rm=TRUE),
                                                                             p=median(p,na.rm=TRUE)),by=byvars,.SDcols=delta.m]
                    data.table::setnames(bootcv.contrasts,c(byvars,delta.m,"lower","upper","p"))
                }else{ 
                    cv.contrasts <- NULL
                    bootcv.contrasts <- NULL
                }
                out <- list(score=bootcv.score,contrasts=bootcv.contrasts)
                if (keep.cv)
                    out <- c(out,list(cv.score=cv.score,cv.contrasts=cv.contrasts))
                out
            })
            names(crossvalPerf) <- metrics
        }
    }
    # }}}
    #------------------output-----------------------------------
    # {{{ enrich the output object
    if (is.null(crossval))
        output <- noSplit
    else{
        ## output <- list(noSplitPerf=noSplit,crossValPerf=crossvalPerf)
        output <- crossvalPerf
    }
    models <- mlevs
    names(models) <- mlabels
    if (null.model==TRUE) nm <- names(models)[1] else nm <- NULL
    if (!keep.splitindex) split.method$index <- "split index was not saved"
    output <- c(output,list(responseType=responseType,
                            dolist=dolist,
                            cause=cause,
                            states=states,
                            null.model=nm,
                            models=models,
                            censType=censType,
                            censoringHandling=cens.method,
                            split.method=split.method,
                            metrics=metrics,
                            plots=plots,
                            summary=summary,
                            call=theCall))
    
    for (p in c(plots)){
        output[[p]]$plotmethod <- p
        class(output[[p]]) <- paste0("score",p)
    }
    for (m in c(metrics,summary)){
        output[[m]]$metrics <- m
        class(output[[m]]) <- paste0("score",m)
    }
    class(output) <- "Score"
    output
# }}}
}

##' @export
Score <- function(object,...){
    UseMethod("Score",object=object)
}

# {{{ Brier score

Brier.binary <- function(DT,se.fit,alpha,N,NT,NF,dolist,keep.residuals=FALSE,DT.residuals,DT.bootcount,...){
    residuals=Brier=risk=model=ReSpOnSe=lower=upper=se=ID=NULL
    DT[,residuals:=(ReSpOnSe-risk)^2,by=model]
    if (se.fit>0L){
        ## data.table::setorder(DT,model,ReSpOnSe)
        data.table::setkey(DT,model,ID)
        score <- DT[,data.table::data.table(Brier=sum(residuals)/N,se=sd(residuals)/sqrt(N)),by=list(model)]
        if (se.fit==1L){   
            score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        }else{
            score[,se:=NULL]
        } 
        data.table::setkey(score,model)
        data.table::setkey(DT,model)
        DT <- DT[score]
        if (length(dolist)>0){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,IF=residuals,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist,se.fit=se.fit)]
            setnames(contrasts.Brier,"delta","delta.Brier")
            output <- list(score=score,contrasts=contrasts.Brier)
        }else{
            output <- list(score=score)
        }
    }else{
        output <- list(score=DT[,list(Brier=mean(residuals)),by=list(model)])
    }
    if (keep.residuals) output <- c(output,list(residuals=DT[,data.table::data.table(ID,ReSpOnSe,model,risk,residuals)]))
    output
}

Brier.survival <- function(DT,MC,se.fit,alpha,N,NT,NF,dolist,keep.residuals=FALSE,DT.residuals,DT.bootcount,...){
    Yt=time=times=Residuals=risk=Brier=ipcwResiduals=WTi=Wt=status=setorder=model=IF.Brier=data.table=sd=lower=qnorm=se=upper=NULL
    ## compute 0/1 outcome:
    DT[,Yt:=1*(time<=times)]
    ## compute residuals
    DT[,Residuals:=(Yt-risk)^2]
    ## apply weights 
    DT[,ipcwResiduals:=Residuals/WTi]
    DT[Yt==0,ipcwResiduals:=Residuals/Wt]
    ## deal with censored observations before t
    DT[Yt==1 & status==0,ipcwResiduals:=0]
    ## DT[,c("Yt","time","WTi","Wt","status","Residuals"):=NULL]
    if (se.fit>0L){
        data.table::setorder(DT,model,times,time,-status)
        DT[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],time=time,Yt=Yt,ipcwResiduals=ipcwResiduals,MC=MC),by=list(model,times)]
        score <- DT[,data.table(Brier=sum(ipcwResiduals)/N,se=sd(IF.Brier)/sqrt(N)),by=list(model,times)]
        if (se.fit==1L){   
            score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        }else{
            score[,se:=NULL]
        } 
        data.table::setkey(score,model,times)
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        if (length(dolist)>0L){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist,se.fit=se.fit),by=list(times)]
            setnames(contrasts.Brier,"delta","delta.Brier")
            output <- list(score=score,contrasts=contrasts.Brier)
        }
        else{
            output <- list(score=score)
        }
    }else{
        output <- list(score=DT[,data.table(Brier=sum(ipcwResiduals)/N),by=list(model,times)])
    }
    if (keep.residuals) output <- c(output,list(residuals=DT))
    output
}

Brier.competing.risks <- function(DT,MC,se.fit,alpha,N,NT,NF,dolist,keep.residuals=FALSE,DT.residuals,DT.bootcount,cause,states,...){
    Yt=time=times=event=Brier=Residuals=risk=ipcwResiduals=WTi=Wt=status=setorder=model=IF.Brier=data.table=sd=lower=qnorm=se=upper=NULL
    ## compute 0/1 outcome:
    DT[,Yt:=1*(time<=times & event==cause)]
    ## compute residuals
    DT[,Residuals:=(Yt-risk)^2]
    ## apply weights 
    DT[,ipcwResiduals:=Residuals/WTi]
    DT[time>times,ipcwResiduals:=Residuals/Wt]
    ## deal with censored observations before t
    DT[time<=times & status==0,ipcwResiduals:=0]
    ## DT[,c("Yt","time","WTi","Wt","status","Residuals"):=NULL]
    if (se.fit>0L){
        data.table::setorder(DT,model,times,time,-status)
        DT[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],time=time,Yt=Yt,ipcwResiduals=ipcwResiduals,MC=MC),by=list(model,times)]
        score <- DT[,data.table(Brier=sum(ipcwResiduals)/N,se=sd(IF.Brier)/sqrt(N)),by=list(model,times)]
        if (se.fit==1L){   
            score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        }else{
            score[,se:=NULL]
        } 
        data.table::setkey(score,model,times)
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        if (length(dolist)>0){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist,se.fit=se.fit),by=list(times)]
            setnames(contrasts.Brier,"delta","delta.Brier")
            output <- list(score=score,contrasts=contrasts.Brier)
        } else{
            output <- list(score=score)
        }
    }else{
        output <- list(score=DT[,data.table(Brier=sum(ipcwResiduals)/N),by=list(model,times)])
    }
    if (keep.residuals) output <- c(output,list(residuals=DT))
    output
}

# }}}

# {{{ AUC

## helper functions 
AireTrap <- function(FP,TP,N){
    N <- length(FP)
    sum((FP-c(0,FP[-N]))*((c(0,TP[-N])+TP)/2))
}

## do not want to depend on Daim as they turn marker to ensure auc > 0.5
delongtest <-  function(risk, score, dolist, response, cause, alpha, se.fit) {
    cov=lower=upper=p=AUC=se=lower=upper=NULL
    Cases <- response == cause
    Controls <- response != cause
    nControls <- sum(!Cases)
    nCases <- sum(Cases)
    ## if (nbCases==0 || nbControls ==0 || length(unique(risk))==1) return(rep(0,n))
    nauc <- ncol(risk)
    auc <- score[["AUC"]]
    modelnames <- score[["model"]]
    riskcontrols <- as.matrix(risk[Controls,])
    riskcases <- as.matrix(risk[Cases,])
    V10 <- matrix(0, nrow = nCases, ncol = nauc)
    V01 <- matrix(0, nrow = nControls, ncol = nauc)
    tmn <- t(riskcontrols)
    tmp <- t(riskcases)
    for (i in 1:nCases) {
        V10[i, ] <- rowSums(tmn < tmp[, i]) + 0.5 * rowSums(tmn == tmp[, i])
    }
    for (i in 1:nControls) {
        V01[i, ] <- rowSums(tmp > tmn[, i]) + 0.5 * rowSums(tmp == tmn[, i])
    }
    V10 <- V10/nControls
    V01 <- V01/nCases
    W10 <- cov(V10)
    W01 <- cov(V01)
    S <- W10/nCases + W01/nControls
    se.auc <- sqrt(diag(S))
    score <- data.table(model=colnames(risk),AUC=auc)
    if (se.fit==1L){
        score[,se:=se.auc]
        score[,lower:=pmax(0,AUC-qnorm(1-alpha/2)*se)]
        score[,upper:=pmin(1,AUC+qnorm(1-alpha/2)*se)]
        setcolorder(score,c("model","AUC","se","lower","upper"))
    }else{ 
        setcolorder(score,c("model","AUC"))
    }
    names(auc) <- 1:nauc
    ## q1 <- auc/(2 - auc)
    ## q2 <- 2 * auc^2/(1 + auc)
    ## aucvar <- (auc * (1 - auc) + (nCases - 1) * (q1 - auc^2) + (nControls - 1) * (q2 - auc^2))/(nCases * nControls)
    if (length(dolist)>0){
        ncomp <- nauc * (nauc - 1)/2
        delta.AUC <- numeric(ncomp) 
        se <- numeric(ncomp)
        model <- numeric(ncomp)
        reference <- numeric(ncomp)
        ctr <- 1
        Qnorm <- qnorm(1 - alpha/2)
        for (i in 1:(nauc - 1)) {
            for (j in (i + 1):nauc) {
                delta.AUC[ctr] <- auc[j]-auc[i]
                ## cor.auc[ctr] <- S[i, j]/sqrt(S[i, i] * S[j, j])
                LSL <- t(c(1, -1)) %*% S[c(j, i), c(j, i)] %*% c(1, -1)
                ## print(c(1/LSL,rms::matinv(LSL)))
                se[ctr] <- sqrt(LSL)
                ## tmpz <- (delta.AUC[ctr]) %*% rms::matinv(LSL) %*% delta.AUC[ctr]
                ## tmpz <- (delta.AUC[ctr]) %*% (1/LSL) %*% delta.AUC[ctr]
                model[ctr] <- modelnames[j]
                reference[ctr] <- modelnames[i]
                ctr <- ctr + 1
            }
        }
        deltaAUC <- data.table(model,reference,delta.AUC=as.vector(delta.AUC),se)
        deltaAUC[,p:=2*pnorm(abs(delta.AUC)/se,lower.tail=FALSE)]
        if (se.fit==1L){ 
            deltaAUC[,lower:=delta.AUC-Qnorm*se]
            deltaAUC[,upper:=delta.AUC+Qnorm*se]
        }
        list(score = score, contrasts = deltaAUC)
    }else{
        list(score = score, contrasts = NULL)
    }
}

auRoc.numeric <- function(X,D,breaks,ROC){
    if (is.null(breaks)) breaks <- rev(sort(unique(X))) ## need to reverse when high X is concordant with {response=1}
    TPR <- c(prodlim::sindex(jump.times=X[D==1],eval.times=breaks,comp="greater",strict=FALSE)/sum(D==1))
    FPR <- c(prodlim::sindex(jump.times=X[D==0],eval.times=breaks,comp="greater",strict=FALSE)/sum(D==0))
    if (ROC==TRUE)
        data.table(risk=breaks,TPR,FPR)
    else
        0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
}
auRoc.factor <- function(X,D,ROC){
    TPR <- (sum(D==1)-table(X[D==1]))/sum(D==1)
    FPR <- table(X[D==0])/sum(D==0)
    if (ROC==TRUE)
        data.table(cbind(risk=c(sort(unique(X))),TPR,FPR))
    else
        0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
}

AUC.binary <- function(DT,breaks=NULL,se.fit,alpha,N,NT,NF,dolist,ROC,...){
    model=risk=ReSpOnSe=FPR=TPR=ID=NULL
    aucDT <- DT[model>0]
    dolist <- dolist[sapply(dolist,function(do){match("0",do,nomatch=0L)})==0]
    data.table::setkey(aucDT,model,ID)
    if (is.factor(DT[["risk"]])){
        score <- aucDT[,auRoc.factor(risk,ReSpOnSe,ROC=ROC),by=list(model)]
    }
    else{
        score <- aucDT[,auRoc.numeric(risk,ReSpOnSe,breaks=NULL,ROC=ROC),by=list(model)]
    }
    if (ROC==FALSE){
        setnames(score,"V1","AUC")
        output <- list(score=score)
    } else{
        AUC <- score[,list(AUC=0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))),by=list(model)]
        ROC <- score
        output <- list(score=AUC,ROC=ROC)
    }
    if (se.fit>0L){
        xRisk <- data.table::dcast.data.table(aucDT[],ID~model,value.var="risk")[,-1,with=FALSE]
        delong.res <- delongtest(risk=xRisk,
                                 score=output$score,
                                 dolist=dolist,
                                 response=aucDT[model==model[1],ReSpOnSe],
                                 cause="1",
                                 se.fit=se.fit,
                                 alpha=alpha)
        output$score <- delong.res$score
        output$contrasts <- delong.res$contrasts
        output
    }else{
        output
    }
}


AUC.survival <- function(DT,MC,se.fit,alpha,N,NT,NF,dolist,ROC,...){
    model=times=risk=Cases=time=status=Controls=TPR=FPR=WTi=Wt=ipcwControls=ipcwCases=IF.AUC=lower=se=upper=AUC=NULL
    cause <- 1
    aucDT <- DT[model>0]
    dolist <- dolist[sapply(dolist,function(do){match("0",do,nomatch=0L)})==0]
    ## assign Weights before ordering
    aucDT[,ipcwControls:=1/(Wt*N)]
    aucDT[,ipcwCases:=1/(WTi*N)]
    ## order data
    data.table::setorder(aucDT,model,times,-risk)
    ## identify cases and controls
    aucDT[,Cases:=(time <= times &  status==cause)]
    aucDT[,Controls:=(time > times)] 
    ## prepare Weights
    aucDT[Cases==0,ipcwCases:=0]
    aucDT[Controls==0,ipcwControls:=0]
    ## compute denominator
    aucDT[,TPR:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    aucDT[,FPR:=(cumsum(ipcwControls))/(sum(ipcwControls)),by=list(model,times)]
    nodups <- aucDT[,c(!duplicated(risk)[-1],TRUE),by=list(model,times)]$V1
    if (ROC==TRUE) {
        output <- list(ROC=aucDT[nodups,c("model","times","risk","TPR","FPR"),with=FALSE])
    }else{
        output <- NULL
    }
    score <- aucDT[nodups,list(AUC=AireTrap(FPR,TPR)),by=list(model,times)]
    data.table::setkey(score,model,times)
    if (se.fit>0L){
        ## compute influence function
        data.table::setorder(aucDT,model,times,time,-status)
        aucDT[,IF.AUC:=getInfluenceCurve.AUC.survival(t=times[1],n=N,time=time,risk=risk,Cases=Cases,Controls=Controls,ipcwControls=ipcwControls,ipcwCases=ipcwCases,MC=MC), by=list(model,times)]
        se.score <- aucDT[,list(se=sd(IF.AUC)/sqrt(N)),by=list(model,times)]
        data.table::setkey(se.score,model,times)
        score <- score[se.score]
        if (se.fit==1L){   
            score[,lower:=pmax(0,AUC-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,AUC+qnorm(1-alpha/2)*se)]
        }else{
            score[,se:=NULL]
        }
        data.table::setkey(aucDT,model,times)
        aucDT <- aucDT[score]
        if (length(dolist)>0){
            contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),
                                                   NF=NF,
                                                   N=N,
                                                   alpha=alpha,
                                                   dolist=dolist,se.fit=se.fit),by=list(times)]
            setnames(contrasts.AUC,"delta","delta.AUC")
            output <- c(list(score=score,contrasts=contrasts.AUC),output)
        }else{
            output <- c(list(score=score),output)
        }
    }else{
        output <- c(list(score=score),output)
    }
    output
}


AUC.competing.risks <- function(DT,MC,se.fit,alpha,N,NT,NF,dolist,cause,states,ROC,...){
    model=times=risk=Cases=time=status=event=Controls1=Controls2=TPR=FPR=WTi=Wt=ipcwControls1=ipcwControls2=ipcwCases=IF.AUC=lower=se=upper=AUC=NULL
    aucDT <- DT[model>0]
    dolist <- dolist[sapply(dolist,function(do){match("0",do,nomatch=0L)})==0]
    ## assign Weights before ordering
    aucDT[,ipcwControls1:=1/(Wt*N)]
    aucDT[,ipcwControls2:=1/(WTi*N)]
    aucDT[,ipcwCases:=1/(WTi*N)]
    aucDT[,ipcwControls2:=1/(WTi*N)]
    ## order data
    data.table::setorder(aucDT,model,times,-risk)
    ## identify cases and controls
    aucDT[,Cases:=(time <=times &  event==cause)]
    aucDT[,Controls1:=(time > times)] 
    aucDT[,Controls2:=(time <=times &  event!=cause & status !=0)]
    ## prepare Weights
    aucDT[Cases==0,ipcwCases:=0]
    aucDT[Controls1==0,ipcwControls1:=0]
    aucDT[Controls2==0,ipcwControls2:=0]
    ## compute denominator
    ## ROC <- aucDT[,list(TPR=c(0,cumsum(ipcwCases)),FPR=c(0,cumsum(ipcwControls1)+cumsum(ipcwControls2))),by=list(model,times)]
    aucDT[,TPR:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    aucDT[,FPR:=(cumsum(ipcwControls1)+cumsum(ipcwControls2))/(sum(ipcwControls2)+sum(ipcwControls1)),by=list(model,times)]
    nodups <- aucDT[,c(!duplicated(risk)[-1],TRUE),by=list(model,times)]$V1
    if (ROC==TRUE) {
        output <- list(ROC=aucDT[nodups,c("model","times","risk","TPR","FPR"),with=FALSE])
    }else{
        output <- NULL
    }
    score <- aucDT[nodups,list(AUC=AireTrap(FPR,TPR)),by=list(model,times)]
    data.table::setkey(score,model,times)
    if (se.fit>0L){
        ## compute influence function
        data.table::setorder(aucDT,model,times,time,-status)
        aucDT[,IF.AUC:=getInfluenceCurve.AUC.competing.risks(t=times[1],n=N,time=time,risk=risk,ipcwControls1=ipcwControls1,ipcwControls2=ipcwControls2,ipcwCases=ipcwCases,Cases=Cases,Controls1=Controls1,Controls2=Controls2,MC=MC), by=list(model,times)]
        se.score <- aucDT[,list(se=sd(IF.AUC)/sqrt(N)),by=list(model,times)]
        data.table::setkey(se.score,model,times)
        score <- score[se.score]
        if (se.fit==1L){
            score[,lower:=pmax(0,AUC-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,AUC+qnorm(1-alpha/2)*se)]
        }else{
            score[,se:=NULL]
        }
        data.table::setkey(aucDT,model,times)
        aucDT <- aucDT[score]
        if (length(dolist)>0){
            contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist,se.fit=se.fit),by=list(times)]
            setnames(contrasts.AUC,"delta","delta.AUC")
            output <- c(list(score=score,contrasts=contrasts.AUC),output)
        }else{
            output <- c(list(score=score),output)
        }
    }else{
        output <- c(list(score=score),output)
    }
    output
}

# }}}


