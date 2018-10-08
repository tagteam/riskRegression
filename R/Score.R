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
##' For uncensored binary outcome the Delong-Delong test is used to contrast AUC of rival models.
##' In right censored survival data (with and without competing risks) 
##' the p-values correspond to Wald tests based on standard errors obtained with an estimate of the influence function
##' as described in detail in the appendix of Blanche et al. (2015).
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
##' @param plots Character vector specifying for which plots to put data into the result.
##'     Currently implemented are \code{"ROC"}, \code{"Calibration"} and \code{"boxplot"}.
##'     In addition, one can plot AUC and Brier score as function of time as soon as
##'     \code{times} has at least two different values.
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
##' @param se.fit Logical or \code{0} or \code{1}. If \code{FALSE} or \code{0} do not calculate standard errors.
##' @param conservative Logical, only relevant in right censored data. If \code{TRUE} ignore
##'        variability of the estimate of the inverse probability of censoring weights when calculating standard
##'        errors for prediction performance parameters. This can potentially reduce computation time and memory usage
##'        at a usually very small expense of a slightly higher standard error. 
##' @param multi.split.test Logical or \code{0} or \code{1}. If \code{FALSE} or \code{0} do not calculate multi-split tests. This argument is ignored when \code{split.method} is "none".
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
##'     data. Either \code{"ipcw"} or \code{"pseudo"}.
##'     Here IPCW refers to inverse probability of censoring weights and \code{pseudo} for jackknife pseudo values.
##'     Right now pseudo values  are only used for calibration curves.
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
##' @param parallel The type of parallel operation to be used (if any). If missing, the default is \code{"no"}.
##' @param ncpus integer: number of processes to be used in parallel operation.
##' @param cl An optional \code{parallel} or \code{snow} cluster for use if \code{parallel = "snow"}. If not supplied, a cluster on the local machine is created for the duration of the \code{Score} call.
##' @param keep list of characters (not case sensitive) which determines additional output.
##' \code{"residuals"} provides Brier score residuals and
##' \code{"splitindex"} provides sampling index used to split the data into training and validation sets.
##' \code{"vcov"} provides the variance-covariance matrix of the estimated parameters. 
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
##' @param debug Logical. If \code{TRUE} indicate landmark in progress of the program.
##' @param useEventTimes obsolete. 
##' @param nullModel obsolete.
##' @param censMethod obsolete.
##' @param censModel obsolete.
##' @param splitMethod obsolete.
##' @param ... Named list containging additional arguments that are passed on to the \code{predictRisk} methods corresponding to object. See examples.
##' @return List with scores and assessments of contrasts, i.e.,
##'     tests and confidence limits for performance and difference in performance (AUC and Brier),
##'     summaries and plots. Most elements are in\code{data.table} format.
##' @details
##' This function works with one or multiple models that predict the risk of an event R(t|X) for a subject
##' characterized by predictors X at time t.
##' With binary endpoints (outcome 0/1 without time component) the risk is simply R(X).
##' In case of a survival object
##' without competing risks the function still works with predicted event probabilities, i.e., R(t|X)=1-S(t|X) where
##' S(t|X) is the predicted survival chance for subject X at time t.
##'
##' The already existing predictRisk methods (see methods(predictRisk)) may not cover all models and methods
##' for predicting risks. But users can quickly extend the package as explained in detail in Mogensen et al. (2012) for
##' the predecessors \code{pec::predictSurvProb} and \code{pec::predictEventProb} which have been unified as
##' \code{riskRegression::predictRisk}.
##'
##' Bootstrap Crossvalidation (see also Gerds & Schumacher 2007 and Mogensen et al. 2012)
##'
##' B=10, M (not specified or M=NROW(data))
##' Training of each of the models in each of 10 bootstrap data sets (learning data sets).
##' Learning data sets are obtained by sampling \code{NROW(data)} subjects of the data set
##' with replacement. There are roughly \code{.632*NROW(data)} subjects in the learning data (inbag)
##' and \code{.368*NROW(data)} subjects not in the validation data sets (out-of-bag).
##' 
##' These are used to estimate the scores: AUC, Brier, etc. Reported are averages across the 10 splits.
##'
##' ## Bootstrap with replacement
##' \code{
##' set.seed(13)
##' N=17
##' data = data.frame(id=1:N, y=rbinom(N,1,.3),x=rnorm(N))
##' boot.index = sample(1:N,size=N,replace=TRUE)
##' boot.index
##' inbag = 1:N %in% boot.index
##' outofbag = !inbag 
##' learn.data = data[inbag]
##' val.data = data[outofbag]
##' riskRegression:::getSplitMethod("bootcv",B=10,N=17)$index
##' }
##' NOTE: the number .632 is the expected probability to draw one subject (for example subject 1) with
##'    replacement from the data, which does not depend on the sample size:
##' \code{B=10000}
##' \code{N=137}
##' \code{mean(sapply(1:B, function(b){match(1,sample(1:N,size=N,replace=TRUE),nomatch=0)}))}
##' \code{N=30}
##' \code{mean(sapply(1:B, function(b){match(1,sample(1:N,size=N,replace=TRUE),nomatch=0)}))}
##' \code{N=300}
##' \code{mean(sapply(1:B, function(b){match(1,sample(1:N,size=N,replace=TRUE),nomatch=0)}))}
##'
##' 
##' ## Bootstrap without replacement (training size set to be 70 percent of data)
##' B=10, M=.7
##' 
##' Training of each of the models in each of 10 bootstrap data sets (learning data sets).
##' Learning data sets are obtained by sampling \code{round(.8*NROW(data))} subjects of the data set
##' without replacement. There are \code{NROW(data)-round(.8*NROW(data))} subjects not in the learning data sets.
##' These are used to estimate the scores: AUC, Brier, etc. Reported are averages across the 10 splits.
##' \code{
##' set.seed(13)
##' N=17
##' data = data.frame(id=1:N, y=rbinom(N,1,.3),x=rnorm(N))
##' boot.index = sample(1:N,size=M,replace=FALSE)
##' boot.index
##' inbag = 1:N %in% boot.index
##' outofbag = !inbag 
##' learn.data = data[inbag]
##' val.data = data[outofbag]
##' riskRegression:::getSplitMethod("bootcv",B=10,N=17,M=.7)$index
##' }
##' 
##' @examples
##' # binary outcome
##' library(lava)
##' set.seed(18)
##' learndat <- sampleData(48,outcome="binary")
##' testdat <- sampleData(40,outcome="binary")
##'
##' ## score logistic regression models
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
##' lr2 = glm(Y~X3+X5,data=learndat,family=binomial)
##' Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5)"=lr2),formula=Y~1,data=testdat)
##'
##' ## ROC curve and calibration plot 
##' xb=Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,
##'         data=testdat,plots=c("calibration","ROC"))
##' plotROC(xb)
##' plotCalibration(xb)
##'
##' ## compute AUC for a list of continuous markers
##' markers = as.list(testdat[,.(X6,X7,X8,X9,X10)])
##' Score(markers,formula=Y~1,data=testdat,metrics=c("auc"))
##'
##' # cross-validation
##' \dontrun{
##' learndat=sampleData(400,outcome="binary")
##' lr1a = glm(Y~X6,data=learndat,family=binomial)
##' lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
##' ## bootstrap cross-validation
##' x1=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=100)
##' x1
##' ## leave-one-out and leave-pair-out bootstrap
##' x2=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,
##'               split.method="loob",
##'               B=100,plots="calibration")
##' x2
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
##' library(survival)
##' set.seed(18)
##' trainSurv <- sampleData(400,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
##' x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
##'       split.method="loob",B=100,plots="calibration")
##' 
##' x2= Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
##'       split.method="bootcv",B=100)
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
##'\dontrun{
##' # reproduce some results of Table IV of Blanche et al. Stat Med 2013
##' data(Paquid)
##' ResPaquid <- Score(list("DSST"=-Paquid$DSST,"MMSE"=-Paquid$MMSE),
##'                    formula=Hist(time,status)~1,
##'                    data=Paquid,
##'                    null.model = FALSE,
##'                    conf.int=TRUE,
##'                    metrics=c("auc"),
##'                    times=c(3,5,10),
##'                    plots="ROC")
##' ResPaquid
##' plotROC(ResPaquid,time=5)
##' }
##' @author Thomas A Gerds \email{tag@@biostat.ku.dk} and Paul Blanche \email{paul.blanche@@univ-ubs.fr}
##' @references
##'
##' Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
##' Evaluating Random Forests for Survival Analysis Using Prediction Error
##' Curves. Journal of Statistical Software, 50(11), 1-23. URL
##' http://www.jstatsoft.org/v50/i11/.
##'
##' Paul Blanche, Cecile Proust-Lima, Lucie Loubere, Claudine Berr, Jean- Francois Dartigues, and
##' Helene Jacqmin-Gadda. Quantifying and comparing dynamic predictive accuracy of joint models
##' for longitudinal marker and time-to-event in presence of censoring and competing risks.
##' Biometrics, 71 (1):102--113, 2015.
##' 
##' P. Blanche, J-F Dartigues, and H. Jacqmin-Gadda. Estimating and comparing
##' time-dependent areas under receiver operating characteristic curves for
##' censored event times with competing risks. Statistics in Medicine,
##' 32(30):5381--5397, 2013.
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
# {{{ Score.list
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
                       se.fit=TRUE,
                       conservative=FALSE,
                       multi.split.test=FALSE,
                       conf.int=.95,
                       contrasts=TRUE,
                       probs=c(0,0.25,0.5,0.75,1),
                       cens.method="ipcw",
                       cens.model="cox",
                       split.method,
                       B,
                       M,
                       seed,
                       trainseeds,
                       parallel=c("no","multicore","snow"),
                       ncpus=1,
                       cl=NULL,
                       keep,
                       predictRisk.args,
                       debug=0L,
                       useEventTimes,
                       nullModel,
                       censMethod,
                       censModel,
                       splitMethod,
                       ...){
    if (!missing(useEventTimes)) {
        warning("(Spelling of) argument 'useEventTimes' is obsolete and will be disabled in future versions of this function.\nUse 'use.event.times' instead.")
        use.event.times <- useEventTimes
    }
    if (!missing(nullModel)) {
        warning("(Spelling of) argument 'nullModel' is obsolete and will be disabled in future versions of this function.\nUse 'null.model' instead.")
        null.model <- nullModel
    }
    if (!missing(splitMethod)) {
        warning("(Spelling of) argument 'splitMethod' is obsolete and will be disabled in future versions of this function.\nUse 'split.method' instead.")
        split.method <- splitMethod
    }
    if (!missing(censMethod)) {
        warning("(Spelling of) argument 'censMethod' is obsolete and will be disabled in future versions of this function.\nUse 'cens.method' instead.")
        cens.method <- censMethod
    }
    if (!missing(censModel)) {
        warning("(Spelling of) argument 'censModel' is obsolete and will be disabled in future versions of this function.\nUse 'cens.model' instead.")
        cens.model <- censModel
    }

    se.conservative=IF.AUC.conservative=IF.AUC0=IF.AUC=IC0=Brier=AUC=casecontrol=se=nth.times=time=status=ID=WTi=ID=risk=IF.Brier=lower=upper=crossval=b=time=status=model=reference=p=model=pseudovalue=ReSpOnSe=residuals=event=j=NULL

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
    if (length(posRR <- grep("^ipa$|^rr$|^r2|rsq",summary,ignore.case=TRUE))>0){
        if (!null.model) stop("Need the null model to compute IPA/R^2 but argument 'null.model' is FALSE.")
        summary <- summary[-posRR]
        if (!("Brier" %in% metrics)) metrics <- c(metrics,"Brier")
        if (!null.model) {
            null.model <- TRUE 
            warning("Value of argument 'null.model' ignored as the null model is needed to compute IPA/R^2.")
        }
        ipa <- TRUE
    }else{
        ipa <- FALSE
    }
    ## Plots <- lapply(plots,grep,c("Roc","Cal"),ignore.case=TRUE,value=TRUE)
    if ("ROC" %in% plots) {
        ## add AUC if needed
        if (!("AUC" %in% metrics)) metrics <- c(metrics,"AUC")
    }
    if ("Calibration" %in% plots) {
        ## add pseudo if needed
        if (!("pseudo" %in% cens.method)) cens.method <- c(cens.method,"pseudo")
    }

    # }}}
    # -----------------parse other arguments and prepare data---------
    # {{{ censoring model arguments
    if (length(grep("^km|^kaplan|^marg",cens.model,ignore.case=TRUE))>0){
        cens.model <- "KaplanMeier"
    } else{
        if (length(attr(terms(formula),"factors"))>0){
            cens.model <- "cox"
        } else{
            cens.model <- "KaplanMeier"
        }
    }
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
    response.type <- attr(response,"model")
    if (response.type %in% c("survival","competing.risks")){
        byvars <- c("model","times")
    } else{
        byvars <- c("model")
    }
    states <- attr(response,"states")
    if (missing(cause)||is.null(cause)){
        if (response.type=="binary"){
            cause <- attr(response,"event")
        }
        else{
            cause <- states[[1]]
        }
    }
    position.cause <- match(cause,states,nomatch=0)
    if (position.cause==0) stop(paste0("Requested cause: ",cause,". Available causes: ", paste(states,collapse=",")))
    # add null model and find names for the object
    if (null.model==TRUE){
        nullobject <- getNullModel(formula=formula,data=data,response.type=response.type)
    } else{
        nullobject <- NULL
    }
    ## put ReSpOnSe for binary and (time, event, status) in the first column(s) 
    ## data[,eval(responsevars):=NULL]
    data <- cbind(response,data)
    N <- NROW(data)
    # data have to be ordered when ipcw is called
    if (response.type %in% c("survival","competing.risks")){    
        neworder <- data[,order(time,-status)]
        data.table::setorder(data,time,-status)
    }else{
        neworder <- 1:N
    }
    # add ID variable for merging purposes and because output has long format
    data[,ID:=1:N]
    if (response.type=="binary")
        data[,event:=factor(ReSpOnSe,levels=1:(length(states)+1),labels=c(states,"censored"))]
    if (response.type=="survival")
        formula <- stats::update(formula,"Hist(time,status)~.")
    if (response.type=="competing.risks")
        formula <- stats::update(formula,"Hist(time,event)~.")
    ## stop("Dont know how to predict response of type ",response.type))
    cens.type <- attr(response,"cens.type")
    if (is.null(cens.type)) cens.type <- "uncensoredData"
    rm(response)
    # }}}
    # {{{ SplitMethod & parallel stuff
    if (!missing(seed)) set.seed(seed)
    split.method <- getSplitMethod(split.method=split.method,B=B,N=N,M=M)
    B <- split.method$B
    splitIndex <- split.method$index
    do.resample <- !(is.null(splitIndex))
    if (split.method$internal.name!="noplan"){
        if (missing(parallel)) parallel <- "no"
        parallel <- match.arg(parallel)
        if (ncpus<=1) {
            parallel <- "no"
        } 
        ## if (ncpus <- pmin(ncpus,parallel::detectCores()))
        switch(parallel,"no"={
            foreach::registerDoSEQ()
        },"multicore"={
            if(.Platform$OS.type == "windows"){stop("multicore does not work on windows. set argument parallel to 'snow'.")}
            loadNamespace("parallel")
            if (is.null(cl)) cl <- parallel::makeForkCluster(nnodes=ncpus)
            doParallel::registerDoParallel(cl=cl,cores=ncpus)
            on.exit(parallel::stopCluster(cl))
        },"snow"={
            loadNamespace("parallel")
            if (is.null(cl)) cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
            on.exit(parallel::stopCluster(cl))
            doParallel::registerDoParallel(cl=cl,cores=ncpus)                
        })
        if (split.method$name=="BootCv" && multi.split.test==TRUE){
            if ("AUC" %in% metrics) {
                warning("Cannot do multi-split test with AUC yet. Forced multi.split.test=FALSE")
                multi.split.test=FALSE
            }else{
                warning("Cannot do deal with conservative=FALSE when also multi.split.test=TRUE. Forced conservative=TRUE.")
                conservative=TRUE
            }
        }
        if (se.fit==TRUE){
            if (response.type=="competing.risks")
                warning("Under construction. Check devtools::install_github('tagteam/riskRegression') for progress.")
        }
        if ("ROC" %in% plots){
            warning("Cannot (not yet) do ROC analysis in combination with internal validation\n. Check devtools::install_github('tagteam/riskRegression') for progress.")
        }
    }
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
            if (split.method$internal.name!="noplan")
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
        if("vcov" %in% tolower(keep)) keep.vcov=TRUE else keep.vcov = FALSE
        if ("splitindex" %in% tolower(keep)) keep.splitindex=TRUE else keep.splitindex = FALSE
        if ("cv" %in% tolower(keep)) keep.cv=TRUE else keep.cv = FALSE
    }else{
        keep.residuals=FALSE
        keep.vcov=FALSE
        keep.cv=FALSE
        keep.splitindex=FALSE
    }
    # }}}
    # {{{ resolve se.fit and contrasts
    if (missing(se.fit)){
        if (is.logical(conf.int) && conf.int==FALSE
            || conf.int<=0
            || conf.int>1) 
            se.fit <- 0L
        else
            se.fit <- 1L
    }else{
        stopifnot(is.logical(se.fit)||se.fit==0||se.fit==1)
    }
    if (split.method$internal.name=="noplan") multi.split.test <- FALSE
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
        if ((is.logical(contrasts) && contrasts==FALSE)){
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
    # {{{ Evaluation landmarks and horizons (times)
    if (response.type %in% c("survival","competing.risks")){
        ## in case of a tie, events are earlier than right censored
        eventTimes <- unique(data[,time])
        maxtime <- eventTimes[length(eventTimes)]
        if (debug==TRUE){
            cat("\nThe maxtime is set at:",maxtime,"\n")
        }
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
    if (response.type %in% c("survival","competing.risks")){
        if (cens.type=="rightCensored"){
            if (se.fit>0L && "AUC" %in% metrics && conservative==TRUE) {
                ## FIXME: need conservative formula for AUC 
                warning("Cannot do conservative==TRUE with AUC yet.")
            }
            if (se.fit>0L && "AUC" %in% metrics && cens.model=="cox"){
                if (!(split.method$name %in% c("LeaveOneOutBoot","BootCv"))){
                    warning("Cannot (not yet) estimate standard errors for AUC with Cox IPCW.\nTherefore, force cens.model to be marginal.")
                    cens.model <- "KaplanMeier"
                }
            }
            Weights <- getCensoringWeights(formula=formula,
                                           data=data,
                                           times=times,
                                           cens.model=cens.model,
                                           response.type=response.type,
                                           ## FIXME: need conservative formula for AUC
                                           influence.curve=(se.fit==TRUE && (conservative==0L || "AUC" %in% metrics)))
            ## 
            ## if cens.model is marginal then IC is a matrix (ntimes,newdata) 
            ## if cens.model is Cox then IC is an array (nlearn, ntimes, newdata)
            ## IC is an array with dimension (nlearn, times, newdata)
            ## IC_G(t,z;x_k)
        }else { 
            if (cens.type=="uncensored"){
                cens.method <- c("none")
                Weights <- NULL
            }
            else{
                stop("Cannot handle this type of censoring.")
            }
        }
        if ("pseudo" %in% cens.method){
            if (cens.type=="rightCensored"
                && (response.type %in% c("survival","competing.risks"))){        
                margForm <- update(formula,paste(".~1"))
                margFit <- prodlim::prodlim(margForm,data=data)
                ## position.cause is the result of match(cause, states)
                jack <- data.table(ID=data[["ID"]],
                                   times=rep(times,rep(N,NT)),
                                   pseudovalue=c(prodlim::jackknife(margFit,cause=position.cause,times=times)))
                if (response.type=="survival") jack[,pseudovalue:=1-pseudovalue]
            }
        }
    } else{ 
        Weights <- NULL
    }
    # }}}
    # -----------------performance program----------------------
    # {{{ define getPerformanceData, setup data in long-format, response, pred, weights, model, times, loop
    # {{{ header
    getPerformanceData <- function(testdata,
                                   testweights,
                                   traindata=NULL,
                                   trainseed=NULL){

        ID=NULL
        # inherit everything else from parent frame: object, nullobject, NF, NT, times, cause, response.type, etc.
        Brier=IPA=NULL
        looping <- !is.null(traindata)
        ## if (!looping) b=0
        N <- NROW(testdata)
        # split data vertically into response and predictors X
        response <- testdata[,1:response.dim,with=FALSE]
        response[,ID:=testdata[["ID"]]]
        setkey(response,ID)
        X <- testdata[,-c(1:response.dim),with=FALSE]
        if (debug) if (looping) message(paste0("Loop round: ",b)) 
        if (debug) message("extracted test set and prepared output object")
        # }}}
# {{{ collect pred as long format data.table

        args <- switch(response.type,"binary"={list(newdata=X)},
                       "survival"={list(newdata=X,times=times)},
                       "competing.risks"={list(newdata=X,times=times,cause=cause)},
                       stop("Unknown response.type."))
        ## remove our cbinded response (see above) from traindata to avoid clash when model uses Hist(time,status)
        ## where status has 0,1,2 but now event history response has status=0,1
        ## the original response is still there
        if(!is.null(traindata)){
            trainX <- traindata[,-c(1:response.dim),with=FALSE]
            trainX[,ID:=NULL]
        }
        pred <- data.table::rbindlist(lapply(mlevs, function(f){
            if (f>0 && (length(extra.args <- unlist(lapply(object.classes[[f]],function(cc){predictRisk.args[[cc]]})))>0)){
                args <- c(args,extra.args)
            }
            ## predictions given as numeric values
            if (f!=0 && any(c("integer","factor","numeric","matrix") %in% object.classes[[f]])){
                ## sort predictions by ID
                if (!is.null(dim(object[[f]]))) {## input matrix
                    if(!is.null(include.times)){ ## remove columns at times beyond max time
                        p <- c(do.call("predictRisk",c(list(object=object[[f]][,include.times,drop=FALSE]),args))[neworder,])
                    } else{ 
                        p <- c(do.call("predictRisk",c(list(object=object[[f]]),args))[neworder,])
                    }
                }
                else{ ## either binary or only one time point
                    p <- do.call("predictRisk", c(list(object=object[[f]]),args))[neworder]
                }
            }else{
                # predictions given as model which needs training in crossvalidation loops
                if (looping){
                    set.seed(trainseed)
                    if (f==0) model.f=nullobject[[1]] else model.f=object[[f]]
                    model.f$call$data <- trainX
                    trained.model <- try(eval(model.f$call),silent=TRUE)
                    if ("try-error" %in% class(trained.model)){
                        message(paste0("Failed to fit model ",f,ifelse(looping,paste0(" in cross-validation step ",b)),":"))
                        try(eval(model.f$call),silent=FALSE)
                        stop()
                    }
                } else{
                    if (f==0)
                        trained.model <- nullobject[[1]]
                    else
                        trained.model <- object[[f]]
                }
                p <- c(do.call("predictRisk", c(list(object=trained.model),args)))
                if (f==0 && response.type!="binary") {## glm predicts the same value for all subjects
                    p <- rep(p,rep(N,NT))
                }
            }
            if (response.type%in%c("survival","competing.risks")){
                out <- data.table(ID=testdata[["ID"]],model=f,risk=p,times=rep(times,rep(N,NT)))
                setkey(out,model,times,ID)
                out
            } else {
                out <- data.table(ID=testdata[["ID"]],model=f,risk=p)
                setkey(out,model,ID)
                out
            }
        }))
        if (debug) message("trained the model(s) and extracted the predictions")
        
        # }}}
# {{{ merge with Weights (IPCW inner loop)

        if (response.type %in% c("survival","competing.risks")){
            if (cens.type=="rightCensored"){
                Weights <- testweights
                ## add subject specific weights 
                response[,WTi:=Weights$IPCW.subject.times]
            } else {
                if (cens.type=="uncensored"){
                    Weights <- list(IPCW.times=rep(1,NT),IPCW.subject.times=matrix(1,ncol=NT,nrow=N))
                    Weights$method <- "marginal"
                    response[,WTi:=1]
                } else{
                    stop("Cannot handle this type of censoring.")
                }
            }
            ## add time point specific weights 
            if (Weights$method=="marginal"){
                Wt <- data.table(times=times,Wt=Weights$IPCW.times)
                ## OBS: many digits in times may cause merge problems
                pred <- merge(pred,Wt,by=c("times"))
            }else{
                Wt <- data.table(times=rep(times,rep(N,NT)),
                                 Wt=c(Weights$IPCW.times),
                                 ID=testdata$ID)
                pred <- merge(pred,Wt,by=c("ID","times"))
            }
            if (debug) message("merged the weights with input for performance metrics")
        } else {
            ## if (response.type=="binary")
            Weights <- NULL
        }
        if (debug) message("added weights to predictions")

        # }}}
        # {{{ merge with response
        DT=merge(response,pred,by="ID")
        DT
    }
    # }}}
    # }}} 
    # {{{ define function to test performance
    computePerformance <- function(DT,
                                   N,
                                   se.fit,
                                   conservative,
                                   cens.model,
                                   multi.split.test,
                                   keep.residuals,
                                   keep.vcov){
        IPA=Brier=NULL
        # inherit everything else from parent frame: summary, metrics, plots, alpha, probs, dolist, et
        out <- vector(mode="list",
                      length=length(c(summary,metrics,plots)))
        names(out) <- c(summary,metrics,plots)
        input <- list(DT=DT,
                      N=N,
                      NT=NT,
                      NF=NF,
                      alpha=alpha,
                      se.fit=se.fit,
                      conservative=conservative,
                      cens.model=cens.model,
                      multi.split.test=multi.split.test,
                      keep.residuals=keep.residuals,
                      keep.vcov=keep.vcov,
                      ## DT.residuals=DT.residuals,
                      dolist=dolist,Q=probs,ROC=FALSE,MC=Weights$IC)
        # {{{ collect data for calibration plots
        if ("Calibration" %in% plots){
            if (response.type=="binary" || cens.type=="uncensored")
                out[["Calibration"]]$plotframe <- DT[model!=0]
            else{
                out[["Calibration"]]$plotframe <- merge(jack,DT[model!=0],by=c("ID","times"))
            }
            out[["Calibration"]]$plotframe[,model:=factor(model,levels=mlevs,mlabels)]
        }
        # }}}
        if (response.type=="competing.risks") {
            input <- c(input,list(cause=cause,states=states))
        }
        ## make sure that Brier score comes first, so that we can remove the null.model afterwards
        for (m in sort(metrics,decreasing=TRUE)){
            if (m=="AUC" && ("ROC" %in% plots)){
                input <- replace(input, "ROC",TRUE)
                ## call AUC method
                out[[m]] <- do.call(paste(m,response.type,sep="."),input)
                out[["ROC"]]$plotframe <- out[[m]]$ROC
                out[["ROC"]]$plotframe[,model:=factor(model,levels=mlevs,mlabels)]
                out[[m]]$ROC <- NULL
            }else{                
                input <- replace(input, "ROC",FALSE)
                ## call Brier or AUC method
                out[[m]] <- do.call(paste(m,response.type,sep="."),input)
            }
            if (!is.null(out[[m]]$score)){
                out[[m]]$score[,model:=factor(model,levels=mlevs,mlabels)]
            }
            ## set model and reference in model comparison results
            if (!is.null(out[[m]]$contrasts)>0){
                out[[m]]$contrasts[,model:=factor(model,levels=mlevs,mlabels)]
                out[[m]]$contrasts[,reference:=factor(reference,levels=mlevs,mlabels)]
            }
        }
        ## summary should be after metrics because IPA/R^2 depends on Brier score
        if (ipa){
            if (response.type=="binary")
                out[["Brier"]][["score"]][,IPA:=1-Brier/Brier[model=="Null model"]]
            else
                out[["Brier"]][["score"]][,IPA:=1-Brier/Brier[model=="Null model"],by=times]
        }
        for (s in summary){
            if (s=="risks") {
                out[[s]] <- list(score=copy(input$DT),contrasts=NULL)
            } else{
                out[[s]] <- do.call(paste(s,response.type,sep="."),input)
            }
            if (NROW(out[[s]]$contrasts)>0){
                out[[s]]$score[,model:=factor(model,levels=mlevs,mlabels)]
            }
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

    missing.predictions <- "Don't know yet"
    if (split.method$internal.name %in% c("noplan",".632+")){
        DT <- getPerformanceData(testdata=data,
                                 testweights=Weights,
                                 traindata=NULL,
                                 trainseed=NULL)
        if (any(is.na(DT[["risk"]]))){
            missing.predictions <- DT[,list("Missing.values"=sum(is.na(risk))),by=byvars]
            warning("Missing values in the predicted risk. See `missing.predictions' in output list.")
        }else{
            missing.predictions <- "None"
        }
        noSplit <- computePerformance(DT,
                                      N=N,
                                      se.fit=se.fit,
                                      conservative=conservative,
                                      cens.model,
                                      multi.split.test=multi.split.test,
                                      keep.residuals=keep.residuals,
                                      keep.vcov=keep.vcov)
        if (debug) message("computed apparent performance")
    }

    # }}}
    # -----------------crossvalidation performance---------------------
    # {{{
    # {{{ bootstrap re-fitting
    if (split.method$name%in%c("BootCv","LeaveOneOutBoot")){
        if (missing(trainseeds)||is.null(trainseeds)){
            if (!missing(seed)) set.seed(seed)
            trainseeds <- sample(1:1000000,size=B,replace=FALSE)
        }
        if (parallel=="snow") exports <- c("data","split.method","Weights","N","trainseeds") else exports <- NULL
        DT.B <- rbindlist(foreach (b=1:B,.export=exports) %dopar%{
            ## DT.B <- rbindlist(lapply(1:B,function(b){
            traindata=data[split.method$index[,b]]
            ## setkey(traindata,ID)
            testids <- (match(1:N,unique(split.method$index[,b]),nomatch=0)==0)
            ## NOTE: subset.data.table preserves order
            testdata <- subset(data,testids)
            if (cens.type=="rightCensored"){
                testweights <- Weights
                testweights$IPCW.subject.times <- subset(testweights$IPCW.subject.times,testids)
                if (Weights$dim>0){
                    testweights$IPCW.times <- subset(testweights$IPCW.times,testids)
                }
            } else {
                testweights <- NULL
            }
            DT.b <- getPerformanceData(testdata=testdata,
                                       testweights=testweights,
                                       traindata=traindata,
                                       trainseed=trainseeds[[b]])
            DT.b[,b:=b]
            DT.b
        })
        if (any(is.na(DT.B[["risk"]]))){
            missing.predictions <- DT.B[,list("Missing.values"=sum(is.na(risk))),by=byvars]
            warning("Missing values in the predicted risk. See `missing.predictions' in output list.")
        }
        ## FIXME: subset influence curves
        ## case se.fit=1 we need only p-values for multi-split tests
        ## se.fit.cv <- se.fit*2
        ## cb <- computePerformance(DT.b,
        ## se.fit=se.fit.cv,
        ## multi.split.test=multi.split.test,
        ## keep.residuals=keep.residuals)
        ## cb
        if (debug) message("setup data for cross-validation performance")
        Response <- data[,c(1:response.dim),with=FALSE]
        Response[,ID:=data[["ID"]]]
        setkey(Response,ID)
        # }}}
        # {{{ Leave-one-out bootstrap
        if (split.method$name=="LeaveOneOutBoot"){
            crossvalPerf <- lapply(metrics,function(m){
                # {{{ AUC LOOB
                if (m=="AUC"){
                    if (response.type=="binary")
                        auc.loob <- data.table(model=mlevs)
                    else
                        auc.loob <- data.table(expand.grid(times=times,model=mlevs))
                    auc.loob[,AUC:=as.numeric(NA)]
                    ## for each pair of individuals sum the concordance of the bootstraps where *both* individuals are out-of-bag 
                    ## divide by number of times the pair is out-of-bag later
                    if (se.fit==TRUE){
                        aucDT <- NULL
                    }
                    if (response.type=="binary"){
                        NT <- 1
                    }
                    for (s in 1:NT){
                        t <- times[s]
                        if (response.type=="binary"){
                            ## the following indices have to be logical!!!
                            cases.index <- Response[,ReSpOnSe==1]
                            controls.index <- !cases.index
                            cc.status <- factor(cases.index,levels=c(TRUE,FALSE),labels=c("case","control"))
                        }else{
                            if (response.type=="survival"){
                                ## event of interest before times
                                ## the following indices have to be logical!!!
                                cases.index <- Response[,time<=t & status==1]
                                controls.index <- Response[,time>t]
                                cc.status <- factor(rep("censored",N),levels=c("censored","case","control"))
                                cc.status[cases.index] <- "case"
                                cc.status[controls.index] <- "control"
                            }
                            else{ ## competing risks
                                ## event of interest before times
                                ## the following indices have to be logical!!!
                                cases.index <- Response[,time<=t & status==1]
                                controls.index <- Response[,time>t | status==2]
                                cc.status <- factor(rep("censored",N),levels=c("censored","case","control"))
                                cc.status[cases.index] <- "case"
                                cc.status[controls.index] <- "control"
                            }
                        }
                        if (cens.type=="rightCensored"){
                            ## IPCW
                            weights.cases <- cases.index/Weights$IPCW.subject.times
                            if (Weights$method=="marginal"){
                                weights.controls <- controls.index/Weights$IPCW.times[s]
                            }else{
                                weights.controls <- controls.index/Weights$IPCW.times[,s]
                            }
                            weightMatrix <- outer(weights.cases[cases.index], weights.controls[controls.index], "*")
                        }else{ ## uncensored
                            weights.cases <- cases.index/1
                            weights.controls <- controls.index/1
                            weightMatrix <- outer(weights.cases[cases.index], weights.controls[controls.index], "*")
                        }
                        Phi <- (1/N^2)*sum(weights.cases[cases.index])*sum(weights.controls[controls.index])
                        which.cases <- (1:N)[cases.index]
                        which.controls <- (1:N)[controls.index]
                        for (mod in mlevs){
                            Ib <- matrix(0, sum(cases.index), sum(controls.index))
                            auc <- matrix(0, sum(cases.index), sum(controls.index))
                            for (u in 1:B){## cannot use b as running index because b==b does not work in data.table
                                ## test <- DT.B[model==mod&times==t&b==u]
                                oob <- match(1:N,unique(split.method$index[,u]),nomatch=0)==0
                                ## to use the cpp function AUCijFun we
                                ## need a vector of length equal to the number of cases (which.cases) for the current time point
                                ## which has arbitrary values in places where subjects are inbag and the predicted risks
                                ## for those out-of-bag. need another vector for controls.
                                riskset <- data.table::data.table(ID=1:N,casecontrol=cc.status,oob=oob)
                                data.table::setkey(riskset,ID)
                                if (is.null(t)){
                                    oob.risk <- DT.B[model==mod&b==u,data.table::data.table(ID,risk)]
                                }else{
                                    oob.risk <- DT.B[model==mod&times==t&b==u,data.table::data.table(ID,risk)]
                                }
                                data.table::setkey(oob.risk,ID)
                                riskset <- oob.risk[riskset]
                                riskset[is.na(risk),risk:=-9]
                                Ib.ij <- outer((cases.index*oob)[which.cases],(controls.index*oob)[which.controls],"*")
                                auc.ij <- AUCijFun(riskset[casecontrol=="case",risk],riskset[casecontrol=="control",risk])*Ib.ij
                                ## Ib.ij is 1 when the pair out of bag
                                ## print(head(oob))
                                ## print(auc.ij[1:5,1:5])
                                auc <- auc+auc.ij
                                Ib <- Ib + Ib.ij
                            }
                            auc <- (auc*weightMatrix)/Ib
                            # FIXME: why are there NA's?
                            auc[is.na(auc)] <- 0
                            ## Leave-one-pair-out bootstrap estimate of AUC
                            aucLPO <- (1/N^2)*sum(colSums(auc))*(1/Phi)
                            if (is.null(t)){
                                auc.loob[model==mod,AUC:=aucLPO]
                            }else{
                                auc.loob[times==t&model==mod,AUC:=aucLPO]
                            }
                            if (se.fit==1L){
                                ## ## First part of influence function
                                ic0Case <- rowSums(auc)
                                ic0Control <- colSums(auc)
                                ic0 <- (1/(Phi*N))*c(ic0Case, ic0Control)-2*aucLPO
                                id.cases <- data[["ID"]][cc.status=="case"]
                                id.controls <- data[["ID"]][cc.status=="control"]
                                id.censored <- data[["ID"]][cc.status=="censored"]
                                if (is.null(t)){
                                    this.aucDT <- data.table(model=mod,ID = c(id.cases,id.controls), IF.AUC0 = ic0)
                                }else{
                                    this.aucDT <- data.table(model=mod,times=t,ID = c(id.cases,id.controls,id.censored), IF.AUC0 = c(ic0, rep(-2*aucLPO,length(id.censored))))
                                }
                                aucDT <- rbindlist(list(aucDT,this.aucDT),use.names=TRUE,fill=TRUE)
                                if (response.type=="binary" || cens.type=="uncensored"){
                                    icPhi <- (aucLPO/Phi)*((weights.cases-(1/N)*sum(weights.cases))*(1/N)*sum(weights.controls)+(weights.controls-(1/N)*sum(weights.controls))*(1/N)*sum(weights.controls))-2*aucLPO
                                    if (is.null(t)){
                                        data.table::setkey(aucDT,model,ID)
                                        aucDT[model==mod,IF.AUC:=IF.AUC0-icPhi]
                                        auc.loob[model==mod,se:=sd(aucDT[["IF.AUC"]])/sqrt(N)]
                                    }else{
                                        data.table::setkey(aucDT,model,times,ID)
                                        aucDT[times==t&model==mod,IF.AUC:=IF.AUC0-icPhi]
                                        auc.loob[times==t&model==mod,se:=sd(aucDT[["IF.AUC"]])/sqrt(N)]
                                    }
                                }else{
                                    if (cens.type=="rightCensored" && (conservative==FALSE)) {
                                        ## ## Influence function for G - i.e. censoring survival distribution 
                                        ic.weights <- matrix(0,N,N)
                                        if  (cens.model=="cox"){
                                            k=0 ## counts subject-times with event before t
                                            for (i in 1:N){
                                                if (i %in% id.cases){
                                                    k=k+1
                                                    ic.weights[i,] <- Weights$IC$IC.subject[i,k,]/(Weights$IPCW.subject.times[i])
                                                }else{
                                                    if (i %in% id.controls){ ## min(T,C)>t
                                                        ic.weights[i,] <- Weights$IC$IC.times[i,s,]/(Weights$IPCW.times[i,s])
                                                    }
                                                }
                                            }
                                        }else{
                                            k=0 ## counts subject-times with event before t
                                            for (i in 1:N){
                                                if (i %in% id.cases){
                                                    ## FIXME: need to check IC
                                                    pos.i <- sindex(jump.times=unique(data[["time"]]),eval.times=data[["time"]][i])
                                                    ic.weights[i,] <- Weights$IC[pos.i,]/(Weights$IPCW.subject.times[i])
                                                }else{
                                                    if (i %in% id.controls){ ## min(T,C)>t
                                                        ic.weights[i,] <- Weights$IC[s,]/(Weights$IPCW.times[s])
                                                    }
                                                }
                                            }
                                        }
                                        ## ## Part of influence function related to Weights
                                        ic.weightsCase <- as.numeric(rowSumsCrossprod(as.matrix(rowSums(auc)), ic.weights[which.cases,], 0))  
                                        ic.weightsControl <- as.numeric(rowSumsCrossprod(as.matrix(colSums(auc)), ic.weights[which.controls,], 0)) 
                                        ic.weightsCC <- (1/(Phi*N^2))*(ic.weightsCase+ic.weightsControl)
                                    }
                                    ## ## Part of influence function related to Phi
                                    ## icPhiCase <- colMeans(ic.weights[which.cases,])
                                    icPhiCase <- as.numeric(rowSumsCrossprod(as.matrix(weights.cases[which.cases]),ic.weights[which.cases,],0))
                                    icPhiControl <- as.numeric(rowSumsCrossprod(as.matrix(weights.controls[which.controls]),ic.weights[which.controls,],0))
                                    icPhi <- (aucLPO/Phi)*((weights.cases-(1/N)*sum(weights.cases))*(1/N)*sum(weights.controls)+(weights.controls-(1/N)*sum(weights.controls))*(1/N)*sum(weights.cases)) - 2*aucLPO    
                                    ## ## Combine all parts of influence function
                                    ## ic1 <- data.table(ID=data[["ID"]], "ic.weightsCC" = ic.weightsCC, "icPhi" = icPhi)
                                    data.table::setkey(aucDT,model,times,ID)
                                    if(conservative==TRUE){
                                        aucDT[model==mod&times==t, IF.AUC:=IF.AUC0-icPhi]
                                    }else{
                                        aucDT[model==mod&times==t, IF.AUC:=IF.AUC0-ic.weightsCC-icPhi]
                                        aucDT[model==mod&times==t, IF.AUC.conservative:=IF.AUC0-icPhi]
                                    }
                                    aucDT[,IF.AUC0:=NULL]
                                    auc.loob[model==mod&times==t,se:= sd(aucDT[model==mod&times==t,IF.AUC])/sqrt(N)]
                                    auc.loob[model==mod&times==t,se.conservative:=sd(aucDT[model==mod&times==t,IF.AUC.conservative])/sqrt(N)]
                                    ## testSE <- sqrt(sum(ic[["ic"]]^2))/N
                                }
                            }
                        }
                    }
                    if (se.fit==1L){
                        auc.loob[,lower:=pmax(0,AUC-qnorm(1-alpha/2)*se)]
                        auc.loob[,upper:=pmin(1,AUC + qnorm(1-alpha/2)*se)]
                    }
                    output <- list(score=auc.loob)
                    if (length(dolist)>0L){
                        if (se.fit==FALSE){
                            if (match("times",byvars,nomatch=0)){
                                contrasts.AUC <- auc.loob[,getComparisons(data.table(x=AUC,model=model),
                                                                          NF=NF,
                                                                          N=N,
                                                                          alpha=alpha,
                                                                          dolist=dolist,
                                                                          se.fit=FALSE),by=list(times)]
                            } else{
                                contrasts.AUC <- auc.loob[,getComparisons(data.table(x=AUC,model=model),
                                                                          NF=NF,
                                                                          N=N,
                                                                          alpha=alpha,
                                                                          dolist=dolist,
                                                                          se.fit=FALSE)]
                            }
                        }else{
                            aucDT <- merge(aucDT,auc.loob,by=byvars)
                            if (match("times",byvars,nomatch=0)){
                                contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),
                                                                       NF=NF,
                                                                       N=N,
                                                                       alpha=alpha,
                                                                       dolist=dolist,
                                                                       se.fit=TRUE),by=list(times)]
                            } else{
                                contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),
                                                                       NF=NF,
                                                                       N=N,
                                                                       alpha=alpha,
                                                                       dolist=dolist,
                                                                       se.fit=TRUE)]
                            }
                        }
                        output <- c(output,list(contrasts=contrasts.AUC))
                    }
                    if (!is.null(output$score)){
                        output$score[,model:=factor(model,levels=mlevs,mlabels)]
                    }
                    ## set model and reference in model comparison results
                    if (!is.null(output$contrasts)>0){
                        output$contrasts[,model:=factor(model,levels=mlevs,mlabels)]
                        output$contrasts[,reference:=factor(reference,levels=mlevs,mlabels)]
                    }
                    return(output)
                }

                # }}}
                # {{{ Brier LOOB
                if (m=="Brier"){
                    ## sum across bootstrap samples where subject i is out of bag
                    if (cens.type=="rightCensored"){
                        if (response.type=="survival"){
                            ## event of interest before times
                            DT.B[time<=times & status==1,residuals:=(1-risk)^2/WTi]
                        }
                        else{ ## competing risks
                            ## event of interest before times
                            DT.B[time<=times & status==1 & event==cause,residuals:=(1-risk)^2/WTi]
                            ## competing event before times
                            DT.B[time<=times & status==1 &event!=cause,residuals:=(0-risk)^2/WTi]
                        }   
                        ## right censored before times
                        DT.B[time<=times & status==0,residuals:=0]
                        ## no event at times
                        DT.B[time>times,residuals:=(risk)^2/Wt]
                    }else{
                        DT.B[,residuals:=(ReSpOnSe-risk)^2]
                    }
                    ## for each individual sum the residuals of the bootstraps where this individual is out-of-bag 
                    ## divide by number of times out-off-bag later
                    ## DT.B <- DT.B[,data.table::data.table(residuals=sum(residuals)),by=c(byvars,"ID")]
                    DT.B <- DT.B[,data.table::data.table(risk=mean(risk),residuals=sum(residuals)),by=c(byvars,"ID")]
                    ## get denominator
                    Ib <- split.method$B-tabulate(unlist(apply(split.method$index,2,unique)))
                    ## REMOVE ME
                    ## Ib <- Ib[order(data$ID)]
                    if (any(Ib==0)) {
                        warning("Some subjects are never out of bag.\nResults are unreliable until You increase the number of bootstrap replications (argument 'B').")
                        Ib.include <- Ib!=0
                        Ib <- Ib[Ib.include]
                        ## don't subset residuals, they are only
                        ## available for those with Ib.include==1
                        ## DT.B <- DT.B[Ib.include]
                    } else Ib.include <- NULL
                    ## Ib is the count of how many times subject i is out of bag
                    ## the order of Ib matches the order of ID
                    ## within groups defined by times and model (byvars)
                    data.table::setkeyv(DT.B,c(byvars,"ID"))
                    DT.B[,residuals:=residuals/Ib]
                    ## leave-one-out bootstrap estimate 
                    DT.B[,Brier:=mean(residuals),by=byvars]
                    ## standard error via influence function
                    if (se.fit==1L){
                        ## influence function when censoring model is known or data are uncensored
                        DT.B[,IC0:=residuals-Brier]
                        ## se.brier <- DT.B[,list(se=sd(IC0, na.rm=TRUE)/sqrt(N)),by=byvars]
                        ## DT.B[,Brier:=NULL]
                        if (cens.type=="rightCensored" && conservative!=TRUE){
                            ## this is a new DT.B
                            DT.B <- cbind(data[,1:response.dim,with=FALSE],DT.B)
                            DT.B[,nth.times:=as.numeric(factor(times))]
                            DT.B[,WTi:=Weights$IPCW.subject.times]
                            if (Weights$method=="marginal"){
                                Wt <- data.table(times=times,Wt=Weights$IPCW.times)
                                ## OBS: many digits in times may cause merge problems
                                DT.B <- merge(DT.B,Wt,by=c("times"))
                            }else{
                                Wt <- data.table(times=rep(times,rep(N,NT)),
                                                 Wt=c(Weights$IPCW.times),
                                                 ID=data$ID)
                                DT.B <- merge(DT.B,Wt,by=c("ID","times"))
                            }
                            DT.B[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                                                    time=time,
                                                                    IC0,
                                                                    residuals=residuals,
                                                                    WTi=WTi,
                                                                    Wt=Wt,
                                                                    IC.G=Weights$IC,
                                                                    cens.model=cens.model,
                                                                    nth.times=nth.times[1]),by=byvars]
                            score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,
                                                           se=sd(IF.Brier)/sqrt(N),
                                                           se.conservative=sd(IC0)/sqrt(N)),by=byvars]
                        }else{
                            ## either conservative == TRUE or binary or uncensored
                            score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,se=sd(IC0)/sqrt(N)),
                                               by=byvars]
                            setnames(DT.B,"IC0","IF.Brier")
                        }
                        score.loob[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
                        score.loob[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
                    } else{
                        score.loob <- DT.B[,data.table(Brier=sum(residuals)/N), by=byvars]
                    }
                    data.table::setkeyv(score.loob,byvars)
                    ## data.table::setkey(DT.B,model,times)
                    ## DT.B <- DT.B[score.loob]
                    if (length(dolist)>0L){
                        if (se.fit==FALSE){
                            if (match("times",byvars,nomatch=0))
                                contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,model=model),
                                                                        NF=NF,
                                                                        N=N,
                                                                        alpha=alpha,
                                                                        dolist=dolist,
                                                                        se.fit=FALSE),by=list(times)]
                            else
                                contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,model=model),
                                                                        NF=NF,
                                                                        N=N,
                                                                        alpha=alpha,
                                                                        dolist=dolist,
                                                                        se.fit=FALSE)]
                        } else{
                            if (match("times",byvars,nomatch=0))
                                contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),
                                                                        NF=NF,
                                                                        N=N,
                                                                        alpha=alpha,
                                                                        dolist=dolist,
                                                                        se.fit=TRUE),by=list(times)]
                            else
                                contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),
                                                                        NF=NF,
                                                                        N=N,
                                                                        alpha=alpha,
                                                                        dolist=dolist,
                                                                        se.fit=TRUE)]
                        }
                        setnames(contrasts.Brier,"delta","delta.Brier")
                        output <- list(score=score.loob,contrasts=contrasts.Brier)
                    } else{
                        output <- list(score=score.loob)
                    }
                    if (keep.residuals) {
                        DT.B[,model:=factor(model,levels=mlevs,mlabels)]
                        output <- c(output,list(residuals=DT.B[,c("ID",names(response),"model","times","risk","residuals"),with=FALSE]))
                    }
                    if (!is.null(output$score)){
                        output$score[,model:=factor(model,levels=mlevs,mlabels)]
                    }
                    ## set model and reference in model comparison results
                    if (!is.null(output$contrasts)>0){
                        output$contrasts[,model:=factor(model,levels=mlevs,mlabels)]
                        output$contrasts[,reference:=factor(reference,levels=mlevs,mlabels)]
                    }
                    return(output)
                }
            
                # }}}
            })
            names(crossvalPerf) <- metrics
            # }}}
        }else{ ## split.method bootcv
            # {{{ bootcv
            if (parallel=="snow") exports <- c("DT.b","N.b","cens.model","multi.split.test") else exports <- NULL
            crossval <- foreach (j=1:B,.export=c("DT.b","N.b","cens.model","multi.split.test")) %dopar%{
                ## crossval <- lapply(1:B,function(j){
                DT.b <- DT.B[b==j]
                N.b <- length(unique(DT.b[["ID"]]))
                computePerformance(DT.b,
                                   N=N.b,
                                   se.fit=FALSE,
                                   conservative=TRUE, ## cannot subset IC yet
                                   cens.model=cens.model,
                                   multi.split.test=multi.split.test,
                                   keep.residuals=FALSE,
                                   keep.vcov=FALSE)
            }
            crossvalPerf <- lapply(metrics,function(m){
                if (response.type %in% c("survival","competing.risks")){
                    byvars <- c("model","times")
                } else{
                    byvars <- c("model")
                }
                ## score
                if (length(crossval[[1]][[m]]$score)>0){
                    cv.score <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$score}))
                    if (se.fit==TRUE){
                        bootcv.score <- cv.score[,data.table::data.table(mean(.SD[[m]],na.rm=TRUE),
                                                                         lower=quantile(.SD[[m]],alpha/2,na.rm=TRUE),
                                                                         upper=quantile(.SD[[m]],(1-alpha/2),na.rm=TRUE)),by=byvars,.SDcols=m]
                        data.table::setnames(bootcv.score,c(byvars,m,"lower","upper"))
                    } else{
                        bootcv.score <- cv.score[,data.table::data.table(mean(.SD[[m]],na.rm=TRUE)),by=byvars,.SDcols=m]
                        data.table::setnames(bootcv.score,c(byvars,m))
                    }
                }else{
                    cv.score <- NULL
                    bootcv.score <- NULL
                }
                ## contrasts and multi-split test
                if (length(crossval[[1]][[m]]$contrasts)>0){
                    cv.contrasts <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$contrasts}))
                    delta.m <- paste0("delta.",m)
                    bootcv.contrasts <- switch(as.character(se.fit+3*multi.split.test),
                                               "4"={
                                                   cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE),
                                                                                        lower=quantile(.SD[[delta.m]],alpha/2,na.rm=TRUE),
                                                                                        upper=quantile(.SD[[delta.m]],(1-alpha/2),na.rm=TRUE),
                                                                                        p=median(.SD[["p"]],na.rm=TRUE)),by=c(byvars,"reference"),.SDcols=c(delta.m,"p")]
                                               },
                                               "1"={cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE),
                                                                                         lower=quantile(.SD[[delta.m]],alpha/2,na.rm=TRUE),
                                                                                         upper=quantile(.SD[[delta.m]],(1-alpha/2),na.rm=TRUE)),by=c(byvars,"reference"),.SDcols=c(delta.m)]
                                               },
                                               "3"={cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE),
                                                                                         p=median(.SD[["p"]],na.rm=TRUE)),by=c(byvars,"reference"),.SDcols=c(delta.m,"p")]
                                               },
                                               "0"={
                                                   bootcv.contrasts <- cv.contrasts[,data.table::data.table(mean(.SD[[delta.m]],na.rm=TRUE)),by=c(byvars,"reference"),.SDcols=delta.m]
                                               })
                    data.table::setnames(bootcv.contrasts,"V1",delta.m)
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

    # }}}
    # {{{ collect data for calibration plots
    if ("Calibration" %in% plots){
        if (keep.residuals && split.method$name=="LeaveOneOutBoot"){
            crossvalPerf[["Calibration"]]$plotframe <- crossvalPerf$Brier$Residuals[model!=0,]
        } else{
            ## there are no residuals in this case. residuals are only available for LOOB!
            if (cens.type=="uncensored"){
                crossvalPerf[["Calibration"]]$plotframe <- DT.B[model!=0,data.table::data.table(response=ReSpOnSe,risk=mean(risk)),by=c(byvars,"ID")]
                setcolorder(crossvalPerf[["Calibration"]]$plotframe,c("ID","ReSpOnSe","model","risk"))
            }else{
                rnames <- names(Response)
                rnames <- rnames[rnames!="ID"]
                DT.mean <- DT.B[model!=0,data.table::data.table(risk=mean(risk)),by=c(byvars,"ID")]
                crossvalPerf[["Calibration"]]$plotframe <- merge(DT.mean,Response,by="ID")
                setcolorder(crossvalPerf[["Calibration"]]$plotframe,c("ID","model","times",rnames,"risk")) 
            }
        }
        crossvalPerf[["Calibration"]]$plotframe[,model:=factor(model,levels=mlevs,mlabels)]
        if (keep.residuals==FALSE && split.method$name=="LeaveOneOutBoot"){
            crossvalPerf$Brier$Residuals <- NULL
        }
        if (cens.type=="rightCensored")
            crossvalPerf[["Calibration"]]$plotframe <- merge(jack,crossvalPerf[["Calibration"]]$plotframe,by=c("ID","times"))
    }
    # }}}
}
    # }}}
    #------------------output-----------------------------------
    # {{{ enrich the output object
    
    if (split.method$internal.name=="noplan"){
        if (keep.residuals==TRUE){
            noSplit$Brier$residuals[,model:=factor(model,levels=mlevs,mlabels)]
        }
        output <- noSplit
    } else{
        output <- crossvalPerf
    }
    models <- mlevs
    names(models) <- mlabels
    if (keep.vcov){
        lab.models <- mlabels
        if (null.model){
            names(lab.models) <- paste0("model=",0:(length(mlabels)-1))
        }
        else{
            names(lab.models) <- paste0("model=",1:(length(mlabels)))
        }
        if (!is.null(output$Brier$vcov))
            attr(output$Brier$vcov,"models") <- lab.models
        if (!is.null(output$AUC$vcov))
            if (response.type=="binary"){
                if (NCOL(output$AUC$vcov)==length(mlabels)){
                    colnames(output$AUC$vcov) <- mlabels
                    rownames(output$AUC$vcov) <- mlabels
                }else{
                    colnames(output$AUC$vcov) <- mlabels[-1]
                    rownames(output$AUC$vcov) <- mlabels[-1]
                }
            }else{
                for (ml in 1:length(mlabels)){
                    colnames(output$AUC$vcov) <- gsub(paste0("model=",ml),paste0("model=",mlabels[ml]),colnames(output$AUC$vcov))
                    rownames(output$AUC$vcov) <- gsub(paste0("model=",ml),paste0("model=",mlabels[ml]),rownames(output$AUC$vcov))}
            }
        attr(output$AUC$vcov,"models") <- lab.models
    }
    if (null.model==TRUE) nm <- names(models)[1] else nm <- NULL
    if (!keep.splitindex) split.method$index <- "split index was not saved"
    output <- c(output,list(response.type=response.type,
                            dolist=dolist,
                            cause=cause,
                            states=states,
                            null.model=nm,
                            models=models,
                            cens.type=cens.type,
                            censoringHandling=cens.method,
                            split.method=split.method,
                            metrics=metrics,
                            alpha=alpha,
                            plots=plots,
                            summary=summary,
                            missing.predictions=missing.predictions,
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

Brier.binary <- function(DT,
                         se.fit,
                         conservative=FALSE,
                         cens.model="none",
                         keep.vcov=FALSE,
                         multi.split.test,
                         alpha,
                         N,
                         NT,
                         NF,
                         dolist,
                         keep.residuals=FALSE,
                         ...){
    residuals=Brier=risk=model=ReSpOnSe=lower=upper=se=ID=IF.Brier=NULL
    DT[,residuals:=(ReSpOnSe-risk)^2]
    if (se.fit==1L){
        data.table::setkey(DT,model,ID)
        score <- DT[,data.table::data.table(Brier=sum(residuals)/N,se=sd(residuals)/sqrt(N)),by=list(model)]
        score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
        score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        if (keep.vcov){
            DT[,Brier:=sum(residuals)/N,by=list(model)]
            DT[,IF.Brier:=residuals-Brier]
        }
    }else{
        ## no se.fit
        score <- DT[,data.table(Brier=sum(residuals)/N),by=list(model)]
    }
    data.table::setkey(score,model)
    if (length(dolist)>0){
        ## merge with Brier score
        data.table::setkey(DT,model)
        data.table::setkey(score,model)
        DT <- DT[score]
        if (se.fit==TRUE || multi.split.test==TRUE){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,
                                                             IF=residuals,
                                                             model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=multi.split.test,
                                                  se.fit=se.fit)]
        }else{
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=FALSE,
                                                  se.fit=FALSE)]
        }
        setnames(contrasts.Brier,"delta","delta.Brier")
        output <- list(score=score,contrasts=contrasts.Brier)
    }else{
        output <- list(score=score)
    }
    if (keep.vcov && se.fit==TRUE){
        output <- c(output,list(vcov=getVcov(DT,"IF.Brier")))
    }
    if (keep.residuals) {
        output <- c(output,list(residuals=DT[,data.table::data.table(ID,ReSpOnSe,model,risk,residuals)]))
    }
    output
}

Brier.survival <- function(DT,MC,se.fit,conservative,cens.model,keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,keep.residuals=FALSE,...){
    IC0=nth.times=ID=time=times=raw.Residuals=risk=Brier=residuals=WTi=Wt=status=setorder=model=IF.Brier=data.table=sd=lower=qnorm=se=upper=NULL
    ## compute 0/1 outcome:
    DT[time<=times & status==1,residuals:=(1-risk)^2/WTi]
    DT[time<=times & status==0,residuals:=0]
    DT[time>times,residuals:=(risk)^2/Wt]
    if (se.fit==1L || multi.split.test==TRUE){
        ## data.table::setorder(DT,model,times,time,-status)
        data.table::setorder(DT,model,times,ID)
        DT[,nth.times:=as.numeric(factor(times))]
        DT[,IC0:=residuals-mean(residuals),by=list(model,times)]
        if (conservative==TRUE){
            score <- DT[,data.table(Brier=sum(residuals)/N,
                                    se=sd(IC0)/sqrt(N)),by=list(model,times)]
            setnames(DT,"IC0","IF.Brier")
        }else{
            DT[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                                  time=time,
                                                  IC0,
                                                  residuals=residuals,
                                                  WTi=WTi,
                                                  Wt=Wt,
                                                  IC.G=MC,
                                                  cens.model=cens.model,
                                                  nth.times=nth.times[1]),by=list(model,times)]
            score <- DT[,data.table(Brier=sum(residuals)/N,
                                    se=sd(IF.Brier)/sqrt(N),
                                    se.conservative=sd(IC0)/sqrt(N)),by=list(model,times)]
        }
        if (se.fit==TRUE){
            score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        }
    }else{
        ## no se.fit
        score <- DT[,data.table(Brier=sum(residuals)/N),by=list(model,times)]
    }
    data.table::setkey(score,model,times)
    if (length(dolist)>0L){
        ## merge with Brier score
        data.table::setkey(DT,model,times)
        ## data.table::setkey(score,model,times)
        DT <- DT[score]
        if (se.fit==TRUE || multi.split.test==TRUE){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=multi.split.test,
                                                  se.fit=se.fit),by=list(times)]
        }else{
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=FALSE,
                                                  se.fit=FALSE),by=list(times)]
        }
        setnames(contrasts.Brier,"delta","delta.Brier")
        output <- list(score=score,contrasts=contrasts.Brier)
    } else{
        output <- list(score=score)
    }
    if (keep.vcov && se.fit==TRUE){
        output <- c(output,list(vcov=getVcov(DT,"IF.Brier",times=TRUE)))
    }
    if (keep.residuals) {
        output <- c(output,list(residuals=DT[,c("ID","time","status","model","times","risk","residuals"),with=FALSE]))
    }
    output
}

Brier.competing.risks <- function(DT,MC,se.fit,conservative,cens.model,keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,keep.residuals=FALSE,cause,states,...){
    IC0=nth.times=ID=time=times=event=Brier=raw.Residuals=risk=residuals=WTi=Wt=status=setorder=model=IF.Brier=data.table=sd=lower=qnorm=se=upper=NULL
    ## compute 0/1 outcome:
    thecause <- match(cause,states,nomatch=0)
    if (length(thecause)==0) stop("Cannot identify cause of interest")
    DT[time<=times & status==1 & event==thecause,residuals:=(1-risk)^2/WTi]
    DT[time<=times & status==1 & event!=thecause,residuals:=(risk)^2/WTi]
    DT[time<=times & status==0,residuals:=0]
    DT[time>times,residuals:=(risk)^2/Wt]
    ## deal with censored observations before t
    DT[time<=times & status==0,residuals:=0]
    if (se.fit==1L || multi.split.test==TRUE){
        ## data.table::setorder(DT,model,times,time,-status)
        data.table::setorder(DT,model,times,ID)
        DT[,nth.times:=as.numeric(factor(times))]
        DT[,IC0:=residuals-mean(residuals),by=list(model,times)]
        if (conservative){
            score <- DT[,data.table(Brier=sum(residuals)/N,
                                    se=sd(IC0)/sqrt(N)),by=list(model,times)]
        }else{
            DT[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                                  time=time,
                                                  IC0,
                                                  residuals=residuals,
                                                  WTi=WTi,
                                                  Wt=Wt,
                                                  IC.G=MC,
                                                  cens.model=cens.model,
                                                  nth.times=nth.times[1]),by=list(model,times)]
            score <- DT[,data.table(Brier=sum(residuals)/N,
                                    se=sd(IF.Brier)/sqrt(N),
                                    se.conservative=sd(IC0)/sqrt(N)),by=list(model,times)]
        }
        if (se.fit==TRUE){
        score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
        score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        }
    }else{
        ## no se.fit
        score <- DT[,data.table(Brier=sum(residuals)/N),by=list(model,times)]
    }
    data.table::setkey(score,model,times)
    if (length(dolist)>0){
        data.table::setkey(DT,model,times)
        ## merge with Brier score
        DT <- DT[score]
        data.table::setkey(score,model,times)
        if (se.fit==TRUE || multi.split.test==TRUE){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=multi.split.test,
                                                  se.fit=se.fit),by=list(times)]
        }else{
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=FALSE,
                                                  se.fit=FALSE),by=list(times)]
        }
        setnames(contrasts.Brier,"delta","delta.Brier")
        output <- list(score=score,contrasts=contrasts.Brier)
    } else{
        output <- list(score=score)
    }
    if (keep.residuals) {
        output <- c(output,list(residuals=DT[,c("ID","time","status","model","times","risk","residuals"),with=FALSE]))
    }
    if (keep.vcov && se.fit==TRUE){
        output <- c(output,list(vcov=getVcov(DT,"IF.Brier",times=TRUE)))
    }
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
delongtest <-  function(risk,
                        score,
                        dolist,
                        response,
                        cause,
                        alpha,
                        multi.split.test,
                        se.fit,
                        keep.vcov) {
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
                if (se.fit==TRUE){
                    ## cor.auc[ctr] <- S[i, j]/sqrt(S[i, i] * S[j, j])
                    LSL <- t(c(1, -1)) %*% S[c(j, i), c(j, i)] %*% c(1, -1)
                    ## print(c(1/LSL,rms::matinv(LSL)))
                    se[ctr] <- sqrt(LSL)
                }
                ## tmpz <- (delta.AUC[ctr]) %*% rms::matinv(LSL) %*% delta.AUC[ctr]
                ## tmpz <- (delta.AUC[ctr]) %*% (1/LSL) %*% delta.AUC[ctr]
                model[ctr] <- modelnames[j]
                reference[ctr] <- modelnames[i]
                ctr <- ctr + 1
            }
        }
        if (se.fit==TRUE||multi.split.test==TRUE){
            deltaAUC <- data.table(model,reference,delta.AUC=as.vector(delta.AUC),se)
            deltaAUC[,p:=2*pnorm(abs(delta.AUC)/se,lower.tail=FALSE)]
        }else{
            deltaAUC <- data.table(model,reference,delta.AUC=as.vector(delta.AUC))
        }
        if (se.fit==TRUE){ 
            deltaAUC[,lower:=delta.AUC-Qnorm*se]
            deltaAUC[,upper:=delta.AUC+Qnorm*se]
        }
        out <- list(score = score, contrasts = deltaAUC)
    }else{
        out <- list(score = score, contrasts = NULL)
    }
    if (keep.vcov) {
        out <- c(out,list(vcov=S))
    }
    out
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

AUC.binary <- function(DT,breaks=NULL,se.fit,conservative=FALSE,cens.model="none",keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,ROC,...){
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
    if (length(dolist)>0 || (se.fit==1L)){
        xRisk <- data.table::dcast(aucDT[],ID~model,value.var="risk")[,-1,with=FALSE]
        delong.res <- delongtest(risk=xRisk,
                                 score=output$score,
                                 dolist=dolist,
                                 response=aucDT[model==model[1],ReSpOnSe],
                                 cause="1",
                                 alpha=alpha,
                                 multi.split.test=multi.split.test,
                                 se.fit=se.fit,
                                 keep.vcov=keep.vcov)
        output$score <- delong.res$score
        output$contrasts <- delong.res$contrasts
        if (keep.vcov){
            output$vcov <- delong.res$vcov
        }
        output
    }else{
        output
    }
}


AUC.survival <- function(DT,MC,se.fit,conservative,cens.model,keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,ROC,...){
    ID=model=times=risk=Cases=time=status=Controls=TPR=FPR=WTi=Wt=ipcwControls=ipcwCases=IF.AUC=lower=se=upper=AUC=NULL
    cause <- 1
    aucDT <- DT[model>0]
    ## remove null model comparisons
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
    if (se.fit==1L || multi.split.test==TRUE){
        ## compute influence function
        ## data.table::setorder(aucDT,model,times,time,-status)
        data.table::setorder(aucDT,model,times,ID)
        aucDT[,IF.AUC:=getInfluenceCurve.AUC.survival(t=times[1],
                                                      n=N,
                                                      time=time,
                                                      risk=risk,
                                                      Cases=Cases,
                                                      Controls=Controls,
                                                      ipcwControls=ipcwControls,
                                                      ipcwCases=ipcwCases,
                                                      MC=MC), by=list(model,times)]
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
        if (keep.vcov){
            output <- c(output,list(vcov=getVcov(aucDT,"IF.AUC",times=TRUE)))
        }
    }
    ## add score to object
    output <- c(list(score=score),output)
    if (length(dolist)>0){
        if (se.fit==TRUE || multi.split.test==TRUE){
            contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),
                                                   NF=NF,
                                                   N=N,
                                                   alpha=alpha,
                                                   dolist=dolist,multi.split.test=multi.split.test,se.fit=se.fit),by=list(times)]
        }else{
            contrasts.AUC <- score[,getComparisons(data.table(x=AUC,model=model),
                                                   NF=NF,
                                                   N=N,
                                                   alpha=alpha,
                                                   dolist=dolist,
                                                   multi.split.test=FALSE,                                                   
                                                   se.fit=FALSE),by=list(times)]
        }
        setnames(contrasts.AUC,"delta","delta.AUC")
        output <- c(list(score=score,contrasts=contrasts.AUC),output)
    }
    output
}

AUC.competing.risks <- function(DT,MC,se.fit,conservative,cens.model,keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,cause,states,ROC,...){
    ID=model=times=risk=Cases=time=status=event=Controls1=Controls2=TPR=FPR=WTi=Wt=ipcwControls1=ipcwControls2=ipcwCases=IF.AUC=lower=se=upper=AUC=NULL
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
    thecause <- match(cause,states,nomatch=0)
    if (length(thecause)==0) stop("Cannot identify cause of interest")
    aucDT[,Cases:=(time <=times &  event==thecause)]
    aucDT[,Controls1:=(time > times)] 
    aucDT[,Controls2:=(time <=times &  event!=thecause & status !=0)]
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
    if (se.fit==1L || multi.split.test==TRUE){
        ## compute influence function
        ## data.table::setorder(aucDT,model,times,time,-status)
        data.table::setorder(aucDT,model,times,ID)
        aucDT[,IF.AUC:=getInfluenceCurve.AUC.competing.risks(t=times[1],
                                                             n=N,
                                                             time=time,
                                                             risk=risk,
                                                             ipcwControls1=ipcwControls1,
                                                             ipcwControls2=ipcwControls2,
                                                             ipcwCases=ipcwCases,
                                                             Cases=Cases,
                                                             Controls1=Controls1,
                                                             Controls2=Controls2,
                                                             MC=MC), by=list(model,times)]
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
        if (keep.vcov){
            output <- c(output,list(vcov=getVcov(aucDT,"IF.AUC",times=TRUE)))
        }
    }
    ## add score to object
    output <- c(list(score=score),output)
    if (length(dolist)>0){
        if (se.fit==TRUE || multi.split.test==TRUE){
            contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),
                                                   NF=NF,
                                                   N=N,
                                                   alpha=alpha,
                                                   dolist=dolist,multi.split.test=multi.split.test,
                                                   se.fit=se.fit),by=list(times)]
        }else{
            contrasts.AUC <- score[,getComparisons(data.table(x=AUC,model=model),
                                                   NF=NF,
                                                   N=N,
                                                   alpha=alpha,
                                                   dolist=dolist,
                                                   multi.split.test=FALSE,
                                                   se.fit=FALSE),by=list(times)]
        }
        setnames(contrasts.AUC,"delta","delta.AUC")
        output <- c(list(score=score,contrasts=contrasts.AUC),output)
    }
    output
}

# }}}
