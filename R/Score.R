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
##' @name Score
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
##'     statistics to apply to the predicted risks. Choices are \code{"risks"}, \code{"IPA"},
##'     \code{"riskQuantile"} and \code{"ibs"}. Can be all \code{c("risks","IPA","riskQuantile","ibs")} or a subset thereof.
##'     \itemize{
##'     \item \code{"risks"} adds the predicted risks to the output.
##'     \item \code{"ipa"} computes the index of prediction accuracy (AKA R-squared) based on Brier scores for model vs null model
##'     \item \code{"riskQuantile"} calculates
##'     time-point specific boxplots for the
##'     predicted risks (or biomarker values) conditional on the outcome at the time-point.
##'     \item \code{"ibs"} calculates integrated Brier scores across the time points at which the Brier score is computed. This works only with
##'     time-to-event outcome and the results depend on the argument \code{times}.
##'     }
##'     Set to \code{NULL} to avoid estimation of summary statistics.
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
##' @param progress.bar Style for \code{txtProgressBar}. Can be 1,2,3 see \code{help(txtProgressBar)} or NULL to avoid the progress bar.
##' @param keep list of characters (not case sensitive) which determines additional output.
##' \code{"residuals"} provides Brier score residuals and
##' \code{"splitindex"} provides sampling index used to split the data into training and validation sets.
##' \code{"vcov"} provides the variance-covariance matrix of the estimated parameters.
##' @param censoring.save.memory Only relevant in censored data where censoring weigths are obtained with
##' Cox regression and argument \code{conservative} is set to \code{FALSE}. If \code{TRUE}, save memory by not storing the influence function
##' of the cumulative hazard of the censoring as a matrix when calculating standard errors
##' with Cox censoring. This can allow one to use \code{Score} on larger test data sets,
##' but may be slower.
##' @param predictRisk.args
##'  A list of argument-lists to control how risks are predicted.
##'  The names of the lists should be the S3-classes of the \code{object}.
##'  The argument-lists are then passed on to the S3-class specific predictRisk method.
##'  For example, if your object contains one or several random forest model fitted with the function randomForestSRC::rfsrc then you can
##'  specify additional arguments for the function riskRegression::predictRisk.rfsrc which will pass
##'  these on to the function randomForestSRC::predict.rfsrc. A specific example in this case would be
##'  \code{list(rfsrc=list(na.action="na.impute"))}.
##'  
##'
##'  A more flexible approach is to write a new predictRisk S3-method. See Details.
##' @param errorhandling Argument passed as \code{.errorhandling} to foreach. Default is \code{"pass"}.
##' @param debug Logical. If \code{TRUE} indicate landmarks in progress of the program.
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
##' replacement from the data, which does not depend on the sample size:
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
##'          data=testdat,plots=c("calibration","ROC"))
##' \dontrun{plotROC(xb)
##' plotCalibration(xb)
##' }
##'
##' ## compute AUC for a list of continuous markers
##' markers = as.list(testdat[,.(X6,X7,X8,X9,X10)])
##' Score(markers,formula=Y~1,data=testdat,metrics=c("auc"))
##'
##' # cross-validation
##' \dontrun{
##'     set.seed(10)
##'     learndat=sampleData(400,outcome="binary")
##'     lr1a = glm(Y~X6,data=learndat,family=binomial)
##'     lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
##'     ## bootstrap cross-validation
##'     x1=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,split.method="bootcv",B=100)
##'     x1
##'     ## leave-one-out and leave-pair-out bootstrap
##'     x2=Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,
##'              split.method="loob",
##'              B=100,plots="calibration")
##'     x2
##' }
##' # survival outcome
##'
##' # Score Cox regression models
##' \dontrun{library(survival)
##' library(rms)
##' library(prodlim)
##' set.seed(18)
##' trainSurv <- sampleData(100,outcome="survival")
##' testSurv <- sampleData(40,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
##' x=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'          formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,times=c(5,8))
##' ## Use Cox to estimate censoring weights
##' y=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'          formula=Surv(time,event)~X1+X8,data=testSurv,conf.int=FALSE,times=c(5,8)) 
##' ## Use GLMnet to estimate censoring weights
##' z=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'          formula=Surv(time,event)~X1+X8,cens.model = "GLMnet",data=testSurv,
##'          conf.int=FALSE,times=c(5,8)) 
##' ## Use hal9001 to estimate censoring weights
##' w=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'          formula=Surv(time,event)~X1+X8,cens.model = "Hal9001",data=testSurv,
##'          conf.int=FALSE,times=c(5,8)) 
##' x
##' y
##' z
##' w
##' }
##'
##' \dontrun{library(survival)
##' library(rms)
##' library(prodlim)
##' set.seed(18)
##' trainSurv <- sampleData(100,outcome="survival")
##' testSurv <- sampleData(40,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
##' xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'          formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,times=c(5,8))
##' xs
##' }
##' # Integrated Brier score
##' \dontrun{
##' xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'          formula=Surv(time,event)~1,data=testSurv,conf.int=FALSE,
##'          summary="ibs",
##'          times=sort(unique(testSurv$time)))
##' }
##'
##' # time-dependent AUC for list of markers
##' \dontrun{survmarkers = as.list(testSurv[,.(X6,X7,X8,X9,X10)])
##' Score(survmarkers,
##'       formula=Surv(time,event)~1,metrics="auc",data=testSurv,
##'       conf.int=TRUE,times=c(5,8))
##'
##' # compare models on test data
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=testSurv,conf.int=TRUE,times=c(5,8))
##' }
##' # crossvalidation models in traindata
##' \dontrun{
##'     library(survival)
##'     set.seed(18)
##'     trainSurv <- sampleData(400,outcome="survival")
##'     cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
##'     cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv, y=TRUE, x = TRUE)
##'     x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'                formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
##'                split.method="loob",B=100,plots="calibration")
##'
##'     x2= Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'               formula=Surv(time,event)~1,data=trainSurv,conf.int=TRUE,times=c(5,8),
##'               split.method="bootcv",B=100)
##' }
##'
##' # restrict number of comparisons
##' \dontrun{
##'     Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'           formula=Surv(time,event)~1,data=trainSurv,contrasts=TRUE,
##'           null.model=FALSE,conf.int=TRUE,times=c(5,8),split.method="bootcv",B=3)
##'
##'     # competing risks outcome
##'     set.seed(18)
##'     trainCR <- sampleData(400,outcome="competing.risks")
##'     testCR <- sampleData(400,outcome="competing.risks")
##'     library(riskRegression)
##'     library(cmprsk)
##'     # Cause-specific Cox regression
##'     csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR)
##'     csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR)
##'     # Fine-Gray regression
##'     fgr1 = FGR(Hist(time,event)~X1+X2+X7+X9,data=trainCR,cause=1)
##'     fgr2 = FGR(Hist(time,event)~X3+X5+X6,data=trainCR,cause=1)
##'     Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2,
##'                "FGR(X1+X2+X7+X9)"=fgr1,"FGR(X3+X5+X6)"=fgr2),
##'           formula=Hist(time,event)~1,data=testCR,se.fit=1L,times=c(5,8))
##' }
##'
##'
##'
##' \dontrun{
##'     # reproduce some results of Table IV of Blanche et al. Stat Med 2013
##'     data(Paquid)
##'     ResPaquid <- Score(list("DSST"=-Paquid$DSST,"MMSE"=-Paquid$MMSE),
##'                        formula=Hist(time,status)~1,
##'                        data=Paquid,
##'                        null.model = FALSE,
##'                        conf.int=TRUE,
##'                        metrics=c("auc"),
##'                        times=c(3,5,10),
##'                        plots="ROC")
##'     ResPaquid
##'     plotROC(ResPaquid,time=5)
##' }
##' \dontrun{
##' # parallel options
##' # by erikvona: Here is a generic example of using future
##' # and doFuture, works great with the current version:
##' library(riskRegression)
##' library(future)
##' library(foreach)
##' library(doFuture)
##' library(survival)
##' # Register all available cores for parallel operation
##' plan(multiprocess, workers = availableCores())
##' registerDoFuture()
##' set.seed(10)
##' trainSurv <- sampleData(400,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv,
##'              y=TRUE, x = TRUE)
##' # Bootstrapping on multiple cores
##' x1 = Score(list("Cox(X1+X2+X7+X9)"=cox1),
##'      formula=Surv(time,event)~1,data=trainSurv, times=c(5,8),
##'      parallel = "as.registered", split.method="bootcv",B=100)
##' }
##'
##'
##'
##' @author Thomas A Gerds \email{tag@@biostat.ku.dk} and Paul Blanche \email{paul.blanche@@univ-ubs.fr}
##' @references
##'
##' Thomas A. Gerds and Michael W. Kattan (2021).
##' Medical Risk Prediction Models: With Ties to Machine Learning (1st ed.)
##' Chapman and Hall/CRC
##' https://doi.org/10.1201/9781138384484
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
#'
#' Michael W Kattan and Thomas A Gerds. The index of prediction accuracy: an
#' intuitive measure useful for evaluating risk prediction models. Diagnostic
#' and Prognostic Research, 2(1):7, 2018.
##'
##' @export Score.list
##' @export
##'
# }}}


##' @export
Score <- function(object,...){
    UseMethod("Score",object=object)
}

# {{{ Score.list
##' @rdname Score
##' @export Score.list
##' @export
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
                       parallel=c("no","multicore","snow","as.registered"),
                       ncpus=1,
                       cl=NULL,
                       progress.bar=3,
                       errorhandling="pass",
                       keep,
                       predictRisk.args,
                       debug=0L,
                       censoring.save.memory = FALSE,
                       ...){
    se.conservative=IPCW=IF.AUC.conservative=IF.AUC0=IF.AUC=IC0=Brier=AUC=casecontrol=se=nth.times=time=status=ID=WTi=risk=IF.Brier=lower=upper=crossval=b=time=status=model=reference=p=model=pseudovalue=ReSpOnSe=residuals=event=j=NULL

    # }}}
    theCall <- match.call()
    # {{{ decide about metrics and plots

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
    if (length(posIBS <- grep("^ibs$|^crps$",summary,ignore.case=TRUE))>0){
        summary <- summary[-posIBS]
        if (!("Brier" %in% metrics)) metrics <- c(metrics,"Brier")
        ibs <- TRUE
    }else{
        ibs <- FALSE
    }
    ## IPA
    if (length(posRR <- grep("^ipa$|^rr$|^r2|rsquared$",summary,ignore.case=TRUE))>0){
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
    # {{{ censoring model arguments
    if (length(grep("^km|^kaplan|^marg",cens.model,ignore.case=TRUE))>0){
        cens.model <- "KaplanMeier"
    } else{
        if (length(attr(terms(formula),"factors"))>0){
            # cens.model <- "cox"
        } else{
            cens.model <- "KaplanMeier"
        }
    }
    ### 
    if (cens.model != "cox" && cens.model != "none" && cens.model != "KaplanMeier" && cens.model != "discrete"){
      if (!conservative[[1]]){
        warning("For this model, we can't calculate the IF of the Survival function of the censoring distribution. \n Therefore, force conservative = TRUE")
        conservative[[1]] <- TRUE
      }
      passed.args <- names(as.list(match.call())[-1])
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
    if (data.table::is.data.table(data))
        data <- copy(data)
    else
        data <- data.table::setDT(data)
    responseFormula <- stats::update(formula,~1)
    ## if (missing(event)) event <- 1
    responsevars <- all.vars(responseFormula)
    if (!all(responsevars%in%names(data)))stop("Response variable(s) ",paste(responsevars,collapse=", ")," not found in data.")
    response <- getResponse(formula=responseFormula,
                            data=data,
                            cause=cause,
                            vars=responsevars)
    if (any(is.na(response))) stop("Missing values in response.
Delete missing values *before* calling `Score', if this is reasonable,
c.f., Chapter 7, Section 5 in Gerds & Kattan 2021. Medical risk prediction models. ")
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
    ## but first rename to avoid problems with pre-existing variables with same name(s)
    ## data[,eval(responsevars):=NULL]
    ## data have to be ordered when ipcw is called
    if ("event" %in% names(data)){
        setnames(data,"event","protectedName.event")
    }
    if (response.type=="binary"){
        data[,event:=factor(ReSpOnSe,
                            levels=1:(length(states)+1),
                            labels=c(states,"event-free"))]
    }
    if (response.type %in% c("survival","competing.risks")){
        test.names <- match(colnames(response),names(data),nomatch=0)
        if (sum(test.names)>0){
            protected.names <- names(data)[test.names]
            setnames(data,protected.names,paste0("protectedName.",protected.names))
        }
        data <- cbind(response,data)
        N <- as.numeric(NROW(data))
        neworder <- data[,order(time,-status)]
        data.table::setorder(data,time,-status)
    }else{
        data <- cbind(response,data)
        N <- as.numeric(NROW(data))
        neworder <- 1:N
    }
    ## add ID variable for merging purposes and because output has long format
    ## data[["ID"]]=1:N
    data[,ID:=1:N]
    if (response.type=="survival")
        formula <- stats::update(formula,"prodlim::Hist(time,status)~.")
    if (response.type=="competing.risks")
        formula <- stats::update(formula,"prodlim::Hist(time,event)~.")
    ## stop("Dont know how to predict response of type ",response.type))
    cens.type <- attr(response,"cens.type")
    if (is.null(cens.type)) cens.type <- "uncensored"
    rm(response)
    # }}}
    # {{{ SplitMethod & parallel stuff

    if (!missing(seed)) {
        ## message("Random seed set to control split of data: seed=",seed)
        set.seed(seed)
    }
    split.method <- getSplitMethod(split.method=split.method,B=B,N=N,M=M,seed=seed)
    B <- split.method$B
    splitIndex <- split.method$index
    do.resample <- !(is.null(splitIndex))
    if (split.method$internal.name!="noplan"){
        if (missing(parallel)) parallel <- "no"
        parallel <- match.arg(parallel)
        ## Additionally, I'd like the process to be adjusted when a cluster is
        ## passed or `as.registered` is specified to ignore the `ncpus` argument
        ## (very counter-intuitive that when I'm setting up the cluster I still
        ## need to specify ncpus > 1):
        if (ncpus<=1 && is.null(cl) && parallel != "as.registered") {
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
        if (split.method$name[1]=="BootCv" && multi.split.test[1]==TRUE){
            if ("AUC" %in% metrics) {
                warning("Cannot do multi-split test with AUC yet. Forced multi.split.test=FALSE")
                multi.split.test=FALSE
            }else{
                if (se.fit==TRUE & conservative==FALSE){
                    warning("Cannot do deal with conservative=FALSE when also multi.split.test=TRUE. Forced conservative=TRUE.")
                }
                conservative=TRUE
            }
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
            if(inherits(x=try(fit$call,silent=TRUE),what="try-error")||is.null(fit$call))
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
        if (is.logical(conf.int)[[1]] && conf.int[[1]]==FALSE
            || conf.int[[1]]<=0
            || conf.int[[1]]>1)
            se.fit <- 0L
        else
            se.fit <- 1L
    }else{
        stopifnot(is.logical(se.fit)||se.fit[[1]]==0||se.fit[[1]]==1)
    }
    if (split.method$internal.name=="noplan") multi.split.test <- FALSE
    if (se.fit[1]==1L) {
        if (is.numeric(conf.int) && conf.int[1]<1 && conf.int[1]>0)
            alpha <- 1-conf.int
        else
            alpha <- .05
    }else{
        alpha <- NA
    }
    # allow user to write contrasts=0 instead of contrasts=FALSE
    if (length(contrasts)==1 && length(contrasts[[1]])==1 && contrasts==0) contrasts = FALSE
    if ((NF+length(nullobject))<=1) dolist <- NULL
    else{
        if ((is.logical(contrasts) && contrasts[1]==FALSE)){
            dolist <- NULL
        } else{
            if (is.logical(contrasts) && contrasts[1]==TRUE){
                if (is.null(nullobject)){
                    dolist <- lapply(1:(NF-1),function(i){c(i,(i+1):NF)})
                } else{
                    dolist <- lapply(0:(NF-1),function(i){c(i,(i+1):NF)})
                }
            }else{
                dolist <- contrasts
                if (!is.list(contrasts)) contrasts <- list(contrasts)
                if (!(all(sapply(dolist,function(x){
                    if (!(length(x)>1)) stop("All elements of the list contrasts must contain at least two elements, i.e., the two models to be compared.")
                    all(x<=NF) && all(x>=0)
                }))))
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
        if (!missing(times) && (!is.null(times)) && (times[1]!=FALSE)) warning("Function 'Score': Response type is not time-to-event: argument 'times' will be ignored.",call.=FALSE)
        times <- NULL
        NT <- 1
    }
    # }}}
    # {{{ Dealing with censored data outside the loop
    if (response.type %in% c("survival","competing.risks")){
        if (cens.type=="rightCensored"){
            getIC <- se.fit[[1]] && !conservative[[1]] && cens.model != "KaplanMeier"
            Weights <- getCensoringWeights(formula=formula,
                                           data=data,
                                           times=times,
                                           cens.model=cens.model,
                                           response.type=response.type,
                                           influence.curve=getIC, 
                                           censoring.save.memory = censoring.save.memory)
            ##split.method$internal.name %in% c("noplan",".632+")
            ## if cens.model is marginal then IC is a matrix (ntimes,newdata)
            ## if cens.model is Cox then IC is an array (nlearn, ntimes, newdata)
            ## IC is an array with dimension (nlearn, times, newdata)
            ## IC_G(t,z;x_k)
        }else {
            if (cens.type=="uncensored"){
                cens.method <- "none"
                cens.model <- "none"
                Weights <- NULL
            }
            else{
                stop("Cannot handle this type of censoring.")
            }
        }
        if ("pseudo" %in% cens.method){
            if (cens.type[1]=="rightCensored"
                && (response.type %in% c("survival","competing.risks"))){
                # need to communicate the censoring code of the event variable
                # produced by Hist via model.frame in case of competing risks
                if (response.type=="competing.risks"){
                    censcode <- data[status==0,event[1]]
                    margForm <- Hist(time,event,cens.code=censcode)~1
                }else{
                    censcode <- data[status==0,status[1]]
                    margForm <- update(formula,".~1")
                }
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
    # {{{ Nosplit performance: external data (hopefully not apparent)
    missing.predictions <- "Don't know yet"
    off.predictions <- "Don't know yet"
    if (split.method$internal.name %in% c("noplan",".632+")){
        DT <- getPerformanceData(testdata=data,
                                 testweights=Weights,
                                 traindata=NULL,
                                 trainseed=NULL,
                                 response.type=response.type,
                                 response.dim=response.dim,
                                 times=times,
                                 cause=cause,
                                 neworder=neworder,
                                 debug=debug,
                                 levs=mlevs,
                                 labels=mlabels,
                                 predictRisk.args=predictRisk.args,
                                 nullobject=nullobject,
                                 cens.type=cens.type,
                                 object=object,
                                 object.classes=object.classes,
                                 NT=NT
)
        if (any(is.na(DT[["risk"]]))){
            missing.predictions <- DT[,list("Missing.values"=sum(is.na(risk))),by=byvars]
            missing.predictions[,model:=factor(model,levels=mlevs,mlabels)]
            warning("Missing values in the predicted risks. See `missing.predictions' in output list.")
        }else{
            missing.predictions <- "None"
        }
        if (("Brier"%in% metrics) && (any(is.na(DT[["risk"]]))|| (max(DT[["risk"]])>1 || min(DT[["risk"]])<0))){
            off.predictions <- DT[,list("missing.values"=sum(is.na(risk)),"negative.values"=sum(risk<0,na.rm=TRUE),"values.above.1"=sum(risk>1,na.rm=TRUE)),by=byvars]
            off.predictions[,model:=factor(model,levels=mlevs,mlabels)]
            warning("Predicted values off the probability scale (negative or above 100%). See `off.predictions' in output list.\nOnly a problem for the Brier score, You can stop this warning by setting metrics='auc'.")
        }else{
            off.predictions <- "None"
        }
        noSplit <- computePerformance(DT=DT,
                                      N=N,
                                      NT=NT,
                                      NF=NF,
                                      models=list(levels=mlevs,labels=mlabels),
                                      response.type=response.type,
                                      times=times,
                                      jack=jack,
                                      cens.type=cens.type,
                                      cause=cause,
                                      states=states,
                                      alpha=alpha,
                                      se.fit=se.fit,
                                      conservative=conservative,
                                      cens.model=cens.model,
                                      multi.split.test=multi.split.test,
                                      keep.residuals=keep.residuals,
                                      keep.vcov=keep.vcov,
                                      dolist=dolist,
                                      probs=probs,
                                      metrics=metrics,
                                      plots=plots,
                                      summary=summary,
                                      ibs=ibs,
                                      ipa = ipa,
                                      ROC=FALSE,
                                      MC=Weights$IC,
                                      IC.data=Weights$IC.data)
        if (debug) message("computed apparent performance")
    }
    # }}}
    # {{{ Crossvalidation
    # {{{ bootstrap re-fitting and k-fold-CV

if (split.method$internal.name%in%c("BootCv","LeaveOneOutBoot","crossval")){
    if (missing(trainseeds)||is.null(trainseeds)){
        if (!missing(seed)) set.seed(seed)
        if (split.method$internal.name == "crossval"){
            trainseeds <- sample(1:1000000,size=B*split.method$k,replace=FALSE)
        }else{
            trainseeds <- sample(1:1000000,size=B,replace=FALSE)
        }
    }
    if (parallel=="snow") exports <- c("data","split.method","Weights","N","trainseeds") else exports <- NULL
    if (!is.null(progress.bar)){
        message("Running crossvalidation algorithm")
        if (!(progress.bar %in% c(1,2,3))) progress.bar <- 3
        if (B==1 && split.method$internal.name == "crossval"){
            message(paste0("Fitting the models in ",split.method$k," learning datasets, then predicting the risks in validation datasets"))
            pb <- txtProgressBar(max = split.method$k, style = progress.bar,width=20)
        } else{
            message(paste0("Fitting the models in ",B," learning datasets, then predicting the risks in validation datasets"))
            pb <- txtProgressBar(max = B, style = progress.bar,width=20)
        }
    }
    `%dopar%` <- foreach::`%dopar%`
    ## k-fold-CV
    if (split.method$internal.name == "crossval"){
        DT.B <- foreach::foreach (b=1:B,.export=exports,.packages="data.table",.errorhandling=errorhandling) %dopar%{
            ## repetitions of k-fold to avoid Monte-Carlo error
            index.b <- split.method$index(b) ## contains a sample of the numbers 1:k with replacement
            if((B>1) && !is.null(progress.bar)){setTxtProgressBar(pb, b)}
            DT.b <- rbindlist(lapply(1:split.method$k,function(fold){
                traindata=data[index.b!=fold]
                testids <- index.b==fold # (1:N)[index.b!=fold]
                if((B==1) && !is.null(progress.bar)){
                    setTxtProgressBar(pb, fold)
                }
                ## NOTE: subset.data.table preserves order ## So we need to use subset??
                testdata <- subset(data,testids)
                if (cens.type=="rightCensored"){
                    testweights <- Weights
                    # Need to check what's expected that testids is here and below:
                    testweights$IPCW.subject.times <- subset(testweights$IPCW.subject.times,testids)
                    if (Weights$dim>0){
                        testweights$IPCW.times <- subset(testweights$IPCW.times,testids)
                    }
                } else {
                    testweights <- NULL
                }
                ## predicted risks of model trained without this fold
                ## evaluated and added to this fold
                DT.fold <- getPerformanceData(testdata=testdata,
                                              testweights=testweights,
                                              traindata=traindata,
                                              trainseed=trainseeds[[(b-1)*split.method$k+fold]],
                                              response.type=response.type,
                                              response.dim=response.dim,
                                              times=times,
                                              cause=cause,
                                              neworder=NULL,
                                              debug=debug,
                                              levs=mlevs,
                                              labels=mlabels,
                                              predictRisk.args=predictRisk.args,
                                              nullobject=nullobject,
                                              cens.type=cens.type,
                                              object=object,
                                              object.classes=object.classes,
                                              NT=NT)
                return(DT.fold)
            }))
            DT.b[,b:=b]
            DT.b
        }
    }else{# either LeaveOneOutBoot or BootCv
        DT.B <- foreach::foreach (b=1:B,.export=exports,.packages="data.table",.errorhandling=errorhandling) %dopar%{
            if(!is.null(progress.bar)){setTxtProgressBar(pb, b)}
            ## DT.B <- rbindlist(lapply(1:B,function(b){
            traindata=data[split.method$index(b)]
            ## setkey(traindata,ID)
            testids <- (match(1:N,unique(split.method$index(b)),nomatch=0)==0)
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
            ## print(class(traindata))
            ## print(names(traindata))
            DT.b <- getPerformanceData(testdata=testdata,
                                       testweights=testweights,
                                       traindata=traindata,
                                       trainseed=trainseeds[[b]],
                                       response.type=response.type,
                                       response.dim=response.dim,
                                       times=times,
                                       cause=cause,
                                       neworder=NULL,
                                       debug=debug,
                                       levs=mlevs,
                                       labels=mlabels,
                                       predictRisk.args=predictRisk.args,
                                       nullobject=nullobject,
                                       cens.type=cens.type,
                                       object=object,
                                       object.classes=object.classes,
                                       NT=NT)
            DT.b[,b:=b]
            DT.b
        }
    }
    trycombine <- try(DT.B <- rbindlist(DT.B),silent=TRUE)
    if (inherits(trycombine,"try-error")){
      # error handling (because rbindlist is useless at reporting errors)
      for (b in 1:B){
        if ("error" %in% class(DT.B[[b]]) || "simpleError" %in% "error" %in% class(DT.B[[b]])){
          stop(paste0("Errors occured in training models for crossvalidation. \n ", "Error found in iteration b = ",b,": \n ", DT.B[[b]]))
        }
      }
    }
    else if (nrow(DT.B)==0){
      stop("All training models failed. ")
    }
    
    if (!is.null(progress.bar)){
        cat("\n")
    }
    if (any(is.na(DT.B[["risk"]]))){
        missing.predictions <- DT.B[,list("Missing.values"=sum(is.na(risk))),by=byvars]
        warning("Missing values in the predicted risk. See `missing.predictions' in output list.")
    }
    if (("Brier"%in% metrics)&& (max(DT.B[["risk"]])>1 || min(DT.B[["risk"]])<0)){
        off.predictions <- DT.B[,list("negative.values"=sum(risk<0),"values.above.1"=sum(risk>1)),by=byvars]
        warning("Values off the scale (either negative or above 100%) in the predicted risk. See `off.predictions' in output list.")
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
    Response.names <- names(Response)
    Response.names <- Response.names[Response.names!="ID"]
    setkey(Response,ID)
    ## ## Show format for the data in DT.B
    ## cat(paste("\nDT.B for method:", split.method$name, "\n"))
    ## print(DT.B)
    # }}}
    # {{{ Leave-one-out bootstrap
    ## start clause split.method$name=="LeaveOneOutBoot
    if (split.method$internal.name =="crossval" && B == 1){
      crossvalPerf<-computePerformance(DT=DT.B,
                                       N=N,
                                       NT=NT,
                                       NF=NF,
                                       models=list(levels=mlevs,labels=mlabels),
                                       response.type=response.type,
                                       times=times,
                                       jack=jack,
                                       cens.type=cens.type,
                                       cause=cause,
                                       states=states,
                                       alpha=alpha,
                                       se.fit=se.fit,
                                       conservative=conservative,
                                       cens.model=cens.model,
                                       multi.split.test=multi.split.test,
                                       keep.residuals=FALSE,
                                       keep.vcov=FALSE,
                                       dolist=dolist,
                                       probs=probs,
                                       metrics=metrics,
                                       plots=plots,
                                       summary=summary,
                                       ibs=ibs,
                                       ipa=ipa,
                                       ROC=FALSE,
                                       MC=Weights$IC,
                                       IC.data=Weights$IC.data)
    }
    else if (split.method$name=="LeaveOneOutBoot" | split.method$internal.name =="crossval"){  ## Testing if the crossval works in this loop
        message(paste0("Calculating the performance metrics in long format\nlevel-1 data with ",
                       NROW(DT.B),
                       " rows.",
                       ifelse(NROW(DT.B)>1000000,
                              " This may take a while ...",
                              " This should be fast ...")))
      crossvalPerf <- lapply(metrics, function(m){crossvalPerf.loob(m,
                                                                    times,
                                                                    mlevs,
                                                                    se.fit,
                                                                    response.type,
                                                                    NT,
                                                                    Response,
                                                                    cens.type,
                                                                    Weights,
                                                                    split.method,
                                                                    N,
                                                                    B,
                                                                    DT.B,
                                                                    data,
                                                                    dolist,
                                                                    alpha,
                                                                    byvars,
                                                                    mlabels,
                                                                    ipa,
                                                                    keep.residuals,
                                                                    conservative,
                                                                    cens.model,
                                                                    response.dim,
                                                                    ID,
                                                                    cause)})
      names(crossvalPerf) <- metrics
    }
    if (split.method$name=="BootCv"){
        # {{{ bootcv
        if (parallel=="snow") exports <- c("DT.B","N.b","cens.model","multi.split.test") else exports <- NULL
        if (!is.null(progress.bar)){
            if (!(progress.bar %in% c(1,2,3))) progress.bar <- 3
            pb1 <- txtProgressBar(max = B, style = progress.bar,width=20)
            message(paste0("Calculating the performance metrics in ",B," validation data sets"))
        }
        crossval <- foreach::foreach(j=1:B,.export=exports,.packages="data.table",.errorhandling=errorhandling) %dopar%{
            DT.b <- DT.B[b==j]
            N.b <- length(unique(DT.b[["ID"]]))
            if(!is.null(progress.bar)){
                setTxtProgressBar(pb1, j)
            }
            computePerformance(DT=DT.b,
                               N=N.b,
                               NT=NT,
                               NF=NF,
                               models=list(levels=mlevs,labels=mlabels),
                               response.type=response.type,
                               times=times,
                               jack=jack,
                               cens.type=cens.type,
                               cause=cause,
                               states=states,
                               alpha=alpha,
                               se.fit=FALSE,
                               conservative=TRUE,
                               cens.model=cens.model,
                               multi.split.test=multi.split.test,
                               keep.residuals=FALSE,
                               keep.vcov=FALSE,
                               dolist=dolist,
                               probs=probs,
                               metrics=metrics,
                               plots=plots,
                               summary=summary,
                               ibs=ibs,
                               ipa=ipa,
                               ROC=FALSE,
                               MC=Weights$IC,
                               IC.data=Weights$IC.data)
        }
        if (!is.null(progress.bar)){
            cat("\n")
        }
        crossvalPerf <- lapply(metrics,function(m) crossvalPerf.bootcv(m,crossval,se.fit,multi.split.test,keep.cv,byvars,alpha))
        names(crossvalPerf) <- metrics
        if (ipa==TRUE){
            if (response.type=="binary")
                crossvalPerf[["Brier"]][["score"]][,IPA:=1-Brier/Brier[model=="Null model"]]
            else
                crossvalPerf[["Brier"]][["score"]][,IPA:=1-Brier/Brier[model=="Null model"],by=times]
        }
    }
    # }}}
    # {{{ collect data for plotRisk
    if ("risks"%in% summary){
        crossvalPerf[["risks"]]$score <- DT.B[,{
            c(.SD[1,Response.names,with=FALSE],
              list(risk=mean(risk),
                   sd.risk=sd(risk),
                   oob=.N))},.SDcols=c(Response.names,"risk"),by=c(byvars,"ID")]
        crossvalPerf[["risks"]]$score[,model:=factor(model,levels=mlevs,mlabels)]
        setcolorder(crossvalPerf[["risks"]]$score,c("ID",byvars,Response.names,"risk","sd.risk"))
    }
    # }}}
    # {{{ collect data for calibration plots
    if ("Calibration" %in% plots){
        if (keep.residuals[[1]]==TRUE && split.method$name[[1]]=="LeaveOneOutBoot"){
            crossvalPerf[["Calibration"]]$plotframe <- crossvalPerf$Brier$Residuals[model!=0,]
        } else{
            ## there are no residuals in this case. residuals are only available for LOOB!
            crossvalPerf[["Calibration"]]$plotframe <- DT.B[model!=0,{
                c(.SD[1,Response.names,with=FALSE],
                  list(risk=mean(risk),
                       oob=.N))},.SDcols=c(Response.names,"risk"),by=c(byvars,"ID")]
            setcolorder(crossvalPerf[["Calibration"]]$plotframe,c("ID",byvars,Response.names,"risk"))
        }
        crossvalPerf[["Calibration"]]$plotframe[,model:=factor(model,levels=mlevs,mlabels)]
        if (keep.residuals[[1]]==FALSE && split.method$name[[1]]=="LeaveOneOutBoot"){
            crossvalPerf$Brier$Residuals <- NULL
        }
        if (cens.type=="rightCensored")
            crossvalPerf[["Calibration"]]$plotframe <- merge(jack,crossvalPerf[["Calibration"]]$plotframe,by=c("ID","times"))
    }
    # }}}
}

    # }}}
    # {{{ the output object
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
                            times=times,
                            alpha=alpha,
                            plots=plots,
                            summary=summary,
                            missing.predictions=missing.predictions,
                            off.predictions=off.predictions,
                            call=theCall))
    # remove null model from AUC
    # if (!is.null(output$AUC) && null.model){
    #   output$AUC$score <- output$AUC$score[model!="Null model"]
    #   output$AUC$contrasts <- output$AUC$contrasts[reference!="Null model"]
    # }
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
