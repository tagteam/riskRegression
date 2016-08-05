##' Method to Score risk markers and risk prediction models
##'
##' We compute the Brier score and the area under the ROC curve. For
##' survival possibly with competing risk both are time-dependent.
##' @title Score risk predictions
##' @aliases Score
##' @param object List of risk predictions (see details and examples).
##' @param formula A formula which identifies the outcome (left hand
##'     side). For right censored outcome, the right hand side of the
##'     formula is used to estimate the IPCW model.
##' @param data Data set or table in which the formula can be
##'     interpreted.
##' @param metrics Character vector specifying which metrics to
##'     apply. Implemented are \code{"auc"} and \code{"Brier"}.
##' @param summary  Character vector specifying which summary statistics to
##'     apply to the predicted risks. Implemented is \code{"riskQuantile"}.
##' @param plots Character vector specifying which plots to prepare.
##' @param cause Event of interest. Used for binary outcome \code{Y}
##'     to specify that risks are risks of the event \code{Y=event}
##'     and for competing risks outcome to specify the cause of
##'     interest.
##' @param times For survival and competing risks outcome: list of
##'     prediction horizons. All times which are greater than the
##'     maximal observed time in the data set are removed.
##' @param landmarks Not yet implemented.
##' @param useEventTimes If \code{TRUE} add all unique event times to
##'     argument \code{times}.
##' @param nullModel If \code{TRUE} add the null model which ignores
##'     the covariates and predicts the prevalence for all subjects.
##' @param test If \code{TRUE} compute confidence intervals. Also do
##'     model comparisons specified in \code{dolist}.
##' @param alpha Level of significance.
##' @param dolist Vector of integers specifying which risks to compare
##'     to all other risks in \code{object}.
##' @param probs Quantiles for retrospective summary statistics of the predicted risks 
##' @param censMethod Method for dealing with right censored
##'     data. Either \code{"ipcw"} or \code{"pseudo"}.
##' @param censModel Model for estimating inverse probability of
##'     censored weights.
##' @param splitMethod Method for cross-validation.
##' @param B Number of cross-validation steps.
##' @param M Size of subsamples for cross-validation. If specified it
##'     has to be an integer smaller than the size of \code{data}.
##' @param seed Super seed for setting training data seeds when randomly splitting the data for cross-validation.
##' @param trainseeds Seeds for training models during cross-validation.
##' @param ... Not used
##' @return Data table with scores and tests.
##' @examples
##' # binary outcome
##' library(lava)
##' set.seed(18)
##' learndat <- sampleData(100,outcome="binary")
##' testdat <- sampleData(40,outcome="binary")
##'
##' # score logistic regression models
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
##' lr2 = glm(Y~X3+X5+X6,data=learndat,family=binomial)
##' Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,data=testdat)
##' 
##' # compute AUC for a list of continuous markers
##' markers = as.list(testdat[,1:5])
##' u=Score(markers,formula=Y~1,data=testdat,metrics=c("auc"))
##'
##' # cross-validation
##' lr1a = glm(Y~X6,data=learndat,family=binomial)
##' lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
##' Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,splitMethod="bootcv",B=3)
##'
##' # survival outcome
##' 
##' # Score Cox regression models
##' library(survival)
##' library(rms)
##' library(prodlim)
##' set.seed(18)
##' trainSurv <- sampleData(100,outcome="survival")
##' testSurv <- sampleData(40,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv)
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=testSurv,test=FALSE,times=c(5,8))
##'
##' # time-dependent AUC for list of markers
##' survmarkers = as.list(testSurv[,1:5])
##' Score(survmarkers,
##'       formula=Surv(time,event)~1,metrics="auc",data=testSurv,
##' test=TRUE,times=c(5,8))
##' 
##' # compare models on test data
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=testSurv,test=TRUE,times=c(5,8))
##'
##' # crossvalidation models in traindata
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=trainSurv,test=TRUE,times=c(5,8),
##' splitMethod="bootcv",B=3)
##'
##' # restrict number of comparisons
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=trainSurv,dolist=1,
##' nullModel=FALSE,test=TRUE,times=c(5,8),splitMethod="bootcv",B=3)
##'
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
##'       formula=Hist(time,event)~1,data=testCR,test=TRUE,times=c(5,8))
##' 
##' @author Thomas A Gerds \email{tag@@biostat.ku.dk} and Paul Blanche \email{paul.blanche@@univ-ubs.fr}
##'
##' @export 
Score.list <- function(object,
                       formula,
                       data,
                       metrics=c("auc","brier"),
                       summary=NULL,
                       plots=c("roc","calibration","pvalues"),
                       cause=1,
                       times,
                       landmarks,
                       useEventTimes=FALSE,
                       nullModel=TRUE,
                       test=TRUE,
                       alpha=0.05,
                       dolist,
                       probs=c(0.05,0.25,0.5,0.75,0.95),
                       censMethod="ipcw",
                       censModel="cox",
                       splitMethod,
                       B,
                       M,
                       seed,
                       trainseeds,
                       ...){
    ## G=c(0,0.25,0.5,0.75,1),
    ## summary=match.arg(summary,c("riskQuantile","riskClass"))
    id=time=status=id=WTi=b=time=status=model=reference=p=model=NULL
    # -----------------parse arguments and prepare data---------
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
    response <- getResponse(formula=responseFormula,
                            data=data,
                            cause=cause,
                            vars=responsevars)
    response.dim <- NCOL(response)
    responseType <- attr(response,"model")
    # add null model and find names for the object
    if (nullModel==TRUE){
        nullobject <- getNullModel(formula=formula,data=data,responseType=responseType)
    } else{
        nullobject <- NULL
    }
    ## put ReSpOnSe for binary and (time, event, status) in the first column(s) 
    ## data[,eval(responsevars):=NULL]
    data <- cbind(response,data)
    if (responseType=="survival")
        formula <- stats::update(formula,"Hist(time,status)~.")
    if (responseType=="competing.risks")
        formula <- stats::update(formula,"Hist(time,event)~.")
    N <- NROW(response)
    ## predictHandlerFun <- switch(responseType,
                                ## "survival"="predictRisk",
                                ## "competing.risks"="predictRisk",
                                ## "binary"="predictRisk",
                                ## stop("Dont know how to predict response of type ",responseType))
    censType <- attr(response,"cens.type")
    if (is.null(censType)) censType <- "uncensoredData"
    # }}}
    # {{{ SplitMethod
    splitMethod <- getSplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
    B <- splitMethod$B
    splitIndex <- splitMethod$index
    do.resample <- !(is.null(splitIndex))
    # }}}
    # {{{ Checking the models
    # for predictHandlerFunction
    allmethods <- utils::methods(predictRisk)
    ## wantedMethods <- lapply(object,function(o){
    ## candidateMethods <- paste(predictHandlerFun,class(o),sep=".")
    ## if (all(match(candidateMethods,allmethods,nomatch=0)==0))
    ## stop(paste("Could not find ",predictHandlerFun," method for ",paste(class(o),collapse=" ,"),sep=""))
    ## })
    # checking the models for compatibility with resampling
    if (is.null(names(object))){
        names(object) <- sapply(object,function(o)class(o)[1])}
    else {
        names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
    }
    names.object <- names(object) <- make.unique(names(object))
    NF <- length(object)
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
            ## else fit$call$data <- NULL
            ## fit
        })
        ## names(object) <- names.object
    }
    if ((NF+length(nullobject))<=1) dolist <- NULL 
    if (test==FALSE) dolist <- NULL
    ## test <- FALSE
    if (test==TRUE && missing(dolist)){
        if (is.null(nullobject)) {
            dolist <- 1:(NF-1)
        } else{
            dolist <- 0:(NF-1)
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
        if (missing(landmarks)){
            start <- 0
            if (missing(times)){
                if (useEventTimes==TRUE)
                    times <- unique(c(start,eventTimes))
                else{
                    ## times <- seq(start,maxtime,(maxtime - start)/100)
                    times <- median(eventTimes)
                }
            } else{
                if (useEventTimes==TRUE) 
                    times <- sort(c(start,unique(times),eventTimes))
                else
                    ## times <- sort(unique(c(start,times)))
                    times <- sort(unique(times))
            }
            stopifnot(sum(times<=maxtime)>0)
            times <- times[times<=maxtime]
            NT <-  length(times)
        }
        else{
            stop("Landmark updating not yet implemented.")
        }
    }
    else{
        if (!missing(times)) warning("Function 'Score': Response type is not time-to-event: argument 'times' will be ignored.",call.=FALSE)
        times <- NULL
        NT <- 1
    }
    # }}}
    # ----------------------------find metrics and plots ----------------------
    # {{{
    
    ## Metrics <- lapply(metrics,grep,c("AUC","Brier"),ignore.case=TRUE,value=TRUE)
    metrics[grep("^auc$",metrics,ignore.case=TRUE)] <- "AUC"
    metrics[grep("^brier$",metrics,ignore.case=TRUE)] <- "Brier"
    plots[grep("^roc$",plots,ignore.case=TRUE)] <- "ROC"
    plots[grep("^cal",plots,ignore.case=TRUE)] <- "Cal"
    ## Plots <- lapply(plots,grep,c("Roc","Cal"),ignore.case=TRUE,value=TRUE)
    if ("ROC" %in% plots) {
        ## add AUC if needed
        if (!("AUC" %in% metrics)) metrics <- c(metrics,"AUC")
    }
    # }}}
    # -----------------IPCW outside loop -----------------------
    # {{{ 
    if (responseType %in% c("survival","competing.risks")){
        if (censType=="rightCensored"){
            if ("outside" %in% censMethod){
                if ("ipcw" %in% censMethod){
                    Weights <- getCensoringWeights(formula=formula,
                                                   data=testdata,
                                                   response=response,
                                                   times=times,
                                                   censModel=censModel,
                                                   responseType=responseType)
                }
                else{
                    censMethod <- "jackknife.pseudo.values"
                    pseudoResponse <- getPseudoValues(formula=formula,data=testdata,responseType=responseType,times=times,cause=cause)
                }
            }else{
                censMethod <- c(censMethod,"inside")
            }
        }
        else{
            if (censType=="uncensored"){
                censMethod <- c("none","inside")
                Weights <- NULL
            }
            else{
                stop("Cannot handle this type of censoring.")
            }
        }
    }
    # }}}
    # -----------------performance program----------------------
    # {{{ 
    trainModel <- function(model,data){
        model$call$data <- data
        try(eval(model$call),silent=TRUE)
    }
    computePerformance <- function(object,
                                   nullobject,
                                   testdata,
                                   traindata=NULL,
                                   trainseed=NULL,
                                   metrics,
                                   plots,
                                   times,
                                   cause,
                                   test,
                                   alpha,
                                   dolist,
                                   NF,
                                   NT){
        N <- NROW(testdata)
        # split into response and predictors
        response <- testdata[,1:response.dim,with=FALSE]
        response[,ID:=testdata[["ID"]]]
        X <- testdata[,-c(1:response.dim),with=FALSE]
        # {{{ -----------------IPCW inner loop-----------------------
        if (responseType %in% c("survival","competing.risks")){
            if ("inside" %in% censMethod){
                if (censType=="rightCensored"){
                    if ("ipcw" %in% censMethod){
                        Weights <- getCensoringWeights(formula=formula,
                                                       data=testdata,
                                                       response=response,
                                                       times=times,
                                                       censModel=censModel,
                                                       responseType=responseType)
                        ## add subject specific weights here, and time specific weights later
                        response[,WTi:=Weights$IPCW.subjectTimes]
                    }else{
                        censMethod <- "jackknife.pseudo.values"
                        pseudoResponse <- getPseudoValues(formula=formula,
                                                          data=testdata,
                                                          responseType=responseType,
                                                          times=times,
                                                          cause=cause)
                    }
                } else{
                    if (censType=="uncensored"){
                        Weights <- list(IPCW.times=rep(1,NT),
                                        IPCW.subjectTimes=matrix(1,ncol=NT,nrow=N))
                        Weights$method <- "marginal"
                        response[,WTi:=1]
                    } else{
                        stop("Cannot handle this type of censoring.")
                    }
                }
            }
        } else{
            Weights <- NULL
        }
        # }}}
        # extract predictions as data.table
        args <- switch(responseType,"binary"={list(newdata=X)},
                       "survival"={list(newdata=X,times=times)},
                       "competing.risks"={list(newdata=X,times=times,cause=cause)},
                       stop("Unknown responseType."))
        pred <- data.table::rbindlist(lapply(mlevs, function(f){
            if (f!=0 && any(c("integer","factor","numeric","matrix") %in% class(object[[f]]))){
                if (is.null(dim(object[[f]])))
                    p <- c(object[[f]][testdata[["ID"]]])
                else
                    p <- c(object[[f]][testdata[["ID"]]])
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
                        ## browser(skipCalls=TRUE)
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
        if (!is.null(Weights)){
            if (Weights$method=="marginal"){
                Wt <- data.table(times=times,Wt=Weights$IPCW.times)
            } else {
                Wt <- data.table(times=rep(times,rep(N,NT)),Wt=Weights$IPCW.times)
            }
            pred <- merge(pred,Wt,by="times")
        }
        # compute and test performance
        ## input <- list("prediction"=pred,response=response)
        input <- list(DT=merge(response,pred,by="ID"),
                      N=N,
                      NT=NT,
                      NF=NF,
                      alpha=alpha,
                      test=test,
                      dolist=dolist,probs=probs,ROC=FALSE)
        if (responseType=="competing.risks")
            input <- c(input,list(cause=cause))
        if (responseType %in% c("survival","competing.risks") &&test==TRUE){
            if (censType=="rightCensored"){
                input <- c(input,list(MC=response[,getInfluenceCurve.KM(time=time,status=status)]))
            }
            else{
                input <- c(input,list(MC=matrix(0,ncol=length(response$time),nrow=length(unique(response$time)))))
            }
        }
        out <- vector(mode="list",length=length(c(summary,metrics,plots)))
        names(out) <- c(summary,metrics,plots)
        for (s in summary){
            out[[s]] <- do.call(paste(s,responseType,sep="."),input)
            out[[s]]$score[,model:=factor(model,levels=mlevs,mlabels)]
            if (NROW(out[[s]]$test)>0){
                out[[s]]$test[,model:=factor(model,levels=mlevs,mlabels)]
                out[[s]]$test[,reference:=factor(reference,levels=mlevs,mlabels)]
            }
        }
        for (m in metrics){
            if (m=="AUC" && ("ROC" %in% plots)){
                input <- replace(input, "ROC",TRUE)
                out[[m]] <- do.call(paste(m,responseType,sep="."),input)
                out[["ROC"]]$plotframe <- out[[m]]$ROC
                out[["ROC"]]$plotframe[,model:=factor(model,levels=mlevs,mlabels)]
                out[[m]]$ROC <- NULL
            }else{
                input <- replace(input, "ROC",FALSE)
                out[[m]] <- do.call(paste(m,responseType,sep="."),input)
            }
            out[[m]]$score[,model:=factor(model,levels=mlevs,mlabels)]
            if (NROW(out[[m]]$test)>0){
                out[[m]]$test[,model:=factor(model,levels=mlevs,mlabels)]
                out[[m]]$test[,reference:=factor(reference,levels=mlevs,mlabels)]
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
                                  metrics=metrics,
                                  plots=plots,
                                  times=times,
                                  cause=cause,
                                  test=test,
                                  alpha=alpha,
                                  dolist=dolist,
                                  NF=NF,
                                  NT=NT)
    # }}}
    # -----------------crossvalidation performance---------------------
    # {{{ 
    crossval <- NULL
    if (splitMethod$name=="BootCv"){
        if (missing(trainseeds)||is.null(trainseeds)){
            if (!missing(seed)) set.seed(seed)
            trainseeds <- sample(1:1000000,size=B,replace=FALSE)
        }
        crossval <- foreach (b=1:B) %dopar%{
            traindata=data[splitMethod$index[,b],,drop=FALSE]
            ## subset.data.table preserves order
            testdata <- subset(data,(match(1:N,unique(splitMethod$index[,b]),nomatch=0)==0),drop=FALSE)
            cb <- computePerformance(object=object,
                                     nullobject=nullobject,
                                     testdata=testdata,
                                     traindata=traindata,
                                     trainseed=trainseeds[b],
                                     metrics=metrics,
                                     plots=plots,
                                     times=times,
                                     cause=cause,
                                     test=test,
                                     alpha=alpha,
                                     dolist=dolist,
                                     NF=NF,
                                     NT=NT)
            cb
        }
        ## browser(skipCalls=1)
        ## names(crossval[[1]])
        crossvalPerf <- lapply(metrics,function(m){
            ## if (test==TRUE){
            if (length(crossval[[1]][[m]]$score)>0){
                bootcv <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$score}))
                if (length(crossval[[1]][[m]]$test)>0){
                    if (responseType %in% c("survival","competing.risks")){
                        multisplit.test <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$test[,data.table(times,model,reference,p)]}))
                    }else{
                        multisplit.test <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$test[,data.table(model,reference,p)]}))
                    }
                }else{ 
                    multisplit.test <- NULL
                }
            }else{
                multisplit.test <- NULL
                bootcv <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]}))
            }
            if (responseType %in% c("survival","competing.risks")){
                bootcv <- bootcv[,data.table::data.table(mean(eval(as.name(m))),lower=quantile(eval(as.name(m)),alpha/2),upper=quantile(eval(as.name(m)),(1-alpha/2))),by=list(model,times)]
                data.table::setnames(bootcv,c("model","times",m,paste0(m,c(".lower",".upper"))))
            } else{
                bootcv <- bootcv[,data.table::data.table(mean(eval(as.name(m))),lower=quantile(eval(as.name(m)),alpha/2),upper=quantile(eval(as.name(m)),(1-alpha/2))),by=list(model)]
                data.table::setnames(bootcv,c("model",m,paste0(m,c(".lower",".upper"))))
            }
            out <- list(bootcv)
            names(out) <- "score"
            ## names(out) <- paste0("Cross-validation (average of ",B," steps)")
            if (!is.null(multisplit.test)){
                if (match("times",colnames(multisplit.test),nomatch=0))
                    ms <- list(test=multisplit.test[,data.table(p=median(p)),by=list(model,reference,times)])
                else
                    ms <- list(test=multisplit.test[,data.table(p=median(p)),by=list(model,reference)])
                out <- c(out,ms)
            }
            out
        })
        ## names(crossvalPerf) <- names(crossval[[1]])
        names(crossvalPerf) <- metrics
        ## Brier <- data.table::rbindlist(lapply(crossval,function(x)x[["Brier"]]))
        ## Brier[,list(looboot=mean(ipcwResiduals)),by=list(model,times,ID)]
        ## bootcv=Brier[,list(bootcv=mean(ipcwResiduals)),by=list(model,times,b)]
    }
    # }}}
    #------------------output-----------------------------------
    if (is.null(crossval))
        output <- noSplit
    else{
        ## output <- list(noSplitPerf=noSplit,crossValPerf=crossvalPerf)
        output <- crossvalPerf
    }
    # -----------------specific task handlers-------------------
    # {{{
    models <- mlevs
    names(models) <- mlabels
    output <- c(output,list(responseType=responseType,
                            dolist=dolist,
                            models=models,
                            censType=censType,
                            censoringHandling=censMethod,
                            splitMethod=splitMethod,metrics=metrics,plots=plots,
                            summary=summary))
    
    for (p in c(plots)){
        output[[p]]$plotmethod <- p
        class(output[[p]]) <- c("scoreROC")
    }
    for (m in c(metrics,summary)){
        output[[m]]$metrics <- m
        class(output[[m]]) <- paste0("score",m)
    }
    class(output) <- "Score"
    output
}

##' @export
Score <- function(object,...){
  UseMethod("Score",object=object)
}


#' Boxplot risk quantiles
#' 
#' Retrospective boxplots of risk quantiles conditional on outcome
#' @param x Score object obtained by calling function \code{Score}.
#' @param model Choice of risk prediction model
#' @param reference Choice of reference risk prediction model for
#'     calculation of risk differences.
#' @param type Either \code{"risk"} for predicted risks or
#'     \code{"diff"} for differences between predicted risks.
#' @param timepoint time point specifying the prediction horizon
#' @param lwd line width
#' @param xlim x-axis limits
#' @param xlab x-axis label
#' @param main title of plot
#' @param ... not used
##' @examples
##' db=sampleData(100,outcome="binary")
##' fitconv=glm(Y~X3+X5,data=db,family=binomial)
##' fitnew=glm(Y~X1+X3+X5+X6+X7,data=db,family=binomial)
##' scoreobj=Score(list(new=fitnew,conv=fitconv),formula=Y~1,
##'                data=db,summary="riskQuantile",nullModel=FALSE)
##' boxplot(scoreobj)
##' 
##' library(survival)
##' ds=sampleData(100,outcome="survival")
##' fitconv=coxph(Surv(time,event)~X3+X5,data=ds)
##' fitnew=coxph(Surv(time,event)~X1+X3+X5+X6+X7,data=ds)
##' scoreobj=Score(list(conv=fitconv,new=fitnew),formula=Hist(time,event)~1,
##'                data=ds,summary="riskQuantile",times=5,nullModel=FALSE)
##' boxplot(scoreobj)
##'
##' library(riskRegression)
##' data(Melanoma)
##' fitconv = CSC(Hist(time,status)~invasion+age+sex+logthick,data=Melanoma)
##' fitnew = CSC(Hist(time,status)~invasion+age+sex,data=Melanoma)
##' scoreobj=Score(list(conv=fitconv,new=fitnew),formula=Hist(time,status)~1,
##'                data=Melanoma,summary="riskQuantile",times=5*365.25,nullModel=FALSE)
##' boxplot(scoreobj)
#' @export
boxplot.Score <- function(x,model,reference,type,timepoint,lwd=3,xlim,xlab,main,...){
    times=cause=models=NULL
    fitted <- x$models
    models <- names(x$models)
    if (missing(type)) {
        if (length(models)==1) 
            type="risk"
        else
            type=ifelse(NROW(x$riskQuantile$test)>0,"diff","risk")
    }
    if (type=="diff"){
        pframe <- x$riskQuantile$test
    } else{
        pframe <- x$riskQuantile$score
    }
    if (x$responseType!='binary'){
        if (missing(timepoint))
            timepoint <- max(pframe[["times"]])
        else ## can only do one timepoint
            timepoint <- timepoint[[1]]
        pframe <- pframe[times==timepoint]
    }
    if(missing(model)) mod <- pframe[,model[1]] else mod <- model
    if (type=="diff"){
        if(missing(reference)) ref <- pframe[,reference[1]] else ref <- reference
        pframe <- pframe[model==model & reference==reference]
    }else{
        pframe <- pframe[model==model ]
    }
    qq.pos <- grep("^Q",colnames(pframe))
    if (missing(xlim))
        if (type=="risk"){xlim=c(0,100) 
        } else {
            max <- ceiling(max(abs(100*(pframe[,qq.pos,with=FALSE]))))
            xlim=c(-max,max)
        }
    if (missing(main))
        if (type=="risk") main=mod else main="Difference in predicted risks"
    if (missing(xlab))
        if (type=="risk") xlab="Predicted risk" else xlab=""
    if (x$responseType!="competing.risks"){
        plot(0,0,type="n",
             main=main,
             xlim = xlim,
             ylim = c(0,NROW(pframe)),
             axes=FALSE,
             xlab = xlab,
             ylab = "")
        axis(1,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),labels=paste(seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),"%"))
        text(x=xlim[1],y=c(0.5,1.5,2.5,3),labels=c("Overall","Event","Event-free",expression(bold(Outcome))),pos=2,xpd=NA)
        if (type=="diff"){
            mtext(paste(ref,"higher risk"),side=1,adj=0,line=par()$mgp[1])
            mtext(paste(mod,"higher risk"),side=1,adj=1,line=par()$mgp[1])
        }
        bxp(list(stats=t(100*pframe[cause=="overall",qq.pos,with=FALSE,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=0.5,axes=FALSE)
        bxp(list(stats=t(100*pframe[cause=="event",qq.pos,with=FALSE,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=1.5,axes=FALSE)
        bxp(list(stats=t(100*pframe[cause=="event-free",qq.pos,with=FALSE,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=2.5,axes=FALSE)
    }else{
        plot(0,0,type="n",
             main=main,
             xlim = xlim,
             ylim = c(0,NROW(pframe)),
             axes=FALSE,
             xlab = xlab,
             ylab = "")
        ## axis(1,at=seq(0,100,25),labels=paste(seq(0,100,25),"%"))
        axis(1,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),labels=paste(seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),"%"))
        causes <- pframe[,cause]
        ypos <- c((1:(length(causes)))-0.5,length(causes))
        text(x=xlim[1],y=ypos,labels=c(causes,expression(bold(Outcome))),pos=2,xpd=NA)
        if (type=="diff"){
            mtext(paste(ref,"higher risk"),side=1,adj=0,line=par()$mgp[1])
            mtext(paste(mod,"higher risk"),side=1,adj=1,line=par()$mgp[1])
        }
        for (i in 1:length(causes)){
            cc <- causes[[i]]
            bxp(list(stats=t(100*pframe[cause==cc,qq.pos,with=FALSE,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=ypos[i],axes=FALSE)
        }
    }
    invisible(x)
}

plot.riskQuantile <- function(x,text.title="Outcome after 10 years",xlab="",text=rownames(x),...){
    plot(0,0,type="n",axes=FALSE,xlim=c(-10,10),ylim=c(1,5),xlab=xlab,ylab="")
    axis(1,at=c(-10,-5,-2.5,0,2.5,5,10))
    bxp(list(stats=t(100*x[4,,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=1.5,axes=FALSE)
    bxp(list(stats=t(100*x[3,,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=2.5,axes=FALSE)
    bxp(list(stats=t(100*x[2,,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=3.5,axes=FALSE)
    bxp(list(stats=t(100*x[1,,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=4.5,axes=FALSE)
    text(x=-10,y=c(1.5,2.5,3.5,4.5),labels=rev(text),xpd=NA)
    text(x=-10,y=5,labels=text.title,xpd=NA)
}




##' Plot Brier curve
##'
##' @title Plot Brier curve
#' @export 
#' @param x Object obtained with \code{Score.list}
##' @param models Choice of models to plot
##' @param lwd Line width
##' @param xlim Limits for x-axis 
##' @param ylim Limits for y-axis 
##' @param axes Logical. If \code{TRUE} draw axes.
##' @param ... Not yet used
plot.score.Brier <- function(x,models,lwd=3,xlim,ylim,axes=TRUE,...){
    times=model=Brier=dimcol=lower.Brier=upper.Brier=NULL
    pframe <- x$Brier$score
    if (missing(xlim)) xlim <- pframe[,range(times)]
    if (missing(ylim)) ylim <- c(0,.25)
    plot(0,0,type="n",ylim = ylim,
         xlim = xlim,
         axes=FALSE,
         xlab = "Time",
         ylab = "Brier score")
    if (axes){
        axis(1)
        prodlim::PercentAxis(2,at=seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5))
    }
    if (!missing(models)) pframe <- pframe[model %in% models]
    pframe[,col:=as.numeric(as.factor(model))]
    pframe[,lwd:=lwd]
    pframe[,lines(times,Brier,type="l",lwd=lwd,col=col),by=model]
    pframe[,dimcol:=prodlim::dimColor(col[[1]],density=55),by=model]
    pframe[,polygon(x=c(times,rev(times)),y=c(lower.Brier,rev(upper.Brier)),col=dimcol,border=NA),by=model]
}




Brier.binary <- function(DT,test=TRUE,alpha,N,NT,NF,dolist,...){
    residuals=risk=model=ReSpOnSe=lower.Brier=upper.Brier=se.Brier=NULL
    DT[,residuals:=(ReSpOnSe-risk)^2,by=model]
    score <- DT[,list(Brier=mean(residuals)),by=model]
    if (test==TRUE){
        data.table::setorder(DT,model,ReSpOnSe)
        score <- DT[,data.table::data.table(Brier=sum(residuals)/N,se.Brier=sd(residuals)/sqrt(N)),by=list(model)]
        score[,lower.Brier:=pmax(0,Brier-qnorm(1-alpha/2)*se.Brier)]
        score[,upper.Brier:=Brier + qnorm(1-alpha/2)*se.Brier]
        data.table::setkey(score,model)
        data.table::setkey(DT,model)
        DT <- DT[score]
        test.Brier <- DT[,getComparisons(data.table(x=Brier,IC=residuals,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist)]
        Brier <- list(score=score,test=test.Brier)
    }else{
         Brier <- list(score=DT[,list(Brier=mean(residuals)),by=list(model)])
     }
    Brier
}




## do not want to depend on Daim as they turn auc > 0.5
delongtest <-  function(risk, score, dolist, response, cause, alpha) {
    cov=lower=upper=p=AUC=se.AUC=lower.AUC=upper.AUC=NULL
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
    q1 <- auc/(2 - auc)
    q2 <- 2 * auc^2/(1 + auc)
    aucvar <- (auc * (1 - auc) + (nCases - 1) * (q1 - auc^2) + (nControls - 1) * (q2 - auc^2))/(nCases * nControls)
    ncomp <- nauc * (nauc - 1)/2
    delta.auc <- numeric(ncomp) 
    se.auc <- numeric(ncomp)
    model <- numeric(ncomp)
    reference <- numeric(ncomp)
    ctr <- 1
    Qnorm <- qnorm(1 - alpha/2)
    for (i in 1:(nauc - 1)) {
        for (j in (i + 1):nauc) {
            delta.auc[ctr] <- auc[j]-auc[i]
            ## cor.auc[ctr] <- S[i, j]/sqrt(S[i, i] * S[j, j])
            LSL <- t(c(1, -1)) %*% S[c(j, i), c(j, i)] %*% c(1, -1)
            ## print(c(1/LSL,rms::matinv(LSL)))
            se.auc[ctr] <- sqrt(LSL)
            ## tmpz <- (delta.auc[ctr]) %*% rms::matinv(LSL) %*% delta.auc[ctr]
            ## tmpz <- (delta.auc[ctr]) %*% (1/LSL) %*% delta.auc[ctr]
            model[ctr] <- modelnames[j]
            reference[ctr] <- modelnames[i]
            ctr <- ctr + 1
        }
    }
    deltaAUC <- data.table(model,reference,delta.auc=as.vector(delta.auc),se.auc)
    deltaAUC[,lower:=delta.auc-Qnorm*se.auc]
    deltaAUC[,upper:=delta.auc+Qnorm*se.auc]
    deltaAUC[,p:=2*pnorm(abs(delta.auc)/se.auc,lower.tail=FALSE)]
    names(auc) <- 1:nauc
    auc <- data.table(auc, sqrt(diag(S)))
    setnames(auc,c("AUC","se.AUC"))
    auc[,model:=colnames(risk)]
    auc[,lower.AUC:=pmax(0,AUC-qnorm(1-alpha/2)*se.AUC)]
    auc[,upper.AUC:=pmin(1,AUC+qnorm(1-alpha/2)*se.AUC)]
    setcolorder(auc,c("model","AUC","se.AUC","lower.AUC","upper.AUC"))
    list(auc = auc, difference = deltaAUC)
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

AUC.binary <- function(DT,breaks=NULL,test,alpha,N,NT,NF,dolist,ROC,...){
    model=risk=ReSpOnSe=FPR=TPR=NULL
    data.table::setorder(DT,model)
    ## data.table::setorder(DT,model,risk)
    if (is.factor(DT[["risk"]])){
        score <- DT[,auRoc.factor(risk,ReSpOnSe,ROC=ROC),by=list(model)]
    }
    else{
        score <- DT[,auRoc.numeric(risk,ReSpOnSe,breaks=NULL,ROC=ROC),by=list(model)]
    }
    if (ROC==FALSE){
        setnames(score,"V1","AUC")
        out <- list(score=score)
    } else{
        AUC <- score[,list(AUC=0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))),by=list(model)]
        ROC <- score
        out <- list(score=AUC,ROC=ROC)
    }
    if (test==TRUE){
        xRisk <- data.table::dcast.data.table(DT,ID~model,value.var="risk")[,-1,with=FALSE]
        delong.res <- delongtest(risk=xRisk,score=out$score,dolist=dolist,response=DT[model==1,ReSpOnSe],cause="1",alpha=alpha)
        test.AUC <- delong.res$difference
        out$score <- delong.res$auc
        c(out,list(test=test.AUC))
    }else{
        out
    }
}

Brier.survival <- function(DT,MC,test,alpha,N,NT,NF,dolist,...){
    Yt=time=times=Residuals=risk=ipcwResiduals=WTi=Wt=status=setorder=model=IC.Brier=data.table=sd=lower.Brier=qnorm=se.Brier=upper.Brier=NULL
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
    if (test==TRUE){
        data.table::setorder(DT,model,times,time,-status)
        DT[,IC.Brier:=getInfluenceCurve.Brier(t=times[1],time=time,Yt=Yt,ipcwResiduals=ipcwResiduals,MC=MC),by=list(model,times)]
        score <- DT[,data.table(Brier=sum(ipcwResiduals)/N,se.Brier=sd(IC.Brier)/sqrt(N)),by=list(model,times)]
        score[,lower.Brier:=pmax(0,Brier-qnorm(1-alpha/2)*se.Brier)]
        score[,upper.Brier:=Brier + qnorm(1-alpha/2)*se.Brier]
        data.table::setkey(score,model,times)
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        test.Brier <- DT[,getComparisons(data.table(x=Brier,IC=IC.Brier,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist),by=list(times)]
        Brier <- list(score=score,test=test.Brier)
    }else{
         Brier <- list(score=DT[,data.table(Brier=sum(ipcwResiduals)/N),by=list(model,times)])
     }
}

Brier.competing.risks <- function(DT,MC,test,alpha,N,NT,NF,dolist,cause,...){
    Yt=time=times=event=Residuals=risk=ipcwResiduals=WTi=Wt=status=setorder=model=IC.Brier=data.table=sd=lower.Brier=qnorm=se.Brier=upper.Brier=NULL
    ## compute 0/1 outcome:
    DT[,Yt:=1*(time<=times & event==cause)]
    ## compute residuals
    DT[,Residuals:=(Yt-risk)^2]
    ## apply weights 
    DT[,ipcwResiduals:=Residuals/WTi]
    ## browser()
    DT[time>times,ipcwResiduals:=Residuals/Wt]
    ## deal with censored observations before t
    DT[time<=times & status==0,ipcwResiduals:=0]
    ## DT[,c("Yt","time","WTi","Wt","status","Residuals"):=NULL]
    if (test==TRUE){
        data.table::setorder(DT,model,times,time,-status)
        DT[,IC.Brier:=getInfluenceCurve.Brier(t=times[1],time=time,Yt=Yt,ipcwResiduals=ipcwResiduals,MC=MC),by=list(model,times)]
        score <- DT[,data.table(Brier=sum(ipcwResiduals)/N,se.Brier=sd(IC.Brier)/sqrt(N)),by=list(model,times)]
        score[,lower.Brier:=pmax(0,Brier-qnorm(1-alpha/2)*se.Brier)]
        score[,upper.Brier:=Brier + qnorm(1-alpha/2)*se.Brier]
        data.table::setkey(score,model,times)
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        test.Brier <- DT[,getComparisons(data.table(x=Brier,IC=IC.Brier,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist),by=list(times)]
        Brier <- list(score=score,test=test.Brier)
    }else{
         Brier <- list(score=DT[,data.table(Brier=sum(ipcwResiduals)/N),by=list(model,times)])
     }
    Brier
}

AireTrap <- function(FP,TP,N){
    N <- length(FP)
    sum((FP-c(0,FP[-N]))*((c(0,TP[-N])+TP)/2))
}

AUC.survival <- function(DT,MC,test,alpha,N,NT,NF,dolist,ROC,...){
    model=times=risk=Cases=time=status=Controls=TPR=FPR=WTi=Wt=ipcwControls=ipcwCases=IC.AUC=lower.AUC=se.AUC=upper.AUC=AUC=NULL
    cause <- 1
    ## assign Weights before ordering
    DT[,ipcwControls:=1/(Wt*N)]
    DT[,ipcwCases:=1/(WTi*N)]
    ## order data
    data.table::setorder(DT,model,times,-risk)
    ## identify cases and controls
    DT[,Cases:=(time <= times &  status==cause)]
    DT[,Controls:=(time > times)] 
    ## prepare Weights
    DT[Cases==0,ipcwCases:=0]
    DT[Controls==0,ipcwControls:=0]
    ## compute denominator
    DT[,TPR:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    DT[,FPR:=(cumsum(ipcwControls))/(sum(ipcwControls)),by=list(model,times)]
    nodups <- DT[,c(!duplicated(risk)[-1],TRUE),by=list(model,times)]$V1
    if (ROC==TRUE) {
        output <- list(ROC=DT[nodups,c("model","times","risk","TPR","FPR"),with=FALSE])
    }else{
         output <- NULL
     }
    score <- DT[nodups,list(AUC=AireTrap(FPR,TPR)),by=list(model,times)]
    data.table::setkey(score,model,times)
    if (test==TRUE){
        ## compute influence function
        data.table::setorder(DT,model,times,time,-status)
        DT[,IC.AUC:=getInfluenceCurve.AUC.survival(t=times[1],n=N,time=time,risk=risk,Cases=Cases,Controls=Controls,ipcwControls=ipcwControls,ipcwCases=ipcwCases,MC=MC), by=list(model,times)]
        se.score <- DT[,list(se.AUC=sd(IC.AUC)/sqrt(N)),by=list(model,times)]
        data.table::setkey(se.score,model,times)
        score <- score[se.score]
        score[,lower.AUC:=pmax(0,AUC-qnorm(1-alpha/2)*se.AUC)]
        score[,upper.AUC:=pmin(1,AUC+qnorm(1-alpha/2)*se.AUC)]
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        test.AUC <- DT[,getComparisons(data.table(x=AUC,IC=IC.AUC,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist),by=list(times)]
        output <- c(list(score=score,test=test.AUC),output)
    }else{
         output <- c(list(score=score),output)
     }
    output
}


AUC.competing.risks <- function(DT,MC,test,alpha,N,NT,NF,dolist,cause,ROC,...){
    model=times=risk=Cases=time=status=event=Controls1=Controls2=TPR=FPR=WTi=Wt=ipcwControls1=ipcwControls2=ipcwCases=IC.AUC=lower.AUC=se.AUC=upper.AUC=AUC=NULL
    ## assign Weights before ordering
    DT[,ipcwControls1:=1/(Wt*N)]
    DT[,ipcwControls2:=1/(WTi*N)]
    DT[,ipcwCases:=1/(WTi*N)]
    DT[,ipcwControls2:=1/(WTi*N)]
    ## order data
    data.table::setorder(DT,model,times,-risk)
    ## identify cases and controls
    DT[,Cases:=(time <=times &  event==cause)]
    DT[,Controls1:=(time > times)] 
    DT[,Controls2:=(time <=times &  event!=cause & status !=0)]
    ## prepare Weights
    DT[Cases==0,ipcwCases:=0]
    DT[Controls1==0,ipcwControls1:=0]
    DT[Controls2==0,ipcwControls2:=0]
    ## compute denominator
    ## ROC <- DT[,list(TPR=c(0,cumsum(ipcwCases)),FPR=c(0,cumsum(ipcwControls1)+cumsum(ipcwControls2))),by=list(model,times)]
    DT[,TPR:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    DT[,FPR:=(cumsum(ipcwControls1)+cumsum(ipcwControls2))/(sum(ipcwControls2)+sum(ipcwControls1)),by=list(model,times)]
    nodups <- DT[,c(!duplicated(risk)[-1],TRUE),by=list(model,times)]$V1
    if (ROC==TRUE) {
        output <- list(ROC=DT[nodups,c("model","times","risk","TPR","FPR"),with=FALSE])
    }else{
         output <- NULL
     }
    score <- DT[nodups,list(AUC=AireTrap(FPR,TPR)),by=list(model,times)]
    data.table::setkey(score,model,times)
    if (test==TRUE){
        ## compute influence function
        data.table::setorder(DT,model,times,time,-status)
        DT[,IC.AUC:=getInfluenceCurve.AUC.competing.risks(t=times[1],n=N,time=time,risk=risk,ipcwControls1=ipcwControls1,ipcwControls2=ipcwControls2,ipcwCases=ipcwCases,Cases=Cases,Controls1=Controls1,Controls2=Controls2,MC=MC), by=list(model,times)]
        se.score <- DT[,list(se.AUC=sd(IC.AUC)/sqrt(N)),by=list(model,times)]
        data.table::setkey(se.score,model,times)
        score <- score[se.score]
        score[,lower.AUC:=pmax(0,AUC-qnorm(1-alpha/2)*se.AUC)]
        score[,upper.AUC:=pmin(1,AUC+qnorm(1-alpha/2)*se.AUC)]
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        test.AUC <- DT[,getComparisons(data.table(x=AUC,IC=IC.AUC,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist),by=list(times)]
        output <- c(list(score=score,test=test.AUC),output)
    }else{
         output <- c(list(score=score),output)
     }
    ## browser(skipCalls=1)
    output
}



