### getPerformanceData.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Feb 27 2022 (09:12) 
## Version: 
## Last-Updated: Jun  5 2024 (18:02) 
##           By: Thomas Alexander Gerds
##     Update #: 66
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getPerformanceData <- function(testdata,
                               testweights,
                               traindata=NULL,
                               trainseed=NULL,
                               response.type,
                               neworder,
                               debug,
                               times,
                               cause,
                               levs,
                               labels,
                               predictRisk.args,
                               nullobject,
                               cens.type,
                               object,
                               object.classes,
                               NT,
                               verbose){
    riskRegression_ID=model = risk = NULL
    # inherit everything else from parent frame: object, nullobject, NF, NT, times, cause, response.type, etc.
    Brier=IPA=IBS=NULL
    looping <- length(traindata)>0
    N <- as.numeric(NROW(testdata))
    # split data vertically into response and predictors X
    rr_vars <- grep("^riskRegression_",names(testdata))
    testresponse <- testdata[,rr_vars,with=FALSE]
    data.table::setkey(testresponse,riskRegression_ID)
    X <- testdata[,-rr_vars,with=FALSE]
    if (debug) message("\nExtracted test set and prepared output object")
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
        rr_vars <- grep("^riskRegression_",names(testdata))
        trainX <- traindata[,-rr_vars,with=FALSE]
    }
    pred <- data.table::rbindlist(lapply(levs, function(f){
        ## add object specific arguments to predictRisk methods
        if (f[1]>0 && (length(extra.args <- unlist(lapply(object.classes[[f]],function(cc){predictRisk.args[[cc]]})))>0)){
            args <- c(args,extra.args)
        }
        ## predictions given as numeric values
        if (f[1]!=0 && any(c("integer","factor","numeric","matrix") %in% object.classes[[f]])){
            ## sort predictions by riskRegression_ID
            if (!is.null(dim(object[[f]]))) {## input matrix
                if (response.type=="binary"){
                    p <- do.call("predictRisk", c(list(object=c(object[[f]])),args))[neworder]
                }else{
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
                if (inherits(x=trained.model,what="try-error")){
                    if (verbose>1)message(paste0("Failed to train the following model:"))
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
            if (f[1]==0 && (response.type[1]!="binary")) {## glm predicts the same value for all subjects
                p <- rep(p,rep(N,NT))
            }
        }
        if (response.type%in%c("survival","competing.risks")){
            out <- data.table(riskRegression_ID=testdata[["riskRegression_ID"]],model=f,risk=p,times=rep(times,rep(N,NT)))
            byvars <- c("model","times")
            data.table::setkey(out,model,times,"riskRegression_ID")
            out
        } else {
            out <- data.table(riskRegression_ID=testdata[["riskRegression_ID"]],model=f,risk=p)
            byvars <- c("model")
            setkey(out,model,riskRegression_ID)
            out
        }
    }))
    if (any(is.na(pred$risk))) {
        if (verbose>1)message("Table of missing values in predicted risks:")
        pred[,model:=factor(model,levels=levs,labels)]
        if (response.type[1] == "binary"){
            print(pred[is.na(risk),data.table::data.table("sum(NA)" = .N),by = list(model)])
            stop("Missing values in predicted risk detected.")
        } else
            print(pred[is.na(risk),data.table::data.table("sum(NA)" = .N),by = list(model,times)])
        stop("Missing values in predicted risk detected.")
    }
    if (debug) message("\nTrained the model(s) and extracted the predictions")
    # }}}
    # {{{ merge with Weights (IPCW inner loop)
    if (response.type %in% c("survival","competing.risks")){
        if (cens.type=="rightCensored"){
            Weights <- testweights
            ## add subject specific weights
            set(testresponse,j="WTi",value=Weights$IPCW.subject.times)
        } else {
            if (cens.type=="uncensored"){
                Weights <- list(IPCW.times=rep(1,NT),IPCW.subject.times=matrix(1,ncol=NT,nrow=N))
                Weights$method <- "marginal"
                set(testresponse,j="WTi",value=1)
            } else{
                stop("Cannot handle this type of censoring.")
            }
        }
        ## add time point specific weights
        if (Weights$method=="marginal"){
            Wt <- data.table(times = times,Wt=c(Weights$IPCW.times))
            pred <- Wt[pred,,on="times"]
            data.table::setkey(pred,model,times,riskRegression_ID)
        }else{ # here as many weights as there are subjects at each element of times
            Wt <- rbindlist(lapply(1:length(times),function(s){
                data.table(riskRegression_ID = testresponse$riskRegression_ID,
                           times=rep(times[[s]],nrow(Weights$IPCW.times)),
                           Wt=Weights$IPCW.times[,s])
            }))
            pred <- Wt[pred,,on=c("riskRegression_ID","times")]
        }
        data.table::setkey(pred,model,times,riskRegression_ID)
        if (debug) message("merged the weights with input for performance metrics")
    } else {
        ## if (response.type=="binary")
        Weights <- NULL
    }
    if (debug) message("added weights to predictions")
    # }}}
    # {{{ merge with response
    DT=merge(testresponse,pred,by="riskRegression_ID")
    data.table::setkey(DT,riskRegression_ID)
    DT
}

#----------------------------------------------------------------------
### getPerformanceData.R ends here
