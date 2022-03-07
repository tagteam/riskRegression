### getPerformanceData.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Feb 27 2022 (09:12) 
## Version: 
## Last-Updated: Mar  7 2022 (08:33) 
##           By: Thomas Alexander Gerds
##     Update #: 10
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
                               response.dim,
                               neworder,
                               debug,
                               times,
                               cause,
                               labels,
                               predictRisk.args,
                               nullobject,
                               cens.type,
                               object,
                               object.classes,
                               NT){
    ID=model=NULL
    # inherit everything else from parent frame: object, nullobject, NF, NT, times, cause, response.type, etc.
    Brier=IPA=IBS=NULL
    looping <- !is.null(traindata)
    ## if (!looping) b=0
    N <- as.numeric(NROW(testdata))
    # split data vertically into response and predictors X
    response <- testdata[,1:response.dim,with=FALSE]
    response[,ID:=testdata[["ID"]]]
    setkey(response,ID)
    X <- testdata[,-c(1:response.dim),with=FALSE]
    ## restore sanity
    setnames(X,sub("^protectedName.","",names(X)))
    ## if (debug) if (looping) message(paste0("Loop round: ",b))
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
        ## restore sanity
        setnames(trainX,sub("^protectedName.","",names(trainX)))
        ## trainX <- copy(traindata)
        trainX[,ID:=NULL]
    }
    pred <- data.table::rbindlist(lapply(labels, function(f){
        ## add object specific arguments to predictRisk methods
        if (f[1]>0 && (length(extra.args <- unlist(lapply(object.classes[[f]],function(cc){predictRisk.args[[cc]]})))>0)){
            args <- c(args,extra.args)
        }
        ## predictions given as numeric values
        if (f[1]!=0 && any(c("integer","factor","numeric","matrix") %in% object.classes[[f]])){
            ## sort predictions by ID
            if (!is.null(dim(object[[f]]))) {## input matrix
                if (response.type=="binary"){
                    p <- do.call("predictRisk", c(list(object=c(object[[f]])),args))[neworder]
                }else{
                    ## if(!is.null(include.times)){ ## remove columns at times beyond max time
                    ## p <- c(do.call("predictRisk",c(list(object=object[[f]][,include.times,drop=FALSE]),args))[neworder,])
                    ## } else{
                    p <- c(do.call("predictRisk",c(list(object=object[[f]]),args))[neworder,])
                    ## }
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
                # browser(skipCalls = 1)
                trained.model <- try(eval(model.f$call),silent=TRUE)
                if (inherits(x=trained.model,what="try-error")){
                    message(paste0("Failed to train the following model:"))
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
            set(response,j="WTi",value=Weights$IPCW.subject.times)
        } else {
            if (cens.type=="uncensored"){
                Weights <- list(IPCW.times=rep(1,NT),IPCW.subject.times=matrix(1,ncol=NT,nrow=N))
                Weights$method <- "marginal"
                set(response,j="WTi",value=1)
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


#----------------------------------------------------------------------
### getPerformanceData.R ends here
