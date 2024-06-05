### crossvalPerf.loob.AUC.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2024 (09:20) 
## Version: 
## Last-Updated: Jun  5 2024 (17:57) 
##           By: Thomas Alexander Gerds
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
crossvalPerf.loob.AUC <- function(times,
                                  mlevs,
                                  se.fit,
                                  NT,
                                  response.type,
                                  response.names,
                                  cens.type,
                                  Weights,
                                  split.method,
                                  N,
                                  B,
                                  DT.B,
                                  dolist,
                                  alpha,
                                  byvars,
                                  mlabels,
                                  conservative,
                                  cens.model,
                                  cause,
                                  # crossvalPerf.loob.Brier has more arguments
                                  ...) {
    Response <- DT.B[,c(response.names,"riskRegression_ID"),with = FALSE][DT.B[,.I[1],by = "riskRegression_ID"]$V1]
    # initializing output
    bfold <- fold <- AUC <- riskRegression_event <- model <- b <- risk <- casecontrol <- IF.AUC <- IF.AUC0 <- se <- IF.AUC.conservative <- se.conservative <- lower <- upper <- NF <- reference  <-  riskRegression_event <- riskRegression_time <- riskRegression_status0 <- riskRegression_ID <- NULL
    if (response.type=="binary") {
        auc.loob <- data.table(expand.grid(times=0,model=mlevs)) #add times to auc.loob; now we can write less code for the same thing!
        times <- 0
    }
    auc.loob <- data.table::data.table(expand.grid(times=times,model=mlevs))
    auc.loob[,AUC:=as.numeric(NA)]
    ## for each pair of individuals sum the concordance of the bootstraps where *both* individuals are out-of-bag
    ## divide by number of times the pair is out-of-bag later
    aucDT <- NULL
    redoSplitIndex <- function(x, N){
        res <- rep(TRUE, N)
        res[unique(x)] <- FALSE
        res
    }
    ind.mat <- do.call("cbind",(lapply(1:B,split.method$index)))
    if (split.method$internal.name == "LeaveOneOutBoot"){
        split.index <- apply(ind.mat,2,function(x) redoSplitIndex(x,N))
    }
    else {
        setorder(DT.B, b, riskRegression_ID)
        DT.B$fold <- rep(c(ind.mat), each=NT*length(mlevs))
        DT.B[, bfold:=as.numeric(interaction(b,fold))] # construct new function
        setorder(DT.B, bfold, riskRegression_ID)
        k <- split.method$k
        split.index <- do.call("cbind",lapply(1:k, function(i) ind.mat==i))
    }
    cause <- as.numeric(cause)
    rm(ind.mat)
    # first index for c++ function
    if (response.type%in%c("survival","competing.risks"))
        first_hits <- sindex(Response[["riskRegression_time"]],times)-1
    for (s in 1:NT){
        if (response.type=="binary"){
            t <- 0
            ## the following indices have to be logical!!!
            cases.index <- Response[,riskRegression_event==1]
            controls.index <- !cases.index
        }else{
            t <- times[s]
            if (response.type == "survival"){
                Response[,riskRegression_status0 := riskRegression_status]
            }
            else { # competing risks
                oldEvent <- Response$riskRegression_event
                causeID <- oldEvent == cause
                cause1ID <- oldEvent == 1
                oldEvent[causeID] <- 1
                oldEvent[cause1ID] <- cause
                Response[,riskRegression_status0 := riskRegression_status * oldEvent]
                # need a status variable with value
                # 0 = censored
                # 1 = cause of interest
                # 2 = other causes
            }
            cases.index <- Response[,riskRegression_time<=t & riskRegression_status0==1]
            controls.index1 <- Response[,riskRegression_time>t]
            controls.index2 <-  Response[,riskRegression_status0==2 & riskRegression_time <= t]
            controls.index <- controls.index1 | controls.index2
        }
        # censoring weights
        if (cens.type=="rightCensored"){ 
            if (Weights$method == "marginal"){
                Wt <- Weights$IPCW.times[s]
            }else{
                Wt <- Weights$IPCW.times[,s]
            }
            weights <- as.numeric(1/Wt * controls.index1  + 1/Weights$IPCW.subject.times * (cases.index | controls.index2))
        }else{ ## uncensored
            weights <- rep(1,N)
        }
        w.cases <- weights[cases.index]
        w.controls <- weights[controls.index]
        n.cases <- sum(cases.index)
        n.controls <- sum(controls.index)
        riskRegression_ID.case <- Response[cases.index,]$riskRegression_ID
        riskRegression_ID.controls <- Response[controls.index,]$riskRegression_ID
        muCase <- sum(w.cases)
        muControls <- sum(w.controls)
        Phi <- (1/N^2)*muCase*muControls
        for (mod in mlevs){
            if (mod == 0){
                aucLPO <- 0.5
                auc.loob[times==t&model==mod,AUC:=aucLPO]
                DT.B <- DT.B[model != 0]
            }
            else {
                if (response.type == "binary"){ # no reason for the below trick with k fold, as level one data set is approximately N*B
                    DTrisk <- DT.B[model == mod][["risk"]]
                }
                else {
                    DTrisk <- DT.B[model == mod & times == t][["risk"]]
                }
                if (split.method$internal.name == "crossval"){
                    n.col <- k*B
                } else {
                    n.col <- B
                }
                risk.mat <- matrix(NA,nrow=N,ncol=n.col)
                risk.mat[split.index] <- DTrisk
                res <- aucLoobFun(riskRegression_ID.case,
                                  riskRegression_ID.controls,
                                  risk.mat,
                                  split.index,
                                  weights)
                ic0Case <- res[["ic0Case"]]
                ic0Control <- res[["ic0Control"]]
                warn <- res[["warn"]]
                nu1tauPm <- (1/N^2)*sum(ic0Case)
                ## Leave-one-pair-out bootstrap estimate of AUC
                aucLPO <- nu1tauPm*(1/Phi)
                # following section should be cleaned up(!!!)
                auc.loob[times==t&model==mod,AUC:=aucLPO]
            }
            if (se.fit==1L){ ## should clean this up when cv has been fixed !
                if (mod == 0){
                    aucDT <- data.table::data.table(model=mod,times=t,IF.AUC=rep(0,N))
                }
                else {
                    IF.AUC0 <- rep(0,N)
                    IF.AUC0[cases.index] <- 1/(Phi*N)*ic0Case
                    IF.AUC0[controls.index] <- 1/(Phi*N)*ic0Control
                    if (!response.type == "binary" && !cens.type == "uncensored" && !conservative[[1]]){
                        IFcalculationList <- list(ic0Case = ic0Case,
                                                  ic0Control = ic0Control,
                                                  weights = weights, 
                                                  muCase = muCase,
                                                  muControls = muControls, 
                                                  nu = nu1tauPm,
                                                  firsthit = first_hits[s],
                                                  cases = cases.index,
                                                  controls =controls.index,
                                                  controls1 = controls.index1,
                                                  controls2 = controls.index2)
                        icPart <- getInfluenceFunction.AUC.censoring.term(time = Response[["riskRegression_time"]],
                                                                          event = Response[["riskRegression_status0"]],
                                                                          t= t,
                                                                          IFcalculationList = IFcalculationList,
                                                                          MC = Weights$IC, 
                                                                          cens.model = cens.model, 
                                                                          Wt = Weights$IPCW.times[s], 
                                                                          auc = aucLPO, 
                                                                          nth.times = s)
                    }
                    else {
                        icPart <- 0
                    }
                    icPhi1 <- weights*(aucLPO/Phi)*(cases.index*(1/N)*muControls+controls.index*(1/N)*muCase)
                    this.aucDT <- data.table::data.table(model=mod,times=t,IF.AUC=IF.AUC0+icPart-icPhi1)
                    aucDT <- rbindlist(list(aucDT,this.aucDT),use.names=TRUE,fill=TRUE)
                }
                auc.loob[model==mod&times==t,se:= sd(aucDT[model==mod&times==t,IF.AUC])/sqrt(N)]
            }
        }
    }
    if (warn){
        warning("Some pairs of subjects are never out of bag at the same time. \n You should increase the number of bootstrap replications (argument 'B') ")
    }
    if (response.type == "binary"){
        # remove times again
        auc.loob[,times:=NULL] 
        if (se.fit){
            aucDT[,times:=NULL]
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
        setnames(contrasts.AUC,"delta","delta.AUC")
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


######################################################################
### crossvalPerf.loob.AUC.R ends here
