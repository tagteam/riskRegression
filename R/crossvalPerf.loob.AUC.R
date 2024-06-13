### crossvalPerf.loob.AUC.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2024 (09:20) 
## Version: 
## Last-Updated: Jun 13 2024 (15:40) 
##           By: Thomas Alexander Gerds
##     Update #: 93
#----------------------------------------------------------------------
## 
### Commentary:
##
## For each pair of individuals sum the concordance of the bootstraps
## where *both* individuals are out-of-bag. Then divide by number of
## times the pair is out-of-bag to estimate the AUC. 
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
                                  N,
                                  B,
                                  DT.B,
                                  dolist,
                                  alpha,
                                  byvars,
                                  mlabels,
                                  ipa,
                                  ibs,
                                  keep.residuals,
                                  conservative,
                                  cens.model,
                                  cause) {
    bfold <- fold <- AUC <- riskRegression_event <- model <- b <- risk <- casecontrol <- IF.AUC <- IF.AUC0 <- se <- IF.AUC.conservative <- se.conservative <- lower <- upper <- NF <- reference  <-  riskRegression_event <- riskRegression_time <- riskRegression_status0 <- riskRegression_status <- riskRegression_ID <- .I <- NULL
    setkeyv(DT.B,c(byvars,"riskRegression_ID"))
    if (cens.type == "rightCensored"){
        Response <- DT.B[,c(response.names,"riskRegression_ID","WTi","Wt"),with = FALSE][DT.B[,.I[1],by = "riskRegression_ID"]$V1]
    }else{
        Response <- DT.B[,c(response.names,"riskRegression_ID"),with = FALSE][DT.B[,.I[1],by = "riskRegression_ID"]$V1]}
    setkey(Response,riskRegression_ID)
    N_has_oob <- length(unique(DT.B$riskRegression_ID))
    if (N_has_oob<N) conservative <- TRUE
    oob_warning <- FALSE
    # initializing output
    if (response.type=="binary") {
        auc.loob <- data.table(expand.grid(times=0,model=mlevs)) #add times to auc.loob; now we can write less code for the same thing!
        times <- 0
    }
    auc.loob <- data.table::data.table(expand.grid(times=times,model=mlevs))
    auc.loob[,AUC:=as.numeric(NA)]
    aucDT <- NULL
    # TRUE FALSE matrix indicating which ID's are out-of-bag
    split.index <- matrix(DT.B[,1:N%in%riskRegression_ID,keyby = "b"]$V1,ncol = B,byrow = FALSE)
    cause <- as.numeric(cause)
    # first index for c++ function
    if (response.type%in%c("survival","competing.risks")) {
        first_hits <- sindex(Response[["riskRegression_time"]],times)-1
    }
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
        if ((length(unique(controls.index)) > 1 && length(unique(controls.index)) > 1)) { # avoid trivial cases without cases or controls
            # censoring weights
            if (cens.type=="rightCensored"){
                weights <- as.numeric(1/Response$Wt * controls.index1  + 1/Response$WTi * (cases.index | controls.index2))
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
                    # split.index(i,b) is TRUE if subject i is out-of-bag in bootstrap b
                    # risk.mat(i,b) is the risk of this subject if it is out-of-bag, else any other value (we use NA)
                    # both have dimention N x B
                    all.subjects <- data.table(riskRegression_ID = 1:N)
                    if (response.type == "binary"){
                        risk.mat <- matrix(DT.B[model == mod,{
                            data.table(riskRegression_ID,risk)[all.subjects,on = "riskRegression_ID"]
                        },keyby = "b"]$risk,
                        ncol = B,
                        byrow = FALSE)
                    }else{
                        risk.mat <- matrix(DT.B[model == mod & times == t,{
                            data.table(riskRegression_ID,risk)[all.subjects,on = "riskRegression_ID"]
                        },keyby = "b"]$risk,
                        ncol = B,
                        byrow = FALSE)
                    }
                    res <- aucLoobFun(riskRegression_ID.case,
                                      riskRegression_ID.controls,
                                      risk.mat,
                                      split.index,
                                      weights)
                    ic0Case <- res[["ic0Case"]]
                    ic0Control <- res[["ic0Control"]]
                    oob_warning <- res[["warn"]]
                    if (oob_warning) conservative <- TRUE
                    nu1tauPm <- (1/N^2)*sum(ic0Case)
                    ## Leave-one-pair-out bootstrap estimate of AUC
                    aucLPO <- nu1tauPm*(1/Phi)
                    # following section should be cleaned up(!!!)
                    auc.loob[times==t&model==mod,AUC:=aucLPO]
                }
                if (se.fit==1L){ ## should clean this up when cv has been fixed !
                    if (mod == 0){
                        aucDT <- data.table::data.table(model=mod,times=t,IF.AUC=rep(0,N_has_oob))
                    }
                    else {
                        IF.AUC0 <- rep(0,N_has_oob) # force the correct length in cases where not all pairs are oob 
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
            if (length(aucDT) == 0){ # this happens when there are either no controls or no cases
                contrasts.AUC = NULL
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
                setnames(contrasts.AUC,"delta","delta.AUC")
                output <- c(output,list(contrasts=contrasts.AUC))
            }
        }
    }
    if (!is.null(output$score)){
        output$score[,model:=factor(model,levels=mlevs,mlabels)]
    }
    ## set model and reference in model comparison results
    if (!is.null(output$contrasts)){
        output$contrasts[,model:=factor(model,levels=mlevs,mlabels)]
        output$contrasts[,reference:=factor(reference,levels=mlevs,mlabels)]
    }
    if (oob_warning[1] == TRUE)
        attr(output,"subjects_never_oob") <- 1L
    else
        attr(output,"subjects_never_oob") <- 0L
    return(output)
}


######################################################################
### crossvalPerf.loob.AUC.R ends here
