### crossvalPerf.loob.Brier.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2024 (09:16) 
## Version: 
## Last-Updated: feb 14 2026 (07:55) 
##           By: Thomas Alexander Gerds
##     Update #: 70
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
crossvalPerf.loob.Brier <- function(times,
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
                                    cause){
    riskRegression_status <- riskRegression_time <- residuals <- risk <- WTi <- riskRegression_event <- riskRegression_event <- Brier <- IC0 <- nth.times <- IF.Brier <- lower <- se <- upper <- model <- NF <- IPCW <- reference <- riskRegression_status0 <- IBS <- Wt <- .I <- response <- NULL
    ## sum across bootstrap samples where subject i is out of bag
    if (cens.type=="rightCensored"){
        if (response.type=="survival"){
            ## event of interest before times
            DT.B[riskRegression_time<=times & riskRegression_status==1,residuals:=(1-risk)^2/WTi]
        }
        else{ ## competing risks
            ## event of interest before times
            DT.B[riskRegression_time<=times & riskRegression_status==1 & riskRegression_event==cause,residuals:=(1-risk)^2/WTi]
            ## competing event before times
            DT.B[riskRegression_time<=times & riskRegression_status==1 &riskRegression_event!=cause,residuals:=(0-risk)^2/WTi]
        }
        ## right censored before times
        DT.B[riskRegression_time<=times & riskRegression_status==0,residuals:=0]
        ## no event at times
        DT.B[riskRegression_time>times,residuals:=(risk)^2/Wt]
    }else{
        DT.B[,residuals:=(riskRegression_event-risk)^2]
    }
    ## for each individual sum the residuals of the bootstraps where this individual is out-of-bag
    ## divide by number of times out-off-bag later
    ## keep response and weights for re-merging later
    if (response.type%in%c("survival","competing.risks")){
        DT.info <- DT.B[,c(response.names,"riskRegression_ID","WTi","Wt"),with = FALSE][DT.B[,.I[1],by = "riskRegression_ID"]$V1]
    }else{
        DT.info <- DT.B[,c(response.names,"riskRegression_ID"),with = FALSE][DT.B[,.I[1],by = "riskRegression_ID"]$V1]
    }
    DT.B <- DT.B[,data.table::data.table(risk=mean(risk),residuals=mean(residuals),n.oob = .N), by=c(byvars,"riskRegression_ID")]
    never.oob <- N-length(unique(DT.B$riskRegression_ID))
    ## n.oob is the count of how many times subject i is out of bag
    ## the order of n.oob matches the order of riskRegression_ID
    ## within groups defined by times and model (byvars)
    data.table::setkeyv(DT.B,c(byvars,"riskRegression_ID"))
    ## DT.B[,residuals:=residuals/Ib]
    ## leave-one-out bootstrap estimate
    DT.B[,Brier:=mean(residuals),by=byvars]
    ## standard error via influence function
    if (se.fit==1L){
        ## influence function when censoring model is known or data are uncensored
        DT.B[,IC0:=residuals-Brier]
        if (never.oob>0) {
            conservative <- TRUE
        }
        if (cens.type[1]=="rightCensored" && !conservative){
            # merge with response, riskRegression_ID and weights
            DT.B <- DT.info[DT.B,,on="riskRegression_ID"]
            DT.B[,nth.times:=as.numeric(factor(times))]
            data.table::setkeyv(DT.B,c(byvars,"riskRegression_ID"))
            if (cens.type=="uncensored"){
                DT.B[,IF.Brier:= residuals]
                score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,se=sd(residuals)/sqrt(N)),by=byvars]
            }else{
                #for small values of B, there is the problem that
                #some individuals might be zero times out of the bag
                #this means that DT.B, which should have a number of rows
                # that is a multiple of the
                #amount of observations in the data, does not fulfill this.
                #the calculations in getInfluenceCurve.Brier cannot accomodate this (for now).
                if (response.type == "survival"){
                    DT.B[,riskRegression_status0:=riskRegression_status]
                }
                else {
                    DT.B[,riskRegression_status0:=riskRegression_status*riskRegression_event]
                }
                DT.B[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                                        time=riskRegression_time,
                                                        IC0 = IC0,
                                                        residuals=residuals,
                                                        IC.G=Weights$IC,
                                                        cens.model=cens.model,
                                                        conservative = conservative,
                                                        nth.times=nth.times[1],
                                                        event = riskRegression_status0),by=list(model,times)]
                score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,
                                               se=sd(IF.Brier)/sqrt(N)),by=byvars]
            }
        }else{
            ## either binary or uncensored
            score.loob <- DT.B[,data.table(Brier=sum(residuals)/N,se=sd(IC0)/sqrt(N)),
                               by=byvars]
            setnames(DT.B,"IC0","IF.Brier")
        }
        score.loob[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
        score.loob[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
    } else{
        score.loob <- DT.B[,list(Brier=sum(residuals)/N), by=byvars]
    }
    if (ipa==TRUE){
        if (se.fit == TRUE){
            ipa_vars <- c("riskRegression_ID",byvars,"IF.Brier")
        }else{
            ipa_vars <- c("riskRegression_ID",byvars)
        }
        IPA <- IPACompute(brier.results = list(score = score.loob,iid.decomp = DT.B[,ipa_vars,with = FALSE]),
                          response.type = response.type,
                          se.fit = se.fit,
                          alpha = alpha)
        if (response.type=="binary")
            score.loob[,IPA:=1-Brier/Brier[model==0]]
        else
            score.loob[,IPA:=1-Brier/Brier[model==0],by=times]
    }
    ## summary should be after metrics because IBS and IPA/R^2 depends on Brier score
    if (ibs){
        Dint <- function(x,y,range,na.omit=FALSE){
            if (is.null(range)) range=c(x[1],x[length(x)])
            ##   integrate a step function f with
            ##   values y=f(x) between range[1] and range[2]
            start <- max(range[1],min(x))
            Stop <- min(range[2],max(x))
            if ((Stop-start)<=0)
                return(0)
            else{
                Y=y[x>=start & x<Stop]
                X=x[x>=start & x<Stop]
                if (na.omit){
                    X=X[!is.na(Y)]
                    Y=Y[!is.na(Y)]
                } else if (any(is.na(Y))|| any(is.na(X))){
                    return(NA)
                }
                return(1/(Stop-start) * sum(Y*diff(c(X,Stop))))
            }
        }
        if (response.type!="binary"){
            score.loob[,IBS:=sapply(times,function(t){
                Dint(x=c(0,times),y=c(0,Brier),range=c(0,t))
            }),by=c("model")]
        }
    }
    
    data.table::setkeyv(score.loob,byvars)
    if (length(dolist)>0L){
        if ("times"%in% byvars){by_statement <- "times"} else{by_statement <- NULL}
        if (se.fit){
            contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,model=model,IF = IF.Brier),NF=NF,N=N,alpha=alpha,dolist=dolist,se.fit=TRUE),by=by_statement]
        }else{
            contrasts.Brier <- DT.B[,getComparisons(data.table(x=Brier,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist,se.fit=FALSE),by=by_statement]
        }
        setnames(contrasts.Brier,"delta","delta.Brier")
        output <- list(score=score.loob,contrasts=contrasts.Brier)
    } else{
        output <- list(score=score.loob)
    }
    if (keep.residuals) {
        if (all(c("Wt","WTi")%in%names(DT.B))){
            DT.B[,IPCW:=1/WTi]
            DT.B[riskRegression_time>=times,IPCW:=1/Wt]
            DT.B[riskRegression_time<times & riskRegression_status==0,IPCW:=0]
            output <- c(output,list(residuals=DT.B[,c(response.names,"model","times","risk","residuals","IPCW"),with=FALSE]))
        }else{
            output <- c(output,list(residuals=DT.B[,c(response.names,"model","times","risk","residuals"),with=FALSE]))
        }
    }
    if (!is.null(output$score)){
        data.table::setkeyv(output$score,byvars)
    }
    ## set model and reference in model comparison results
    if (!is.null(output$contrasts)>0){
        data.table::setkeyv(output$score,byvars)
    }
    if (ipa==TRUE){
        output[["IPA"]] <- IPA
    }
    attr(output,"n_subjects_never_oob") = never.oob
    return(output)
}

######################################################################
### crossvalPerf.loob.Brier.R ends here
