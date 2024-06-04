### AUC.competing.risks.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:06)
## Version:
## Last-Updated: Jun  4 2024 (14:35) 
##           By: Thomas Alexander Gerds
##     Update #: 36
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

AUC.competing.risks <- function(DT,
                                breaks=NULL,
                                MC,
                                se.fit,
                                conservative,
                                cens.model,
                                keep.vcov=FALSE,
                                keep.iid=FALSE,
                                multi.split.test,
                                alpha,
                                N,
                                NT,
                                NF,
                                dolist,
                                cause,
                                states,
                                ROC,
                                IC.data,
                                cutpoints,
                                ...){
    riskRegression_ID=model=times=risk=Cases=riskRegression_time=riskRegression_status=riskRegression_event=Controls1=Controls2=TPR=FPR=WTi=Wt=ipcwControls1=ipcwControls2=ipcwCases=IF.AUC=lower=se=upper=AUC=nth.times=NULL
    aucDT <- DT[model>0]
    dolist <- dolist[sapply(dolist,function(do){match("0",do,nomatch=0L)})==0]
    ## assign Weights before ordering
    aucDT[,ipcwControls1:=1/(Wt*N)]
    aucDT[,ipcwControls2:=1/(WTi*N)]
    aucDT[,ipcwCases:=1/(WTi*N)]
    ## order data
    data.table::setorder(aucDT,model,times,-risk)
    ## identify cases and controls
    thecause <- match(cause,states,nomatch=0)
    if (length(thecause)==0) stop("Cannot identify cause of interest")
    aucDT[,Cases:=(riskRegression_time <=times &  riskRegression_event==thecause)]
    aucDT[,Controls1:=(riskRegression_time > times)]
    aucDT[,Controls2:=(riskRegression_time <=times &  riskRegression_event!=thecause & riskRegression_status !=0)]
    ## prepare Weights
    aucDT[Cases==0,ipcwCases:=0]
    aucDT[Controls1==0,ipcwControls1:=0]
    aucDT[Controls2==0,ipcwControls2:=0]
    ## compute denominator
    ## ROC <- aucDT[,list(TPR=c(0,cumsum(ipcwCases)),FPR=c(0,cumsum(ipcwControls1)+cumsum(ipcwControls2))),by=list(model,times)]
    aucDT[,TPR:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    aucDT[,FPR:=(cumsum(ipcwControls1)+cumsum(ipcwControls2))/(sum(ipcwControls2)+sum(ipcwControls1)),by=list(model,times)]
    nodups <- aucDT[,c(!duplicated(risk)[-1],TRUE),by=list(model,times)]$V1
    if (!is.null(cutpoints)){
        breaks <- sort(cutpoints,decreasing = TRUE)
        aucDT[,nth.times:=as.numeric(factor(times))]
        cutpoint.helper.fun <- function(FPR,
                                        TPR,
                                        risk,
                                        ipcwCases,
                                        ipcwControls1,
                                        ipcwControls2,
                                        N,
                                        riskRegression_time,
                                        times,
                                        riskRegression_event,
                                        cens.model,
                                        nth.times,
                                        conservative,
                                        IC.G,
                                        cutpoints,
                                        se.fit){
            den_TPR<-sum(ipcwCases) ## estimate the cumulative incidence via IPCW
            den_FPR<-sum(ipcwControls1+ipcwControls2)
            indeces <- sindex(risk,cutpoints,comp = "greater",TRUE)
            res <- list()
            # FIXME
            ordered <- order(riskRegression_time) ## can probably move this outside to improve computation time, for now keep it
            for (i in 1:length(cutpoints)){
                den_PPV <- sum(ipcwCases[risk > cutpoints[i]]+ipcwControls1[risk > cutpoints[i]] + ipcwControls2[risk > cutpoints[i]])
                den_NPV <- 1-den_PPV
                if (indeces[i] != 0){
                    TPRi <- TPR[indeces[i]]
                    FPRi <- FPR[indeces[i]]
                    if (se.fit){
                        IC0.TPR <- ipcwCases*N*((risk > cutpoints[i])-TPRi)/den_TPR
                        IC0.FPR <- (ipcwControls1+ipcwControls2)*N*((risk > cutpoints[i])-FPRi)/(1-den_TPR)
                        SE.TPR <- sd(getInfluenceCurve.Brier(t = times,
                                                             time = riskRegression_time[ordered],
                                                             IC0 = IC0.TPR[ordered],
                                                             residuals = IC0.TPR[ordered],
                                                             IC.G = IC.G,
                                                             cens.model = cens.model,
                                                             nth.times = nth.times,
                                                             conservative = conservative,
                                                             event = riskRegression_event[ordered]))/sqrt(N)
                        SE.FPR <- sd(getInfluenceCurve.Brier(t = times,
                                                             time = riskRegression_time[ordered],
                                                             IC0 = IC0.FPR[ordered],
                                                             residuals = IC0.FPR[ordered],
                                                             IC.G = IC.G,
                                                             cens.model = cens.model,
                                                             nth.times = nth.times,
                                                             conservative = conservative,
                                                             event = riskRegression_event[ordered]))/sqrt(N)
                    }
                    else {
                        SE.TPR <- SE.FPR <- NA
                    }
                }
                else {
                    TPRi <- FPRi <- 0
                    SE.TPR <- SE.FPR <- NA
                }
                if (den_PPV > 1e-10){
                    PPV <- (TPRi*den_TPR)/den_PPV
                    if (se.fit){
                        IC0.PPV <- (risk > cutpoints[i])/den_PPV*(((ipcwCases+ipcwControls2)*N)*(1*(riskRegression_event==1)-1*(riskRegression_event!=0)*PPV)-ipcwControls1*N*PPV) #OBS, check other causes, paul's implementation
                        SE.PPV <- sd(getInfluenceCurve.Brier(t = times,
                                                             time = riskRegression_time[ordered],
                                                             IC0 = IC0.PPV[ordered],
                                                             residuals = IC0.PPV[ordered],
                                                             IC.G = IC.G,
                                                             cens.model = cens.model,
                                                             nth.times = nth.times,
                                                             conservative = conservative,
                                                             event = riskRegression_event[ordered]))/sqrt(N)
                    }
                    else {
                        SE.PPV <- NA
                    }
                }
                else {
                    PPV <- NA
                }
                if (den_NPV > 1e-10){
                    NPV <- ((1-FPRi)*den_FPR)/den_NPV
                    if (se.fit){
                        IC0.NPV <- (risk <= cutpoints[i])/den_NPV*(((ipcwCases+ipcwControls2)*N)*(1*(riskRegression_event!=1 & riskRegression_event!=0)-1*(riskRegression_event!=0)*NPV)+ipcwControls1*N*(1-NPV)) #OBS, check other causes, paul's implementation
                        SE.NPV <- sd(getInfluenceCurve.Brier(t = times,
                                                             time = riskRegression_time[ordered],
                                                             IC0 = IC0.NPV[ordered],
                                                             residuals = IC0.NPV[ordered],
                                                             IC.G = IC.G,
                                                             cens.model = cens.model,
                                                             nth.times = nth.times,
                                                             conservative = conservative,
                                                             event = riskRegression_event[ordered]))/sqrt(N)
                    }
                    else {
                        SE.NPV <- NA
                    }
                }
                else {
                    NPV <- NA
                }
                res[[i]] <- data.table(risk = cutpoints[i],
                                       TPR=TPRi,
                                       SE.TPR=SE.TPR,
                                       FPR=FPRi,
                                       SE.FPR=SE.FPR,
                                       PPV=PPV,
                                       SE.PPV=SE.PPV,
                                       NPV=NPV,
                                       SE.NPV=SE.NPV)
            }
            do.call("rbind",res)
        }
        output <- list(res.cut=aucDT[, cutpoint.helper.fun(FPR = FPR,
                                                           TPR = TPR,
                                                           risk = risk,
                                                           ipcwCases = ipcwCases,
                                                           ipcwControls1 = ipcwControls1,
                                                           ipcwControls2 = ipcwControls2,
                                                           N = N,
                                                           riskRegression_time = riskRegression_time,
                                                           times = times[1],
                                                           riskRegression_event = riskRegression_status*riskRegression_event,
                                                           cens.model = cens.model,
                                                           nth.times = nth.times[1],
                                                           conservative = conservative,
                                                           IC.G = MC,
                                                           cutpoints = cutpoints,
                                                           se.fit = se.fit),by=list(model,times)])
    }
    else if (ROC) {
        if (is.null(breaks)){
            output <- list(ROC=aucDT[nodups,c("model","times","risk","TPR","FPR"),with=FALSE])
        }
        else {
            breaks <- sort(breaks,decreasing = TRUE)
            helper.fun <- function(FPR,TPR,risk, breaks){
                indeces <- sindex(risk,breaks,comp = "greater",FALSE)
                data.table(risk = breaks, TPR = c(rep(0,sum(indeces==0)),TPR[indeces[indeces!=0]]), FPR = c(rep(0,sum(indeces==0)),FPR[indeces[indeces!=0]]))
            }
            output <- list(ROC=aucDT[, helper.fun(FPR,TPR,risk,breaks=breaks),by=list(model,times)])
        }
    }else{
        output <- NULL
    }
    AireTrap <- function(FP,TP,N){
        N <- length(FP)
        sum((FP-c(0,FP[-N]))*((c(0,TP[-N])+TP)/2))
    }
    score <- aucDT[nodups,list(AUC=AireTrap(FPR,TPR)),by=list(model,times)]
    data.table::setkey(score,model,times)
    aucDT <- merge(score,aucDT,all=TRUE)
    if (se.fit[[1]]==1L || multi.split.test[[1]]==TRUE){
        aucDT[,nth.times:=as.numeric(factor(times))]
        
        ## compute influence function
        ## data.table::setorder(aucDT,model,times,time,-riskRegression_status)
        data.table::setorder(aucDT,model,times,riskRegression_ID)
        aucDT[,IF.AUC:=getInfluenceCurve.AUC(t = times[1],
                                             time = riskRegression_time,
                                             event = riskRegression_status*riskRegression_event,
                                             WTi = WTi,
                                             Wt = Wt,
                                             risk = risk,
                                             MC = MC,
                                             auc = AUC[1],
                                             nth.times = nth.times[1],
                                             conservative = conservative[[1]],
                                             cens.model = cens.model), by=list(model,times)]
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
        if (keep.vcov[[1]] == TRUE){
            output <- c(output,list(vcov=getVcov(aucDT,"IF.AUC",times=TRUE)))
        }
        if (keep.iid[[1]] == TRUE && se.fit[[1]] == TRUE) {
            output <- c(output,
                        list(iid.decomp = aucDT[,data.table::data.table(riskRegression_ID,model,cause,times,IF.AUC)]))
        }
    }
    ## add score to object
    output <- c(list(score=score),output)
    if (length(dolist)>0){
        if (se.fit[[1]]==TRUE || multi.split.test[[1]]==TRUE){
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

#----------------------------------------------------------------------
### AUC.competing.risks.R ends here
