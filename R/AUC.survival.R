### AUC.survival.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:06)
## Version:
## Last-Updated: Jun 25 2024 (09:47) 
##           By: Thomas Alexander Gerds
##     Update #: 52
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
AUC.survival <- function(DT,
                         breaks=NULL,
                         MC,
                         se.fit,
                         conservative,
                         cens.model,
                         keep.vcov=FALSE,
                         keep.iid=FALSE,
                         alpha,
                         N,
                         NT,
                         NF,
                         dolist,
                         ROC,
                         IC.data,
                         cutpoints,
                         ...){
    riskRegression_ID=model=times=risk=Cases=riskRegression_time=riskRegression_status=Controls=TPR=FPR=WTi=Wt=ipcwControls=ipcwCases=IF.AUC=lower=se=upper=AUC=nth.times=NULL
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
    aucDT[,Cases:=(riskRegression_time <= times &  riskRegression_status==cause)]
    aucDT[,Controls:=(riskRegression_time > times)]
    ## prepare Weights
    aucDT[Cases==0,ipcwCases:=0]
    aucDT[Controls==0,ipcwControls:=0]
    ## compute denominator
    aucDT[,TPR:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)] # technically sum_i I(M_i >= M_j) not M_i > M_j
    aucDT[,FPR:=(cumsum(ipcwControls))/(sum(ipcwControls)),by=list(model,times)]
    nodups <- aucDT[,c(!duplicated(risk)[-1],TRUE),by=list(model,times)]$V1
    AireTrap <- function(FP,TP){
        N <- length(FP)
        sum((FP-c(0,FP[-N]))*((c(0,TP[-N])+TP)/2))
    }
    score <- aucDT[nodups,list(AUC=AireTrap(FPR,TPR)),by=list(model,times)]
    if (!is.null(cutpoints)){
        breaks <- sort(cutpoints,decreasing = TRUE)
        aucDT[,nth.times:=as.numeric(factor(times))]
        cutpoint.helper.fun <- function(FPR,
                                        TPR,
                                        risk,
                                        ipcwCases,
                                        ipcwControls,
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
            den_FPR<-sum(ipcwControls)
            indeces <- sindex(risk,cutpoints,comp = "greater",TRUE)
            res <- list()
            # FIXME
            ordered <- order(riskRegression_time) ## can probably move this outside to improve computation time, for now keep it
            SE.TPR <- SE.FPR <- SE.PPV <- SE.NPV <- NA
            for (i in 1:length(cutpoints)){
                den_PPV <- sum(ipcwCases[risk > cutpoints[i]]+ipcwControls[risk > cutpoints[i]])
                den_NPV <- 1-den_PPV
                if (indeces[i] != 0){
                    TPRi <- TPR[indeces[i]]
                    FPRi <- FPR[indeces[i]]
                    if (se.fit){
                        IC0.TPR <- ipcwCases*N*((risk > cutpoints[i])-TPRi)/den_TPR
                        SE.TPR <- sd(getInfluenceCurve.Brier(t = times,
                                                             time = riskRegression_time[ordered],
                                                             IC0 = IC0.TPR[ordered],
                                                             residuals = IC0.TPR[ordered],
                                                             IC.G = IC.G,
                                                             cens.model = cens.model,
                                                             nth.times = nth.times,
                                                             conservative = conservative,
                                                             event = riskRegression_event[ordered]))/sqrt(N)
                        IC0.FPR <- (ipcwControls)*N*((risk > cutpoints[i])-FPRi)/(1-den_TPR)
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
                }
                else {
                    TPRi <- FPRi <- 0
                }
                if (den_PPV > 1e-10){
                    PPV <- (TPRi*den_TPR)/den_PPV
                    if (se.fit){
                        #OBS, check other causes, paul's implementation
                        IC0.PPV <- (risk > cutpoints[i])/den_PPV*(((ipcwCases)*N)*(1*(riskRegression_event==1)-1*(riskRegression_event!=0)*PPV)-ipcwControls*N*PPV)
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
                }
                else {
                    PPV <- NA
                }
                if (den_NPV > 1e-10){
                    NPV <- ((1-FPRi)*den_FPR)/den_NPV
                    if (se.fit){
                        IC0.NPV <- (risk <= cutpoints[i])/den_NPV*(((ipcwCases)*N)*(1*(riskRegression_event!=1 & riskRegression_event!=0)-1*(riskRegression_event!=0)*NPV)+ipcwControls*N*(1-NPV)) #OBS, check other causes, paul's implementation
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
        output <- list(cutpoints=aucDT[,
                                     cutpoint.helper.fun(FPR = FPR,
                                                         TPR = TPR,
                                                         risk = risk,
                                                         ipcwCases = ipcwCases,
                                                         ipcwControls = ipcwControls,
                                                         N = N,
                                                         riskRegression_time = riskRegression_time,
                                                         times = times[1],
                                                         riskRegression_event = riskRegression_status,
                                                         cens.model = cens.model,
                                                         nth.times = nth.times[1],
                                                         conservative = conservative,
                                                         IC.G = MC,
                                                         cutpoints = cutpoints,
                                                         se.fit = se.fit),by=list(model,times)])
    }
    else if (ROC==TRUE) {
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
    aucDT <- merge(score,aucDT,by = c("model","times"),all=TRUE)
    data.table::setkey(aucDT,model,times)
    if (se.fit[[1]]==1L){
        aucDT[,nth.times:=as.numeric(factor(times))]
        ## compute influence function
        ## data.table::setorder(aucDT,model,times,riskRegression_time,-riskRegression_status)
        data.table::setorder(aucDT,model,times,riskRegression_ID)
        aucDT[,IF.AUC:=getInfluenceCurve.AUC(t = times[1],
                                             time = riskRegression_time,
                                             event = riskRegression_status,
                                             WTi = WTi,
                                             Wt = Wt,
                                             risk = risk,
                                             MC = MC,
                                             auc = AUC[1],
                                             nth.times = nth.times[1],
                                             conservative = conservative[[1]],
                                             cens.model = cens.model), by=list(model,times)]
        se.score <- aucDT[,list(se=sd(IF.AUC)/sqrt(N)),by=list(model,times)]
        score <- score[se.score,,on = c("model","times")]
        data.table::setkey(score,model,times)
        if (se.fit==1L){
            score[,lower:=pmax(0,AUC-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,AUC+qnorm(1-alpha/2)*se)]
        }else{
            score[,se:=NULL]
        }
        # join with auc and (se, lower, upper) if se.fit
        aucDT <- score[aucDT,,on = c("model","times")]
        data.table::setkey(aucDT,model,times)
        if (keep.vcov[[1]] == TRUE){
            output <- c(output,list(vcov=getVcov(aucDT,"IF.AUC",times=TRUE)))
        }
        if (keep.iid[[1]] == TRUE && se.fit[[1]] == TRUE) {
            output <- c(output,
                        list(iid.decomp = aucDT[,data.table::data.table(riskRegression_ID,model,times,IF.AUC)]))
        }
        
    }
    ## add score to object
    output <- c(list(score=score),output)
    if (length(dolist)>0){
        if (se.fit[[1]]==TRUE){
            contrasts.AUC <- aucDT[,getComparisons(data.table(x=AUC,IF=IF.AUC,model=model),
                                                   NF=NF,
                                                   N=N,
                                                   alpha=alpha,
                                                   dolist=dolist,se.fit=se.fit),by=list(times)]
        }else{
            contrasts.AUC <- score[,getComparisons(data.table(x=AUC,model=model),
                                                   NF=NF,
                                                   N=N,
                                                   alpha=alpha,
                                                   dolist=dolist,
                                                   se.fit=FALSE),by=list(times)]
        }
        setnames(contrasts.AUC,"delta","delta.AUC")
        output <- c(list(score=score,contrasts=contrasts.AUC),output)
    }
    output
}

#----------------------------------------------------------------------
### AUC.survival.R ends here
