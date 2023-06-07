### AUC.survival.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:06)
## Version:
## Last-Updated: Nov 21 2022 (13:55) 
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

AUC.survival <- function(DT,breaks=NULL,MC,se.fit,conservative,cens.model,keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,ROC,IC.data,...){
    ID=model=times=risk=Cases=time=status=Controls=TPR=FPR=WTi=Wt=ipcwControls=ipcwCases=IF.AUC=lower=se=upper=AUC=nth.times=NULL
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
        ## data.table::setorder(aucDT,model,times,time,-status)
        data.table::setorder(aucDT,model,times,ID)
        aucDT[,IF.AUC:=getInfluenceCurve.AUC(times[1],
                                             time,
                                             status,
                                             WTi,
                                             Wt,
                                             risk,
                                             MC,
                                             AUC[1],
                                             nth.times[1],
                                             conservative[[1]],
                                             cens.model), by=list(model,times)]
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
        if (se.fit[[1]]==TRUE || multi.split.test[[1]]==TRUE){
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

#----------------------------------------------------------------------
### AUC.survival.R ends here
