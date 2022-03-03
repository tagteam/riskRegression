### AUC.survival.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:06) 
## Version: 
## Last-Updated: Jan 11 2022 (17:23) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

AUC.survival <- function(DT,MC,se.fit,conservative,cens.model,keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,ROC,...){
    ID=model=times=risk=Cases=time=status=Controls=TPR=FPR=WTi=Wt=ipcwControls=ipcwCases=IF.AUC=lower=se=upper=AUC=NULL
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
        ## compute influence function
        ## data.table::setorder(aucDT,model,times,time,-status)
        data.table::setorder(aucDT,model,times,ID)
        if (conservative==TRUE && (cens.model == "KaplanMeier" || cens.model == "none")){
            aucDT[,IF.AUC:=getInfluenceFunctionAUCSurvival(time,status,times[1],risk,WTi,Wt[1],AUC[1]), by=list(model,times)]
        }
        else {
            aucDT[,IF.AUC:=getInfluenceCurve.AUC.survival(t=times[1],
                                                           n=N,
                                                           time=time,
                                                           risk=risk,
                                                           Cases=Cases,
                                                           Controls=Controls,
                                                           ipcwControls=ipcwControls,
                                                           ipcwCases=ipcwCases,
                                                           MC=MC), by=list(model,times)]
        }
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
