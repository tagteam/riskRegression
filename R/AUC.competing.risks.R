### AUC.competing.risks.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:06)
## Version:
## Last-Updated: May 31 2022 (11:37) 
##           By: Thomas Alexander Gerds
##     Update #: 19
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

AUC.competing.risks <- function(DT,MC,se.fit,conservative,cens.model,keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,cause,states,ROC,old.ic.method,IC.data,...){
    ID=model=times=risk=Cases=time=status=event=Controls1=Controls2=TPR=FPR=WTi=Wt=ipcwControls1=ipcwControls2=ipcwCases=IF.AUC=lower=se=upper=AUC=NULL
    aucDT <- DT[model>0]
    dolist <- dolist[sapply(dolist,function(do){match("0",do,nomatch=0L)})==0]
    ## assign Weights before ordering
    aucDT[,ipcwControls1:=1/(Wt*N)]
    aucDT[,ipcwControls2:=1/(WTi*N)]
    aucDT[,ipcwCases:=1/(WTi*N)]
    aucDT[,ipcwControls2:=1/(WTi*N)]
    ## order data
    data.table::setorder(aucDT,model,times,-risk)
    ## identify cases and controls
    thecause <- match(cause,states,nomatch=0)
    if (length(thecause)==0) stop("Cannot identify cause of interest")
    aucDT[,Cases:=(time <=times &  event==thecause)]
    aucDT[,Controls1:=(time > times)]
    aucDT[,Controls2:=(time <=times &  event!=thecause & status !=0)]
    ## prepare Weights
    aucDT[Cases==0,ipcwCases:=0]
    aucDT[Controls1==0,ipcwControls1:=0]
    aucDT[Controls2==0,ipcwControls2:=0]
    ## compute denominator
    ## ROC <- aucDT[,list(TPR=c(0,cumsum(ipcwCases)),FPR=c(0,cumsum(ipcwControls1)+cumsum(ipcwControls2))),by=list(model,times)]
    aucDT[,TPR:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    aucDT[,FPR:=(cumsum(ipcwControls1)+cumsum(ipcwControls2))/(sum(ipcwControls2)+sum(ipcwControls1)),by=list(model,times)]
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
        if (cens.model == "KaplanMeier" || cens.model == "none"){
            if (old.ic.method){
                aucDT[,IF.AUC:={
                    if (sum(Controls2)==0){
                        getInfluenceCurve.AUC.survival(t=times[1],
                                                       n=N,
                                                       time=time,
                                                       risk=risk,
                                                       Cases=Cases,
                                                       Controls=Controls1,
                                                       ipcwControls=ipcwControls1,
                                                       ipcwCases=ipcwCases,
                                                       MC=MC)
                    }else{
                        getInfluenceCurve.AUC.competing.risks(t=times[1],
                                                              n=N,
                                                              time=time,
                                                              risk=risk,
                                                              ipcwControls1=ipcwControls1,
                                                              ipcwControls2=ipcwControls2,
                                                              ipcwCases=ipcwCases,
                                                              Cases=Cases,
                                                              Controls1=Controls1,
                                                              Controls2=Controls2,
                                                              MC=MC)
                    }
                }, by=list(model,times)]
                # aucDT[,IF.AUC:=getInfluenceCurve.AUC.slow(times[1],N,time,status*event,risk,WTi,Wt[1],MC,AUC[1],FALSE)+0.5*getInfluenceCurve.AUC.slow(times[1],N,time,status*event,risk,WTi,Wt[1],MC,AUC[1],TRUE), by=list(model,times)]
                # browser()
            }
            else {
                aucDT[,IF.AUC:={
                    if (sum(Controls2)==0){
                        getInfluenceCurveHelper(time,status*event,times[1],risk,WTi,Wt[1],AUC[1])
                    }
                    else {
                        # getInfluenceFunctionAUCConservative(time,status*event,times[1],risk,WTi,Wt[1],score$AUC,FALSE)
                        getInfluenceCurveHelper(time,status*event,times[1],risk,WTi,Wt[1],AUC[1])#+0.5*getInfluenceFunctionAUC(time,status*event,times[1],risk,WTi,Wt[1],AUC[1],FALSE,TRUE,FALSE)
                    }
                }, by=list(model,times)]
                # shoud be zero with cont. covariates, the ties part that is
                # browser()
            }
        }
        else {
            # for now does not support ties
            if (conservative){
                aucDT[,IF.AUC:=getInfluenceCurve.AUC.covariates.conservative(times[1],N,time,status*event,risk,WTi,Wt,AUC[1]), by=list(model,times)]
            }
            else {
                aucDT[,IF.AUC:=getInfluenceCurve.AUC.covariates(times[1],N,time,status*event,risk,WTi,Wt,AUC[1],IC.data), by=list(model,times)]
            }
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
