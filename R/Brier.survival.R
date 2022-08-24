### Brier.survival.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:04)
## Version:
## Last-Updated: Jan 11 2022 (17:04)
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

Brier.survival <- function(DT,MC,se.fit,conservative,cens.model,keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,keep.residuals=FALSE,old.ic.method,IC.data,...){
    IC0=IPCW=nth.times=ID=time=times=raw.Residuals=risk=Brier=residuals=WTi=Wt=status=setorder=model=IF.Brier=data.table=sd=lower=qnorm=se=upper=NULL
    ## compute 0/1 outcome:
    DT[time<=times & status==1,residuals:=(1-risk)^2/WTi]
    DT[time<=times & status==0,residuals:=0]
    DT[time>times,residuals:=(risk)^2/Wt]
    
    if (se.fit[[1]]==1L || multi.split.test[[1]]==TRUE){
        ## data.table::setorder(DT,model,times,time,-status)
        data.table::setorder(DT,model,times,ID)
        DT[,nth.times:=as.numeric(factor(times))]
        DT[,IC0:=residuals-mean(residuals),by=list(model,times)]
        if (conservative==TRUE){
            score <- DT[,data.table(Brier=sum(residuals)/N,
                                    se=sd(IC0)/sqrt(N)),by=list(model,times)]
            setnames(DT,"IC0","IF.Brier")
        }else{
            if (cens.model=="none"){
                DT[,IF.Brier:=residuals]
                score <- DT[,data.table(Brier=sum(residuals)/N,
                                        se=sd(residuals)/sqrt(N),
                                        se.conservative=sd(residuals)),by=list(model,times)]
            }else if (cens.model=="KaplanMeier"){
                DT[,Brier := sum(residuals)/N,by=list(model,times)]
                if (old.ic.method){
                    DT[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                                          time=time,
                                                          IC0,
                                                          residuals=residuals,
                                                          WTi=WTi,
                                                          Wt=Wt,
                                                          IC.G=MC,
                                                          cens.model=cens.model,
                                                          nth.times=nth.times[1]),by=list(model,times)]
                  DT[,IF.Brier2 := 0]
                }
                else {
                    browser()
                    DT[,IF.Brier := getInfluenceFunctionBrierKMCensoringUseSquared(times[1],time,residuals,status),by=list(model,times)]
                    # DT[,IF.Brier := getInfluenceFunctionBrierKMCensoring(times[1],time,risk,status,WTi,Brier[1]),by=list(model,times)]
                }
                score <- DT[,data.table(Brier=sum(residuals)/N,
                                        #se2  = sd(IF.Brier2)/sqrt(N),
                                        se=sd(IF.Brier)/sqrt(N),
                                        se.conservative=sd(IC0)/sqrt(N)),by=list(model,times)]
            }
            else {
                if (!old.ic.method){
                  browser()
                  DT[,IF.Brier:=getInfluenceCurve.Brier.covariates(times[1],time,risk,status,WTi,sum(residuals)/N,IC.data), by=list(model,times)]
                  DT[,IF.Brier2:=getInfluenceCurve.Brier.covariates.use.squared(times[1],time,residuals,risk,status,WTi,Wt,IC.data), by=list(model,times)]
                }
                else {
                  DT[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                                        time=time,
                                                        IC0,
                                                        residuals=residuals,
                                                        WTi=WTi,
                                                        Wt=Wt,
                                                        IC.G=MC,
                                                        cens.model=cens.model,
                                                        nth.times=nth.times[1]),by=list(model,times)]
                }
                score <- DT[,data.table(Brier=sum(residuals)/N,
                                        se=sd(IF.Brier)/sqrt(N)),by=list(model,times)]
            }
        }
        if (se.fit==TRUE){
            score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        }
    }else{
        ## no se.fit
        score <- DT[,data.table(Brier=sum(residuals)/N),by=list(model,times)]
    }
    data.table::setkey(score,model,times)
    if (length(dolist)>0L){
        ## merge with Brier score
        data.table::setkey(DT,model,times)
        ## data.table::setkey(score,model,times)
        DT <- DT[score]
        if (se.fit[[1]]==TRUE || multi.split.test[[1]]==TRUE){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=multi.split.test,
                                                  se.fit=se.fit),by=list(times)]
        }else{
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=FALSE,
                                                  se.fit=FALSE),by=list(times)]
        }
        setnames(contrasts.Brier,"delta","delta.Brier")
        output <- list(score=score,contrasts=contrasts.Brier)
    } else{
        output <- list(score=score)
    }
    if (keep.vcov[1] && se.fit[1]==TRUE){
        output <- c(output,list(vcov=getVcov(DT,"IF.Brier",times=TRUE)))
    }
    if (keep.residuals) {
        if (all(c("Wt","WTi")%in%names(DT))){
            DT[,IPCW:=1/WTi]
            DT[time>=times,IPCW:=1/Wt]
            DT[time<times & status==0,IPCW:=0]
            output <- c(output,list(residuals=DT[,c("ID","time","status","model","times","risk","residuals","IPCW"),with=FALSE]))
        }else{
            output <- c(output,list(residuals=DT[,c("ID","time","status","model","times","risk","residuals"),with=FALSE]))
        }
    }
    output
}

#----------------------------------------------------------------------
### Brier.survival.R ends here
