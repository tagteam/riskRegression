### Brier.competing.risks.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:04)
## Version:
## Last-Updated: Jun  5 2024 (07:24) 
##           By: Thomas Alexander Gerds
##     Update #: 9
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

Brier.competing.risks <- function(DT,
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
                                  keep.residuals=FALSE,
                                  cause,
                                  states,
                                  IC.data,
                                  ...){
    IC0=nth.times=riskRegression_ID=riskRegression_time=times=riskRegression_event=Brier=raw.Residuals=risk=residuals=WTi=Wt=riskRegression_status=setorder=model=IF.Brier=data.table=sd=lower=qnorm=se=upper=NULL
    ## compute 0/1 outcome:
    thecause <- match(cause,states,nomatch=0)
    if (length(thecause)==0) stop("Cannot identify cause of interest")
    DT[riskRegression_time<=times & riskRegression_status==1 & riskRegression_event==thecause,residuals:=(1-risk)^2/WTi]
    DT[riskRegression_time<=times & riskRegression_status==1 & riskRegression_event!=thecause,residuals:=(risk)^2/WTi]
    DT[riskRegression_time<=times & riskRegression_status==0,residuals:=0]
    DT[riskRegression_time>times,residuals:=(risk)^2/Wt]
    ## deal with censored observations before 
    DT[riskRegression_time<=times & riskRegression_status==0,residuals:=0]
    if (se.fit[[1]]==1L){
        ## data.table::setorder(DT,model,times,riskRegression_time,-riskRegression_status)
        data.table::setorder(DT,model,times,riskRegression_ID)
        DT[,nth.times:=as.numeric(factor(times))]
        DT[,IC0:=residuals-mean(residuals),by=list(model,times)]
        DT[,IF.Brier:=getInfluenceCurve.Brier(t=times[1],
                                              time=riskRegression_time,
                                              IC0,
                                              residuals=residuals,
                                              IC.G=MC,
                                              cens.model=cens.model,
                                              conservative = conservative,
                                              nth.times=nth.times[1],
                                              event = riskRegression_status*riskRegression_event),by=list(model,times)]
        score <- DT[,data.table(Brier=sum(residuals)/N,
                                se=sd(IF.Brier)/sqrt(N)),by=list(model,times)]
        if (se.fit==TRUE){
            score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
            score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        }
    }else{
        ## no se.fit
        score <- DT[,data.table(Brier=sum(residuals)/N),by=list(model,times)]
    }
    data.table::setkey(score,model,times)
    if (length(dolist)>0){
        data.table::setkey(DT,model,times)
        ## merge with Brier score
        DT <- DT[score]
        data.table::setkey(score,model,times)
        if (se.fit[[1]]==TRUE){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,IF=IF.Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  se.fit=se.fit),by=list(times)]
        }else{
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  se.fit=FALSE),by=list(times)]
        }
        setnames(contrasts.Brier,"delta","delta.Brier")
        output <- list(score=score,contrasts=contrasts.Brier)
    } else{
        output <- list(score=score)
    }
    if (keep.residuals) {
        output <- c(output,list(residuals=DT[,c("riskRegression_ID","riskRegression_time","riskRegression_status","model","times","risk","residuals"),with=FALSE]))
    }
    if (keep.vcov[[1]] && se.fit[[1]]==TRUE){
        output <- c(output,list(vcov=getVcov(DT,"IF.Brier",times=TRUE)))
    }
    if (keep.iid[[1]] && se.fit[[1]] == TRUE) {
        output <- c(output,
                    list(iid.decomp = DT[,data.table::data.table(riskRegression_ID,model,cause,times,IF.Brier)]))
    }
    output
}

#----------------------------------------------------------------------
### Brier.competing.risks.R ends here
