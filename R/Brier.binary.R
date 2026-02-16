### Brier.binary.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:03) 
## Version: 
## Last-Updated: feb  6 2026 (12:08) 
##           By: Thomas Alexander Gerds
##     Update #: 14
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

Brier.binary <- function(DT,
                         se.fit,
                         conservative=FALSE,
                         cens.model="none",
                         keep.vcov=FALSE,
                         keep.iid=FALSE,
                         alpha,
                         N,
                         NT,
                         NF,
                         dolist,
                         keep.residuals=FALSE,
                         ...){
    residuals=Brier=risk=model=riskRegression_event=lower=upper=se=riskRegression_ID=IF.Brier=NULL
    DT[,residuals:=(riskRegression_event-risk)^2]
    if (se.fit==1L){
        data.table::setkey(DT,model,riskRegression_ID)
        score <- DT[,data.table::data.table(Brier=sum(residuals)/N,se=sd(residuals)/sqrt(N)),by=list(model)]
        score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
        score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        if (keep.iid | keep.vcov){
            DT[,Brier:=sum(residuals)/N,by=list(model)]
            DT[,IF.Brier:=residuals-Brier]
        }
    }else{
        ## no se.fit
        score <- DT[,data.table(Brier=sum(residuals)/N),by=list(model)]
    }
    data.table::setkey(score,model)
    if (length(dolist)>0){
        ## merge with Brier score
        data.table::setkey(DT,model)
        DT <- DT[score,,on = c("model")]
        data.table::setkey(DT,model)
        if (se.fit[[1]]==TRUE){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,
                                                             IF=residuals,
                                                             model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  se.fit=se.fit)]
        }else{
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  se.fit=FALSE)]
        }
        setnames(contrasts.Brier,"delta","delta.Brier")
        output <- list(score=score,contrasts=contrasts.Brier)
    }else{
        output <- list(score=score)
    }
    if (keep.vcov[[1]] == TRUE && se.fit[[1]]==TRUE){
        output <- c(output,list(vcov=getVcov(DT,"IF.Brier")))
    }
    if (keep.iid[[1]] && se.fit[[1]] == TRUE) {
        output <- c(output,
                    list(iid.decomp = DT[,data.table::data.table(riskRegression_ID,model,IF.Brier)]))
    }
    if (keep.residuals[[1]] == TRUE) {
        output <- c(output,list(residuals=DT[,data.table::data.table(riskRegression_ID,riskRegression_event,model,risk,residuals)]))
    }
    output
}

#----------------------------------------------------------------------
### Brier.binary.R ends here
