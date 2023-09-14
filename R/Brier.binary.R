### Brier.binary.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:03) 
## Version: 
## Last-Updated: Jun 30 2023 (10:59) 
##           By: Thomas Alexander Gerds
##     Update #: 4
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
                         multi.split.test,
                         alpha,
                         N,
                         NT,
                         NF,
                         dolist,
                         keep.residuals=FALSE,
                         ...){
    residuals=Brier=risk=model=ReSpOnSe=lower=upper=se=ID=IF.Brier=NULL
    DT[,residuals:=(ReSpOnSe-risk)^2]
    if (se.fit==1L){
        data.table::setkey(DT,model,ID)
        score <- DT[,data.table::data.table(Brier=sum(residuals)/N,se=sd(residuals)/sqrt(N)),by=list(model)]
        score[,lower:=pmax(0,Brier-qnorm(1-alpha/2)*se)]
        score[,upper:=pmin(1,Brier + qnorm(1-alpha/2)*se)]
        ## Influence function IPA
        ## Brier.null <- score[model==0][["Brier"]]
        ## Brier.model <- score[model!=0][["Brier"]]
        ## IC.ipa <- -(1/Brier.null)* DT[model!=0][["residuals"]] + DT[model==0][["residuals"]]*Brier.model/(Brier.null)^2
        if (keep.vcov){
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
        data.table::setkey(score,model)
        DT <- DT[score]
        if (se.fit[[1]]==TRUE || multi.split.test[[1]]==TRUE){
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,
                                                             IF=residuals,
                                                             model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=multi.split.test,
                                                  se.fit=se.fit)]
        }else{
            contrasts.Brier <- DT[,getComparisons(data.table(x=Brier,model=model),
                                                  NF=NF,
                                                  N=N,
                                                  alpha=alpha,
                                                  dolist=dolist,
                                                  multi.split.test=FALSE,
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
    if (keep.iid[1] && se.fit[1] == TRUE) {
        output <- c(output,
                    list(iid.decomp = DT[,data.table::data.table(ID,model,IF.Brier)]))
    }
    if (keep.residuals[[1]] == TRUE) {
        output <- c(output,list(residuals=DT[,data.table::data.table(ID,ReSpOnSe,model,risk,residuals)]))
    }
    output
}

#----------------------------------------------------------------------
### Brier.binary.R ends here
