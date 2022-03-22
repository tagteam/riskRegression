### AUC.binary.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:04) 
## Version: 
## Last-Updated: Mar 22 2022 (13:27) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

AUC.binary <- function(DT,breaks=NULL,se.fit,conservative=FALSE,cens.model="none",keep.vcov=FALSE,multi.split.test,alpha,N,NT,NF,dolist,ROC,...){
    model=risk=ReSpOnSe=FPR=TPR=ID=NULL

    # {{{  helper functions
    AireTrap <- function(FP,TP,N){
        N <- length(FP)
        sum((FP-c(0,FP[-N]))*((c(0,TP[-N])+TP)/2))
    }

    ## do not want to depend on Daim as they turn marker to ensure auc > 0.5
    delongtest <-  function(risk,
                            score,
                            dolist,
                            response,
                            cause,
                            alpha,
                            multi.split.test,
                            se.fit,
                            keep.vcov) {
        cov=lower=upper=p=AUC=se=lower=upper=NULL
        auc <- score[["AUC"]]
        nauc <- ncol(risk)
        modelnames <- score[["model"]]
        score <- data.table(model=colnames(risk),AUC=auc)
        if (se.fit==1L){
            Cases <- response == cause
            #Controls <- response != cause
            #riskcontrols <- as.matrix(risk[Controls,])
            riskcontrols <- as.matrix(risk[!Cases,])
            riskcases <- as.matrix(risk[Cases,])

        # new method, uses a fast implementation of delongs covariance matrix
        # Fast Implementation of DeLong’s Algorithm for Comparing the Areas Under Correlated Receiver Operating Characteristic Curves
        # article can be found here:
        # https://ieeexplore.ieee.org/document/6851192
        S <- calculateDelongCovarianceFast(riskcases,riskcontrols)
        se.auc <- sqrt(diag(S))
        score[,se:=se.auc]
        score[,lower:=pmax(0,AUC-qnorm(1-alpha/2)*se)]
        score[,upper:=pmin(1,AUC+qnorm(1-alpha/2)*se)]
        setcolorder(score,c("model","AUC","se","lower","upper"))
    }else{
        setcolorder(score,c("model","AUC"))
    }
    names(auc) <- 1:nauc
    ## q1 <- auc/(2 - auc)
    ## q2 <- 2 * auc^2/(1 + auc)
    ## aucvar <- (auc * (1 - auc) + (nCases - 1) * (q1 - auc^2) + (nControls - 1) * (q2 - auc^2))/(nCases * nControls)
    if (length(dolist)>0){
        ## ncomp <- nauc * (nauc - 1)/2
        ncomp <- length(dolist)
        delta.AUC <- numeric(ncomp)
        se <- numeric(ncomp)
        model <- numeric(ncomp)
        reference <- numeric(ncomp)
        ctr <- 1

        Qnorm <- qnorm(1 - alpha/2)
        for (d in dolist){
            i <- d[1]
            ## for (i in 1:(nauc - 1)) {
            ## for (j in (i + 1):nauc) {
            for (j in d[-1]) {
              delta.AUC[ctr] <- auc[j]-auc[i]
              if (se.fit[[1]]){
                ## cor.auc[ctr] <- S[i, j]/sqrt(S[i, i] * S[j, j])
                LSL <- t(c(1, -1)) %*% S[c(j, i), c(j, i)] %*% c(1, -1)
                ## print(c(1/LSL,rms::matinv(LSL)))
                se[ctr] <- sqrt(LSL)
              }
              ## tmpz <- (delta.AUC[ctr]) %*% rms::matinv(LSL) %*% delta.AUC[ctr]
              ## tmpz <- (delta.AUC[ctr]) %*% (1/LSL) %*% delta.AUC[ctr]
              model[ctr] <- modelnames[j]
              reference[ctr] <- modelnames[i]
              ctr <- ctr + 1
              ## }
            }
        }
        deltaAUC <- data.table(model,reference,delta.AUC=as.vector(delta.AUC))
        if (se.fit[[1]]){
          deltaAUC[,se:=se]
          deltaAUC[,lower:=delta.AUC-Qnorm*se]
          deltaAUC[,upper:=delta.AUC+Qnorm*se]
          deltaAUC[,p:=2*pnorm(abs(delta.AUC)/se,lower.tail=FALSE)]
        }
        out <- list(score = score, contrasts = deltaAUC)
        ##if (se.fit[[1]]==TRUE||multi.split.test[[1]]==TRUE){
        ##    deltaAUC <- data.table(model,reference,delta.AUC=as.vector(delta.AUC),se)
        ##}else{
        ##    deltaAUC <- data.table(model,reference,delta.AUC=as.vector(delta.AUC))
        ##}
    }else{
        out <- list(score = score, contrasts = NULL)
    }
    ## should only be kept if se.fit is true
    if (keep.vcov && se.fit[[1]]==TRUE) {
        out <- c(out,list(vcov=S))
    }
    out
}

    auRoc.numeric <- function(X,D,breaks,ROC){
        if (is.null(breaks)) breaks <- rev(sort(unique(X))) ## need to reverse when high X is concordant with {response=1}
        TPR <- c(prodlim::sindex(jump.times=X[D==1],eval.times=breaks,comp="greater",strict=FALSE)/sum(D==1))
        FPR <- c(prodlim::sindex(jump.times=X[D==0],eval.times=breaks,comp="greater",strict=FALSE)/sum(D==0))
        if (ROC==TRUE)
            data.table(risk=breaks,TPR,FPR)
        else
            0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
    }
    auRoc.factor <- function(X,D,ROC){
        TPR <- (sum(D==1)-table(X[D==1]))/sum(D==1)
        FPR <- table(X[D==0])/sum(D==0)
        if (ROC==TRUE)
            data.table(cbind(risk=c(sort(unique(X))),TPR,FPR))
        else
            0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
    }
    # }}}

    aucDT <- DT[model>0]
    dolist <- dolist[sapply(dolist,function(do){match("0",do,nomatch=0L)})==0]
    data.table::setkey(aucDT,model,ID)
    if (is.factor(DT[["risk"]])){
        score <- aucDT[,auRoc.factor(risk,ReSpOnSe,ROC=ROC),by=list(model)]
    }
    else{
        score <- aucDT[,auRoc.numeric(risk,ReSpOnSe,breaks=NULL,ROC=ROC),by=list(model)]
    }
    if (ROC==FALSE){
        setnames(score,"V1","AUC")
        output <- list(score=score)
    } else{
        AUC <- score[,list(AUC=0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))),by=list(model)]
        ROC <- score
        output <- list(score=AUC,ROC=ROC)
    }
    if (length(dolist)>0 || (se.fit[[1]]==1L)){
        xRisk <- data.table::dcast(aucDT[],ID~model,value.var="risk")[,-1,with=FALSE]
        delong.res <- delongtest(risk=xRisk,
                                 score=output$score,
                                 dolist=dolist,
                                 response=aucDT[model==model[1],ReSpOnSe],
                                 cause="1",
                                 alpha=alpha,
                                 multi.split.test=multi.split.test,
                                 se.fit=se.fit,
                                 keep.vcov=keep.vcov)
        output$score <- delong.res$score
        output$contrasts <- delong.res$contrasts
        if (keep.vcov){
            output$vcov <- delong.res$vcov
        }
        output
    }else{
        output
    }
}



#----------------------------------------------------------------------
### AUC.binary.R ends here
