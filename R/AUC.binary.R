### AUC.binary.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds and Johan Sebastian Ohlendorff
## Created: Jan 11 2022 (17:04) 
## Version: 
## Last-Updated: May 15 2025 (08:49) 
##           By: Thomas Alexander Gerds
##     Update #: 82
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# breaks are used to approximate the ROC curve and AUC
# cutpoints are used to estimate Sens, Spec, PPV, NPV
AUC.binary <- function(DT,
                       breaks=NULL,
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
                       ROC,
                       cutpoints,
                       ...){
    prob_predicted_positive = prob_predicted_negative = PPV=NPV=count_predicted_positive=count_predicted_negative=model=risk=riskRegression_event=FPR=TPR=riskRegression_ID=NULL
    auRoc.numeric <- function(X,D,breaks,ROC){
        ## need to reverse when high X is concordant with {response=1}
        if (is.null(breaks)) breaks <- rev(sort(unique(X))) 
        TPR <- c(prodlim::sindex(jump.times=X[D==1],
                                 eval.times=breaks,
                                 comp="greater",
                                 strict=FALSE)/sum(D==1))
        FPR <- c(prodlim::sindex(jump.times=X[D==0],
                                 eval.times=breaks,
                                 comp="greater",
                                 strict=FALSE)/sum(D==0))
        if (ROC) {
            data.table(risk=breaks,TPR,FPR)
        } else {
            0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
        }
    }
    auRoc.factor <- function(X,D,ROC){
        TPR <- (sum(D==1)-table(X[D==1]))/sum(D==1)
        FPR <- table(X[D==0])/sum(D==0)
        if (ROC==TRUE)
            data.table(cbind(risk=c(sort(unique(X))),TPR,FPR))
        else
            0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
    }
    aucDT <- DT[model>0]
    dolist <- dolist[sapply(dolist,function(do){match("0",do,nomatch=0L)})==0]
    data.table::setkey(aucDT,model,riskRegression_ID)
    if (is.factor(DT[["risk"]])){
        score <- aucDT[,auRoc.factor(risk,
                                     riskRegression_event,
                                     ROC=ROC),by=list(model)]
    } else{
        score <- aucDT[,auRoc.numeric(X = risk,
                                      D = riskRegression_event,
                                      breaks=breaks,
                                      ROC=ROC),by=list(model)]
    }
    if (ROC==FALSE){
        setnames(score,"V1","AUC")
        output <- list(score=score)
    } else {
        AUC <- score[,list(AUC=0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))),by=list(model)]
        ROC <- score
        output <- list(score=AUC,ROC=ROC)
    }
    if(!is.null(cutpoints)){
        cutpoint_results <- aucDT[,{
            # predicted positive: Marker >= cutoff
            # predicted negative: Marker < cutoff
            # if X = c(7,10,77) and breaks = 7 then 1 is predicted
            TPR <- c(prodlim::sindex(jump.times=risk[riskRegression_event==1],eval.times=cutpoints,comp="greater",strict=FALSE)/sum(riskRegression_event==1))
            FPR <- c(prodlim::sindex(jump.times=risk[riskRegression_event==0],eval.times=cutpoints,comp="greater",strict=FALSE)/sum(riskRegression_event==0))
            count_predicted_positive <- prodlim::sindex(jump.times=risk,eval.times=cutpoints,comp="greater",strict=FALSE)
            count_predicted_negative <- prodlim::sindex(jump.times=risk,eval.times=cutpoints,comp="smaller",strict=TRUE)
            PPV <- prodlim::sindex(jump.times=risk[riskRegression_event==1],eval.times=cutpoints,comp="greater",strict=FALSE)/count_predicted_positive
            # if no one is predicted positive then the positive predictive rate is 1
            PPV[count_predicted_positive == 0] <- 1
            NPV <- prodlim::sindex(jump.times=risk[riskRegression_event==0],eval.times=cutpoints,comp="smaller",strict=TRUE)/count_predicted_negative
            # if no one is predicted negative then the negative predictive rate is 0
            NPV[count_predicted_negative == 0] <- 1
            prob_predicted_positive <- count_predicted_positive/.N
            prob_predicted_negative <- count_predicted_negative/.N
            # standard errors
            meanY <- mean(riskRegression_event)
            se.TPR = sapply(1:length(cutpoints),function(i){sd(riskRegression_event/meanY * ((risk > cutpoints[[i]])-TPR[[i]]))/sqrt(.N)})
            se.FPR = sapply(1:length(cutpoints),function(i){sd((1-riskRegression_event)/(1-meanY) * ((risk > cutpoints[[i]])-FPR[[i]]))/sqrt(.N)})
            se.PPV = sapply(1:length(cutpoints),function(i){sd((risk > cutpoints[[i]])/prob_predicted_positive[[i]] * (riskRegression_event - PPV[[i]]))/sqrt(.N)})
            se.NPV = sapply(1:length(cutpoints),function(i){sd((risk <= cutpoints[[i]])/prob_predicted_negative[[i]] * ((1-riskRegression_event) - NPV[[i]]))/sqrt(.N)})
            data.table(cutpoint=cutpoints,TPR,se.TPR,FPR,se.FPR,PPV,se.PPV,NPV,se.NPV)
        },by=list(model)]
        output <- c(output, list(cutpoints=cutpoint_results[]))
    }
    if (length(dolist)>0 || (se.fit[[1]]==1L)){
        ## do not want to depend on Daim as they turn marker to ensure auc > 0.5
        delongtest <-  function(risk,score,dolist,response,cause,alpha,se.fit,keep.vcov) {
            cov=lower=upper=p=AUC=se=lower=upper=NULL
            if (keep.iid == TRUE){warning("Argument 'keep.iid' is ignored. This function does not explicitely calculate the estimated influence function for AUC with binary outcome.")}
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
            if (length(dolist)>0){
                ncomp <- length(dolist)
                delta.AUC <- numeric(ncomp)
                se <- numeric(ncomp)
                model <- numeric(ncomp)
                reference <- numeric(ncomp)
                ctr <- 1
                Qnorm <- qnorm(1 - alpha/2)
                for (d in dolist){
                    i <- d[1]
                    for (j in d[-1]) {
                        delta.AUC[ctr] <- auc[j]-auc[i]
                        if (se.fit[[1]]){
                            LSL <- t(c(1, -1)) %*% S[c(j, i), c(j, i)] %*% c(1, -1)
                            se[ctr] <- sqrt(LSL)
                        }
                        model[ctr] <- modelnames[j]
                        reference[ctr] <- modelnames[i]
                        ctr <- ctr + 1
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
            }else{
                out <- list(score = score, contrasts = NULL)
            }
            ## should only be kept if se.fit is true
            if (keep.vcov && se.fit[[1]]==TRUE) {
                out <- c(out,list(vcov=S))
            }
            out
        }
        xRisk <- data.table::dcast(aucDT[],riskRegression_ID~model,value.var="risk")[,-1,with=FALSE]
        delong.res <- delongtest(risk=xRisk,
                                 score=output$score,
                                 dolist=dolist,
                                 response=aucDT[model==model[1],riskRegression_event],
                                 cause="1",
                                 alpha=alpha,
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
