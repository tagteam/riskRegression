### riskQuantile.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  9 2016 (19:31) 
## Version: 
## last-updated: Jun  6 2016 (10:39) 
##           By: Thomas Alexander Gerds
##     Update #: 145
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
riskQuantile.binary <- function(DT,test,alpha,N,NT,NF,dolist,Q,...){
    reference=model=ReSpOnSe=risk=cause=X=NULL
    models <- unique(DT[,model])
    if (missing(Q)) Q <- c(0.05,0.25,0.5,0.75,0.95)
    else Q <- sort(Q)
    score.event <- DT[ReSpOnSe==1,data.table(t(quantile(risk,probs=Q))),by=list(model)]
    score.event[,cause:="event"]
    score.eventfree <- DT[ReSpOnSe==0,data.table(t(quantile(risk,probs=Q))),by=list(model)]
    score.eventfree[,cause:="event-free"]
    score.overall <- DT[,data.table(t(quantile(risk,probs=Q))),by=list(model)]
    score.overall[,cause:="overall"]
    score <- rbindlist(list(score.overall,score.event,score.eventfree))
    setcolorder(score,c("model","cause",names(score)[-c(1,length(names(score)))]))
    setnames(score,c("model","cause",gsub("^0","Q",as.character(Q))))
    if (length(dolist)>0){
        test <- data.table::rbindlist(lapply(dolist,function(g){
            ## FIXME: when dolist is 0:1 and models are 0:2 this does not work 
            DTdiff <- DT[model>g]
            ## from all models, substract risk of model g 
            DTdiff[,X:=DT[model==g,risk]-risk]
            N <- NROW(DTdiff)
            Xrange <- DTdiff[,range(X)]
            Xmed <- DTdiff[,median(X)]
            changedist.event <- DTdiff[ReSpOnSe==1,data.table(t(quantile(X,probs=Q))),by=list(model)]
            changedist.event[,cause:="event"]
            changedist.eventfree <- DTdiff[ReSpOnSe==0,data.table(t(quantile(X,probs=Q))),by=list(model)]
            changedist.eventfree[,cause:="event-free"]
            changedist.overall <- DTdiff[,data.table(t(quantile(X,probs=Q))),by=list(model)]
            changedist.overall[,cause:="overall"]
            changedist <- rbindlist(list(changedist.overall,changedist.event,changedist.eventfree))
            setcolorder(changedist,c("model","cause",names(changedist)[-c(1,length(names(changedist)))]))
            setnames(changedist,c("model","cause",gsub("^0","Q",as.character(Q))))
            changedist[,reference:=g]
            data.table::setcolorder(changedist,c("reference",colnames(changedist)[-length(colnames(changedist))]))
            changedist
        }))
    }else test <- NULL
    list(score=score,test=test)
}



riskQuantile.survival <- function(DT,test,alpha,N,NT,NF,dolist,Q,...){
    model=event=X=reference=status=times=cause=risk=Wt=WTi=NULL
    models <- unique(DT[,model])
    if (missing(Q)) Q <- c(0.05,0.25,0.5,0.75,0.95)
    else Q <- sort(Q)
    ##
    ## retrospectively (looking back from time t)
    ## we compute conditional on outcome the quantiles of predicted risks 
    ## and quantiles of changes of predicted risks
    ##
    ## Q=c(L=0.25,M=0.5,U=0.75)){
    #######
    ## Let X denote the difference between two 10-year predictions for the same subjects.
    ## This function estimates quantiles of the conditional distibutions of X given
    ## - event until time 10 years
    ## - event free until time 10 years
    ##
    ## For event the estimate is based on the following formula:
    ##
    ## P(X<=x|T<=t) =  P(X<=x, T<=t)    /   P(T<=t)
    ##              =: W(t,x)           /   F(t)
    ##
    ## 1. Direct estimate:
    ##
    ## We estimate F(t) with Kaplan-Meier using all data. For W
    ## we consider a plug-in estimate based on:
    ##
    ##  W(t,x) = \int F(t|X<=x) H(dx)
    ##
    ##  where we substitute the conditional Kaplan-Meier estimate for F(t|X<=x)
    ##  and the empirical distribution for H, ie, Hn= 1/n \sum_i 1{X_i==x}
    ##
    ## 2. IPCW estimate:
    ##
    ##  W_n(d,t) = (1/N) * \sum_i 1(X_i<d, T_i<t) * 1(Delta_i=1)/G(T_i-,X_i)
    ##
    ## where G is the reverse Kaplan-Meier evaluated at T_i- and possibly conditional on X_i.
    ##
    ## For 'event-free analyses' P(X<=x|T>t) is estimated by P(T>t|X<=x) P(X<=x)/P(T>t)
    #######
    surv <- DT[model==models[[1]],1/N*sum((time>times)/Wt)]
    cuminc <- 1-surv
    getQ.event <- function(Q,times,X,time,event,WTi,cuminc){
        Wx <- function(x,q,times,X,Xmed,time,event,WTi,cuminc){
            out <- 1/N*sum((X<=x & time<=times)/WTi)/cuminc-q
            if (is.na(out)) ifelse(x<Xmed,0,1) else out
        }
        Xrange <- range(X)
        Xmed <- median(X)
        qRisk <- sapply(Q,function(q){
            ## FIXME: not appropriate when X is discrete
            x <- try(u <- uniroot(Wx,interval=Xrange,q=q,times=times,X=X,Xmed=Xmed,time=time,event=event,WTi=WTi,cuminc=cuminc)$root,silent=TRUE)
            if ("try-error"%in%class(x)) as.numeric(NA) else u
        })
        qR <- data.table(t(qRisk))
        qR[,cause:="event"]
        qR
    }
    getQ.eventFree <- function(Q,times,X,Xmed,time,Wt,surv){
        Wx0 <- function(x,q,times,X,Xmed,time,Wt,surv){
            out <- 1/N*sum((X<=x & time>times)/Wt)/surv-q
            if (is.na(out)) ifelse(x<Xmed,0,1) else out
        }
        Xrange <- range(X)
        Xmed <- median(X)
        qRisk <- sapply(Q,function(q){
            ## FIXME: not appropriate when X is discrete
            x <- try(u <- uniroot(Wx0,interval=Xrange,q=q,times=times,X=X,Xmed=Xmed,time=time,Wt=Wt,surv=surv)$root,silent=TRUE)
            if ("try-error"%in%class(x)) as.numeric(NA) else u
        })
        qR <- data.table(t(qRisk))
        qR[,cause:="event-free"]
        qR
    }
    score.eventfree <- DT[,getQ.eventFree(Q=Q,times=times,X=risk,time=time,Wt=Wt,surv=surv),by=list(model,times)]
    score.event <- DT[,getQ.event(Q=Q,times=times,X=risk,time=time,event=event,WTi=WTi,cuminc=cuminc),by=list(model,times)]
    score.overall <- DT[,data.table(t(quantile(risk,probs=Q))),by=list(model,times)]
    score.overall[,cause:="overall"]
    colnames(score.overall) <- colnames(score.event)
    score <- rbindlist(list(score.overall,score.event,score.eventfree))
    setcolorder(score,c("model","times","cause",paste0("V",1:length(Q))))
    setnames(score,c("model","times","cause",gsub("^0","Q",as.character(Q))))
    if (length(dolist)>0){
        test <- data.table::rbindlist(lapply(dolist,function(g){
            DTdiff <- DT[model>g]
            ## from all models, substract risk from risk of model g 
            DTdiff[,X:=DT[model==g,risk]-risk]
            setorder(DTdiff, time,-status)
            N <- NROW(DTdiff)
            Xrange <- DTdiff[,range(X)]
            Xmed <- DTdiff[,median(X)]
            changedist.eventfree <- DTdiff[,getQ.eventFree(Q=Q,times=times,X=X,time=time,Wt=Wt,surv=surv),by=list(model,times)]
            changedist.event <- DTdiff[,getQ.event(Q=Q,times=times,X=X,time=time,event=event,WTi=WTi,cuminc=cuminc),by=list(model,times)]
            changedist.overall <- DTdiff[,data.table(t(quantile(X,probs=Q))),by=list(model,times)]
            changedist.overall[,cause:="overall"]
            colnames(changedist.overall) <- colnames(changedist.event)
            changedist <- rbindlist(list(changedist.overall,changedist.event,changedist.eventfree))
            setcolorder(changedist,c("model","times","cause",paste0("V",1:length(Q))))
            setnames(changedist,c("model","times","cause",gsub("^0","Q",as.character(Q))))
            changedist[,reference:=g]
            data.table::setcolorder(changedist,c("reference",colnames(changedist)[-length(colnames(changedist))]))
            changedist
        }))
    }else test <- NULL
    list(score=score,test=test)
}

riskQuantile.competing.risks <- function(DT,test,alpha,N,NT,NF,dolist,cause,Q,...){
    model=event=reference=risk=X=status=times=Wt=WTi=cause=NULL
    models <- unique(DT[,model])
    if (missing(Q)) Q <- c(0.05,0.25,0.5,0.75,0.95)
    else Q <- sort(Q)
    ##
    ## retrospectively (looking back from time t)
    ## we compute conditional on outcome the quantiles of predicted risks 
    ## and quantiles of changes of predicted risks
    ##
    ## if (missing(causes)) 
    causes <- DT[model==models[[1]] & status!=0,sort(unique(event))]
    ## (predrisk,
    ## formula,
    ## data,
    ## cause,
    ## times,
    ## ipcw=FALSE,
    ## Q=c(L=0.25,M=0.5,U=0.75)){
    #######
    ## Let X denote the difference between two 10-year predictions for the same subjects.
    ## This function estimates quantiles of the conditional distibutions of X given
    ## - event of cause 1 until time 10 years
    ## - event of cause 2 until time 10 years
    ## - event of cause 3 until time 10 years
    ## etc
    ## - event free until time 10 years
    ##
    ## For cause j the estimate is based on the following formula:
    ##
    ## P(X<=x|T<=t, cause=j) = P(X<=x, T<=t, cause=j) / P(T<=t, cause=j)
    ##                       =:      W(t,j,x)         /    F(t,j)
    ##
    ## 1. Direct estimate:
    ##
    ## We estimate F(t,j) with Aalen-Johansen using all data. For W
    ## we consider a plug-in estimate based on:
    ##
    ##  W(t,j,x) = \int F(t,j|X<=x) H(dx)
    ##
    ##  where we substitute the conditional Aalen-Johansen estimate for F(t,j|X<=x)
    ##  and the empirical distribution for H, ie, Hn= 1/n \sum_i 1{X_i==x}
    ##
    ##
    ## 2. IPCW estimate:
    ##
    ##  W_n(d,t,j) = (1/N) * \sum_i 1(X_i<d, cause_i=j, T_i<t) * 1(Delta_i=1)/G(T_i-,X_i)
    ##
    ## where G is the reverse Kaplan-Meier evaluated at T_i- and possibly conditional on X_i.
    ##
    ## For 'event-free analyses' P(X<=x|T>t) is estimated by P(T>t|X<=x) P(X<=x)/P(T>t)
    #######
    surv <- DT[model==models[[1]],1/N*sum((time>times)/Wt)]
    cuminc <- lapply(causes,function(cc)DT[model==models[[1]],1/N*sum((event==cc & time<=times)/WTi)])
    names(cuminc) <- causes
    getQ.causes <- function(Q,times,X,time,event,WTi,cuminc,causes){
        Wx <- function(x,q,times,X,Xmed,time,event,WTi,cuminc,cause){
            out <- 1/N*sum((X<=x & event==cause & time<=times)/WTi)/cuminc[[cause]]-q
            if (is.na(out)) ifelse(x<Xmed,0,1) else out
        }
        Xrange <- range(X)
        Xmed <- median(X)
        rbindlist(lapply(causes,function(cause){
            Q.cause <- sapply(Q,function(q){
                ## FIXME: not appropriate when X is discrete
                x <- try(u <- uniroot(Wx,interval=Xrange,q=q,times=times,X=X,Xmed=Xmed,time=time,event=event,WTi=WTi,cuminc=cuminc,cause=cause)$root,silent=TRUE)
                if ("try-error"%in%class(x)) as.numeric(NA) else u
            })
            qR <- data.table(t(Q.cause))
            qR[,cause:=cause]
            qR
            ## data.table(t(c(cause=cause,Q.cause)))
        }))
    }
    getQ.eventFree <- function(Q,times,X,Xmed,time,Wt,surv){
        Wx0 <- function(x,q,times,X,Xmed,time,Wt,surv){
            out <- 1/N*sum((X<=x & time>times)/Wt)/surv-q
            if (is.na(out)) ifelse(x<Xmed,0,1) else out
        }
        Xrange <- range(X)
        Xmed <- median(X)
        qRisk <- sapply(Q,function(q){
            ## FIXME: not appropriate when X is discrete
            x <- try(u <- uniroot(Wx0,interval=Xrange,q=q,times=times,X=X,Xmed=Xmed,time=time,Wt=Wt,surv=surv)$root,silent=TRUE)
            if ("try-error"%in%class(x)) as.numeric(NA) else u
        })
        qR <- data.table(t(qRisk))
        qR[,cause:="event-free"]
        qR
    }
    score.eventfree <- DT[,getQ.eventFree(Q=Q,times=times,X=risk,time=time,Wt=Wt,surv=surv),by=list(model,times)]
    score.causes <- DT[,getQ.causes(Q=Q,times=times,X=risk,time=time,event=event,WTi=WTi,cuminc=cuminc,causes=causes),by=list(model,times)]
    score.overall <- DT[,data.table(t(quantile(risk,probs=Q))),by=list(model,times)]
    score.overall[,cause:="overall"]
    colnames(score.overall) <- colnames(score.causes)
    score <- rbindlist(list(score.overall,score.causes,score.eventfree))
    setcolorder(score,c("model","times","cause",paste0("V",1:length(Q))))
    setnames(score,c("model","times","cause",gsub("^0","Q",as.character(Q))))
    ## browser(skipCalls=1)
    if (length(dolist)>0){
        test <- data.table::rbindlist(lapply(dolist,function(g){
            DTdiff <- DT[model>g]
            ## from all models, substract risk of model g 
            DTdiff[,X:=risk-DT[model==g,risk]]
            setorder(DTdiff, time,-status)
            N <- NROW(DTdiff)
            Xrange <- DTdiff[,range(X)]
            Xmed <- DTdiff[,median(X)]
            changedist.eventfree <- DTdiff[,getQ.eventFree(Q=Q,times=times,X=X,time=time,Wt=Wt,surv=surv),by=list(model,times)]
            changedist.causes <- DTdiff[,getQ.causes(Q=Q,times=times,X=X,time=time,event=event,WTi=WTi,cuminc=cuminc,causes=causes),by=list(model,times)]
            changedist.overall <- DTdiff[,data.table(t(quantile(X,probs=Q))),by=list(model,times)]
            changedist.overall[,cause:="overall"]
            colnames(changedist.overall) <- colnames(changedist.causes)
            changedist <- rbindlist(list(changedist.overall,changedist.causes,changedist.eventfree))
            setcolorder(changedist,c("model","times","cause",paste0("V",1:length(Q))))
            setnames(changedist,c("model","times","cause",gsub("^0","Q",as.character(Q))))
            changedist[,reference:=g]
            data.table::setcolorder(changedist,c("reference",colnames(changedist)[-length(colnames(changedist))]))
            changedist
        }))
    }else test <- NULL
    list(score=score,test=test)
}




#----------------------------------------------------------------------
### riskQuantile.R ends here
