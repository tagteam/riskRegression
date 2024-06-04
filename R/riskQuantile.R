### riskQuantile.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  9 2016 (19:31) 
## Version: 
## last-updated: Jun  4 2024 (07:21) 
##           By: Thomas Alexander Gerds
##     Update #: 307
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getQuantile <- function(x,Fx,Q){
    ## Note: x has to be sorted in ascending order
    ## FIXME: not appropriate when X is discrete?, see help(quantile)
    q.index <- pmin(length(x),1+prodlim::sindex(jump.times=Fx,eval.times=Q,strict=TRUE))
    ## quantile type = 1 
    quant <- x[q.index]
    quant
}

riskQuantile.binary <- function(DT,N,NT,NF,dolist,Q,...){
    reference=model=ReSpOnSe=risk=cause=X=riskRegression_ID=NULL
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
    qnames <- paste0("Q",".",as.character(round(100*Q)))
    setnames(score,c("model","cause",qnames))
    if (length(dolist)>0){
        contrasts <- data.table::rbindlist(lapply(dolist,function(g){
            ## from all models g[-1], substract risk of model g[1] 
            setorder(DT,model,riskRegression_ID)
            DTdiff <- DT[model%in%g[-1]]
            DTref <- rep(DT[model==g[1],risk],length(unique(DTdiff$model)))
            DTdiff[,X:=risk-DTref]
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
            setnames(changedist,c("model","cause",qnames))
            changedist[,reference:=g[1]]
            data.table::setcolorder(changedist,c("reference",colnames(changedist)[-length(colnames(changedist))]))
            changedist
        }))
    }else contrasts <- NULL
    list(score=score,contrasts=contrasts)
}



riskQuantile.survival <- function(DT,N,NT,NF,dolist,Q,...){
    model=event=X=reference=status=times=cause=risk=Wt=WTi=cuminc=riskRegression_ID=NULL
    models <- unique(DT[,model])
    if (missing(Q)) Q <- c(0.05,0.25,0.5,0.75,0.95) else Q <- sort(Q) ##
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
    surv <- DT[model==models[[1]],data.table::data.table("surv"=(1/N*sum((time>times)/Wt))),by=times]
    surv[,cuminc:=1-surv]
    getQ.event <- function(Q,tp,X,time,status,WTi,surv){
        uX <- sort(unique(X[time<=tp & status==1]))
        Wx <- sapply(uX,function(x){sum((X<=x & time<=tp & status==1)/WTi)})/(N*surv[times==tp,cuminc])
        qRisk <- getQuantile(x=uX,Fx=Wx,Q=Q)
        qR <- data.table(t(qRisk))
        qR[,cause:="event"]
        qR
    } 
    getQ.eventFree <- function(Q,tp,X,time,Wt,surv){
        uX <- sort(unique(X[time>tp]))
        Wx <- sapply(uX,function(x){sum((X<=x & time>tp)/Wt)})/(N*surv[times==tp,surv])
        qRisk <- getQuantile(x=uX,Fx=Wx,Q=Q)
        qR <- data.table(t(qRisk))
        qR[,cause:="event-free"]
        qR
    }
    ## a <- DT[model==1]
    ## system.time(a[,getQ.eventFree(Q=Q,tp=times[1],X=risk,time=time,Wt=Wt,surv=surv)])
    score.eventfree <- DT[,getQ.eventFree(Q=Q,tp=times,X=risk,time=time,Wt=Wt,surv=surv),by=list(model,times)]
    ## setkey(DT,model,times)
    ## save(surv,file="~/tmp/surv.rda")
    ## save(DT,file="~/tmp/DT.rda")
    score.event <- DT[,getQ.event(Q=Q,tp=times,X=risk,time=time,status=status,WTi=WTi,surv=surv),by=list(model,times)]
    score.overall <- DT[,data.table(t(quantile(risk,probs=Q))),by=list(model,times)]
    score.overall[,cause:="overall"]
    colnames(score.overall) <- colnames(score.event)
    score <- rbindlist(list(score.overall,score.event,score.eventfree))
    setcolorder(score,c("model","times","cause",paste0("V",1:length(Q))))
    qnames <- paste0("Q",".",as.character(round(100*Q)))
    setnames(score,c("model","times","cause",qnames))
    if (length(dolist)>0){
        contrasts <- data.table::rbindlist(lapply(dolist,function(g){
            ## from all models g[-1], substract risk of model g[1] 
            setorder(DT,model,times,riskRegression_ID)
            DTdiff <- DT[model%in%g[-1]]
            DTref <- rep(DT[model==g[1],risk],length(unique(DTdiff$model)))
            DTdiff[,X:=risk-DTref]
            N <- NROW(DTdiff)
            Xrange <- DTdiff[,range(X)]
            Xmed <- DTdiff[,median(X)]
            changedist.eventfree <- DTdiff[,getQ.eventFree(Q=Q,tp=times,X=X,time=time,Wt=Wt,surv=surv),by=list(model,times)]
            changedist.event <- DTdiff[,getQ.event(Q=Q,tp=times,X=X,time=time,status=status,WTi=WTi,surv=surv),by=list(model,times)]
            changedist.overall <- DTdiff[,data.table(t(quantile(X,probs=Q))),by=list(model,times)]
            changedist.overall[,cause:="overall"]
            colnames(changedist.overall) <- colnames(changedist.event)
            changedist <- rbindlist(list(changedist.overall,changedist.event,changedist.eventfree))
            setcolorder(changedist,c("model","times","cause",paste0("V",1:length(Q))))
            setnames(changedist,c("model","times","cause",qnames))
            changedist[,reference:=g[1]]
            data.table::setcolorder(changedist,c("reference",colnames(changedist)[-length(colnames(changedist))]))
            changedist
        }))
    }
    else
        contrasts <- NULL 
    list(score=score,contrasts=contrasts)
}

riskQuantile.competing.risks <- function(DT,N,NT,NF,dolist,cause,states,Q,...){
    model=event=reference=risk=X=status=times=Wt=WTi=cause=cuminc=riskRegression_ID=NULL
    models <- unique(DT[,model])
    if (missing(Q)) Q <- c(0.05,0.25,0.5,0.75,0.95)
    else Q <- sort(Q)
    ##
    ## retrospectively (looking back from time t)
    ## we compute conditional on outcome the quantiles of predicted risks 
    ## and quantiles of changes of predicted risks
    ##
    ## if (missing(states))
    ## states <- DT[model==models[[1]] & status!=0,sort(unique(event))]
    states.code <- 1:length(states)
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
    surv <- DT[model==models[[1]],data.table::data.table("surv"=(1/N*sum((time>times)/Wt))),by=times]
    cuminc <- lapply(states.code,function(cc){DT[model==models[[1]],data.table::data.table("cuminc"=1/N*sum((event==cc & time<=times)/WTi)),by=times]})
    names(cuminc) <- states
    getQ.states <- function(Q,tp,X,time,event,WTi,cuminc,states.code){
        uX <- sort(unique(X))
        rbindlist(lapply(states.code,function(cause){
            # Note that event==cause implies status==1
            ## the variable event has values 1,2,..,k,k+1 where k+1 is censored
            Wx <- sapply(uX,function(x){sum((X<=x & event==cause & time<=tp)/WTi)})/(N*cuminc[[cause]][times==tp,cuminc])
            qRisk <- getQuantile(x=uX,Fx=Wx,Q=Q)
            qR <- data.table(t(qRisk))
            qR[,cause:=cause]
            qR
        }))
    }
    getQ.eventFree <- function(Q,tp,X,Xmed,time,Wt,surv){
        uX <- sort(unique(X))
        Wx <- sapply(uX,function(x){sum((X<=x & time>tp)/Wt)})/(N*surv[times==tp,surv])
        qRisk <- getQuantile(x=uX,Fx=Wx,Q=Q)
        qR <- data.table(t(qRisk))
        qR[,cause:="event-free"]
        qR
    }
    score.eventfree <- DT[,getQ.eventFree(Q=Q,tp=times,X=risk,time=time,Wt=Wt,surv=surv),by=list(model,times)]
    score.states <- DT[,getQ.states(Q=Q,tp=times,X=risk,time=time,event=event,WTi=WTi,cuminc=cuminc,states.code=states.code),by=list(model,times)]
    score.overall <- DT[,data.table(t(quantile(risk,probs=Q))),by=list(model,times)]
    score.overall[,cause:="overall"]
    colnames(score.overall) <- colnames(score.states)
    score <- rbindlist(list(score.overall,score.states,score.eventfree))
    setcolorder(score,c("model","times","cause",paste0("V",1:length(Q))))
    qnames <- paste0("Q",".",as.character(round(100*Q)))
    setnames(score,c("model","times","cause",qnames))
    if (length(dolist)>0){
        contrasts <- data.table::rbindlist(lapply(dolist,function(g){
            ## from all models g[-1], substract risk of model g[1] 
            setorder(DT,model,times,riskRegression_ID)
            DTdiff <- DT[model%in%g[-1]]
            DTref <- rep(DT[model==g[1],risk],length(unique(DTdiff$model)))
            DTdiff[,X:=risk-DTref]
            N <- NROW(DTdiff)
            Xrange <- DTdiff[,range(X)]
            Xmed <- DTdiff[,median(X)]
            changedist.eventfree <- DTdiff[,getQ.eventFree(Q=Q,tp=times,X=X,time=time,Wt=Wt,surv=surv),by=list(model,times)]
            changedist.states <- DTdiff[,getQ.states(Q=Q,tp=times,X=X,time=time,event=event,WTi=WTi,cuminc=cuminc,states.code=states.code),by=list(model,times)]
            changedist.overall <- DTdiff[,data.table(t(quantile(X,probs=Q))),by=list(model,times)]
            changedist.overall[,cause:="overall"]
            colnames(changedist.overall) <- colnames(changedist.states)
            changedist <- rbindlist(list(changedist.overall,changedist.states,changedist.eventfree))
            setcolorder(changedist,c("model","times","cause",paste0("V",1:length(Q))))
            setnames(changedist,c("model","times","cause",qnames))
            changedist[,reference:=g[1]]
            data.table::setcolorder(changedist,c("reference",colnames(changedist)[-length(colnames(changedist))]))
            changedist
        }))
    }else contrasts <- NULL
    if (!is.null(score))
        score[,cause:=as.character(factor(cause,levels=c("overall",states.code,"event-free"),labels=c("overall",states,"event-free")))]
    if (!is.null(contrasts))
        contrasts[,cause:=as.character(factor(cause,levels=c("overall",states.code,"event-free"),labels=c("overall",states,"event-free")))]
    list(score=score,contrasts=contrasts)
}




#----------------------------------------------------------------------
### riskQuantile.R ends here
