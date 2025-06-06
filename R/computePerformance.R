### computePerformance.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Feb 27 2022 (09:12)
## Version:
## Last-Updated: Dec  6 2024 (12:58) 
##           By: Thomas Alexander Gerds
##     Update #: 26
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
computePerformance <- function(DT,
                               N,
                               NT=NT,
                               NF=NF,
                               models,
                               response.type,
                               times,
                               jack,
                               cens.type,
                               cause,
                               states,
                               alpha,
                               se.fit,
                               conservative,
                               cens.model,
                               keep.residuals,
                               keep.vcov,
                               keep.iid,
                               dolist,
                               probs,
                               metrics,
                               plots,
                               summary,
                               ibs,
                               ipa,
                               ROC=FALSE,
                               MC,
                               IC.data,
                               breaks=NULL,
                               cutpoints=NULL){
    # {{{ input
    IPA=IBS=Brier= model = reference = NULL
    # inherit everything else from parent frame: summary, metrics, plots, alpha, probs, dolist, et
    out <- vector(mode="list",length=length(c(summary,metrics,plots)))
    names(out) <- c(summary,metrics,plots)
    input <- list(DT=DT,
                  N=N,
                  NT=NT,
                  NF=NF,
                  alpha=alpha,
                  se.fit=se.fit,
                  conservative=conservative,
                  cens.model=cens.model,
                  keep.residuals=keep.residuals,
                  keep.vcov=keep.vcov,
                  keep.iid=keep.iid,                  
                  dolist=dolist,Q=probs,ROC=FALSE,MC=MC,IC.data=IC.data,breaks=breaks,cutpoints=cutpoints)
    if (response.type=="competing.risks") {
        input <- c(input,list(cause=cause,states=states))
    }
    # }}}
    # {{{ collect data for summary statistics
    for (s in summary){
        if (s=="risks") {
            out[[s]] <- list(score=copy(input$DT)[,model:=factor(model,levels=models$levels,models$labels)],
                             contrasts=NULL)
        } else{
            out[[s]] <- do.call(paste(s,response.type,sep="."),input)
            if (NROW(out[[s]]$score)>0){
                out[[s]]$score[,model:=factor(model,levels=models$levels,labels=models$labels)]
            }
            if (NROW(out[[s]]$contrasts)>0){
                out[[s]]$contrasts[,model:=factor(model,levels=models$levels,labels=models$labels)]
                out[[s]]$contrasts[,reference:=factor(reference,levels=models$levels,labels=models$labels)]
            }
        }
    }
    # }}}
    # {{{ collect data for calibration plots
    if ("Calibration" %in% plots){
        out[["Calibration"]]$plotframe <- DT[model!=0]
        out[["Calibration"]]$plotframe[,model:=factor(model,levels=models$levels,labels=models$labels)]
        if (length(jack)>0)
            out[["Calibration"]]$plotframe <- merge(jack,DT[model!=0],by=c("riskRegression_ID","times"))
    }
    # }}}
    ## make sure that Brier score comes first, so that we can remove the null.model afterwards
    # {{{ calculating the other metrics
    for (m in sort(metrics,decreasing=TRUE)){
        if (m=="AUC" && ("ROC" %in% plots)){
            input <- replace(input, "ROC",TRUE)
            ## call AUC method
            out[[m]] <- do.call(paste(m,response.type,sep="."),input)
            out[["ROC"]]$plotframe <- out[[m]]$ROC
            out[["ROC"]]$plotframe[,model:=factor(model,levels=models$levels,labels=models$labels)]
            out[[m]]$ROC <- NULL
        }else{
            input <- replace(input, "ROC",FALSE)
            ## call Brier or AUC method
            out[[m]] <- do.call(paste(m,response.type,sep="."),input)
        }
        if (!is.null(out[[m]]$score)){
            out[[m]]$score[,model:=factor(model,levels=models$levels,labels=models$labels)]
        }
        ## set model and reference in model comparison results
        if (!is.null(out[[m]]$contrasts)>0){
            out[[m]]$contrasts[,model:=factor(model,levels=models$levels,labels=models$labels)]
            out[[m]]$contrasts[,reference:=factor(reference,levels=models$levels,labels=models$labels)]
        }
    }
    ## summary should be after metrics because IBS and IPA/R^2 depends on Brier score
    if (ibs){
        Dint <- function(x,y,range,na.omit=FALSE){
            if (is.null(range)) range=c(x[1],x[length(x)])
            ##   integrate a step function f with
            ##   values y=f(x) between range[1] and range[2]
            start <- max(range[1],min(x))
            Stop <- min(range[2],max(x))
            if ((Stop-start)<=0)
                return(0)
            else{
                Y=y[x>=start & x<Stop]
                X=x[x>=start & x<Stop]
                if (na.omit){
                    X=X[!is.na(Y)]
                    Y=Y[!is.na(Y)]
                } else if (any(is.na(Y))|| any(is.na(X))){
                    return(NA)
                }
                return(1/(Stop-start) * sum(Y*diff(c(X,Stop))))
            }
        }
        if (response.type!="binary"){
            out[["Brier"]][["score"]][,IBS:=sapply(times,function(t){
                Dint(x=c(0,times),y=c(0,Brier),range=c(0,t))
            }),by=c("model")]
        }
    }
    if (ipa == TRUE){
        if (response.type=="binary")
            out[["Brier"]][["score"]][,IPA:=1-Brier/Brier[model=="Null model"]]
        else
            out[["Brier"]][["score"]][,IPA:=1-Brier/Brier[model=="Null model"],by=times]
    }
    # }}}
    out[]
}


#----------------------------------------------------------------------
### computePerformance.R ends here
