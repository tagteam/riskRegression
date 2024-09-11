### summary.Score.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Dec 26 2019 (08:58) 
## Version: 
## Last-Updated: Jul  2 2024 (11:50) 
##           By: Thomas Alexander Gerds
##     Update #: 91
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Summarizing a Score object
##'
##' The AUC and the Brier score are put into tables
##' @title Summary of prediction performance metrics
##' @param object Object obtained with \code{Score}.
##' @param times Select time points 
##' @param what Either \code{"score"}, \code{"contrasts"} or both, i.e., \code{c("score","contrasts")}
##' @param models Select which models to summarize. Need to be a subset of \code{object$models}
##' @param digits For rounding everything but p-values
##' @param pvalue.digits For rounding p-values
##' @param ... not used
##' @return List of tables
##' @seealso Score
##' @export summary.Score
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
summary.Score <- function(object,
                          times,
                          what=c("score","contrasts"),
                          models,
                          digits=1,
                          pvalue.digits=4,
                          ...){
    ## score
    out=AUC=model=lower=upper=Brier=Model=x=reference=delta.AUC=p=delta.Brier=Reference=NULL
    what <- sapply(tolower(what),function(w){match.arg(w,c("score","contrasts"))})
    # check models 
    fitted.models <- names(object$models)
    if (missing(models)){
        models <- fitted.models
    } else{
        if (!all(models %in% fitted.models))
            stop(paste0("\nRequested models:",paste0(models,collapse=", "),"\nFitted models:",paste0(fitted.models,collapse=", ")))
    }
    # {{{ binary
    if (object$response.type=="binary"){
        # {{{ score
        if ("score"%in% tolower(what)){
            if ("upper"%in%names(object$AUC$score)){
                tab.AUC <- object$AUC$score[(model %in% models),
                                            data.table::data.table(Model=factor(as.character(model),levels=fitted.models)[,drop=TRUE],
                                                                   "AUC (%)"=Publish::formatCI(x=100*AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
            }else{
                tab.AUC <- object$AUC$score[(model %in% models),
                                            data.table::data.table(Model=factor(as.character(model),levels=fitted.models)[,drop=TRUE],
                                                                   "AUC (%)"=Publish::pubformat(x=100*AUC,digits=digits))]
            }
            if ("upper"%in%names(object$Brier$score)){
                if ("IPA"%in% names(object$Brier$score))
                    tab.Brier <- object$Brier$score[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],"Brier (%)"=Publish::formatCI(x=100*Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits),IPA=sprintf(paste0("%1.",digits,"f"),100*IPA))]
                else
                    tab.Brier <- object$Brier$score[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],"Brier (%)"=Publish::formatCI(x=100*Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
            }else{
                if ("IPA"%in% names(object$Brier$score))
                    tab.Brier <- object$Brier$score[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],"Brier (%)"=Publish::pubformat(x=100*Brier,digits=digits),IPA=sprintf(paste0("%1.",digits,"f"),100*IPA))]
                else
                    tab.Brier <- object$Brier$score[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],"Brier (%)"=Publish::pubformat(x=100*Brier,digits=digits))]
            }
            ## merging the two tables if there are two
            if (length(tab.Brier)>0) {
                setkey(tab.Brier,Model)
                if (length(tab.AUC)>0) {
                    setkey(tab.AUC,Model)
                    tab <- tab.AUC[tab.Brier]
                }else{
                    tab <- tab.Brier
                }
            }else{
                setkey(tab.AUC,Model)
                tab <- tab.AUC
            }
            out <- list(score=tab)
        }
        # }}}
        # {{{ contrasts
        if ("contrasts"%in% tolower(what)){
            if ("upper"%in%names(object$AUC$contrasts)){
                if ("p"%in%names(object$AUC$contrasts)){
                    tab.deltaAUC <- object$AUC$contrasts[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta AUC (%)"=Publish::formatCI(x=100*delta.AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits),"p-value"=format.pval(p,eps=0.001,digits=pvalue.digits))]
                } else{
                    tab.deltaAUC <- object$AUC$contrasts[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta AUC (%)"=Publish::formatCI(x=100*delta.AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
                }
            }else{
                if (length(object$AUC$contrasts)>0){
                    tab.deltaAUC <- object$AUC$contrasts[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta AUC (%)"=Publish::pubformat(x=100*delta.AUC,digits=pvalue.digits))]
                }else{
                    tab.deltaAUC <- NULL
                }
            }
            if ("upper"%in%names(object$Brier$contrasts)){
                if ("p"%in%names(object$Brier$contrasts)){
                    tab.deltaBrier <- object$Brier$contrasts[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta Brier (%)"=Publish::formatCI(x=100*delta.Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits),"p-value"=format.pval(p,eps=0.001,digits=pvalue.digits))]
                } else{
                    tab.deltaBrier <- object$Brier$contrasts[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta Brier (%)"=Publish::formatCI(x=100*delta.Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
                }
            }else{
                tab.deltaBrier <- object$Brier$contrasts[(model %in% models),data.table::data.table(Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta Brier (%)"=Publish::pubformat(x=100*delta.Brier,digits=digits))]
            }
            if (length(tab.deltaBrier)>0) {
                setkey(tab.deltaBrier,Model,Reference)
                if (length(tab.deltaAUC)>0) {
                    setkey(tab.deltaAUC,Model,Reference)
                    tab.delta <- tab.deltaAUC[tab.deltaBrier]
                    if (any(is.na(tmp <- tab.delta[[3]]))){
                        tmp[is.na(tmp)] <- ""
                        set(tab.delta,j=3L,value=tmp)
                    }
                    if (any(is.na(tmp <- tab.delta[[4]]))){
                        tmp[is.na(tmp)] <- ""
                        set(tab.delta,j=4L,value=tmp)
                    }
                    rm(tmp)
                    if ("i.p-value"%in%names(tab.delta))
                        setnames(tab.delta,"i.p-value","p-value")
                }else{
                    tab.delta <- tab.deltaBrier
                }
            }else{
                if (length(tab.deltaAUC)>0){
                    setkey(tab.deltaAUC,Model,Reference)
                }
                tab.delta <- tab.deltaAUC
            }
            if (length(tab.delta)>0)
            out <- c(out,list(contrasts=tab.delta))
        }
        # }}}
        # }}}
    } else{
        # {{{ time-to-event
        if (missing(times)) ttt = object$times else ttt <- times
        # {{{ score
        if ("score"%in% tolower(what)){
            if (length(object$AUC$score)>0){
                setkey(object$AUC$score,times)
                tab.NULLAUC <- data.table(times=ttt,
                                          Model=rep("Null model",length(ttt)),
                                          "AUC (%)"=rep(Publish::pubformat(50,digits=digits),length(ttt)))
                if ("upper"%in%names(object$AUC$score)){
                    tab.AUC <- rbindlist(list(tab.NULLAUC,object$AUC$score[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=as.character(model),"AUC (%)"=Publish::formatCI(x=100*AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]))
                }else{
                    tab.AUC <- rbindlist(list(tab.NULLAUC,object$AUC$score[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=as.character(model),"AUC (%)"=Publish::pubformat(x=100*AUC,digits=digits))]))
                }
            }else{tab.AUC <- NULL}
            if (length(object$Brier$score)>0){
                if ("upper"%in%names(object$Brier$score)){
                    if ("IPA"%in% names(object$Brier$score))
                        tab.Brier <- object$Brier$score[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],"Brier (%)"=Publish::formatCI(x=100*Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits),IPA=sprintf(paste0("%1.",digits,"f"),100*IPA))]
                    else
                        tab.Brier <- object$Brier$score[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],"Brier (%)"=Publish::formatCI(x=100*Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
                }else{
                    if ("IPA"%in% names(object$Brier$score))
                        tab.Brier <- object$Brier$score[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],"Brier (%)"=Publish::pubformat(x=100*Brier,digits=digits),IPA=sprintf(paste0("%1.",digits,"f"),100*IPA))]
                    else
                        tab.Brier <- object$Brier$score[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],"Brier (%)"=Publish::pubformat(x=100*Brier,digits=digits))]
                }
            }else{
                tab.Brier <- NULL
            }
            if (length(tab.Brier)>0) {
                ## setnames(tab.Brier,"Brier","Brier (%)")
                setkey(tab.Brier,times,Model)
                if (length(tab.AUC)>0) {
                    setkey(tab.AUC,times,Model)
                    tab <- tab.AUC[tab.Brier]
                }else{
                    tab <- tab.Brier
                }
            }else{
                setkey(tab.AUC,times,Model)
                tab <- tab.AUC
            }
            out <- list(score=tab)
        }
        # }}}
        # {{{ contrasts
        if ("contrasts"%in% tolower(what)){
            if (length(object$AUC$contrasts)>0){
                if ("upper"%in%names(object$AUC$contrasts)){
                    if ("p"%in%names(object$AUC$contrasts))
                        tab.deltaAUC <- object$AUC$contrasts[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta AUC (%)"=Publish::formatCI(x=100*delta.AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits),"p-value"=format.pval(p,eps=0.001,digits=pvalue.digits))]
                    else
                        tab.deltaAUC <- object$AUC$contrasts[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta AUC (%)"=Publish::formatCI(x=100*delta.AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
                }else{
                    tab.deltaAUC <- object$AUC$contrasts[(times%in%ttt)&(model%in%models),data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta AUC (%)"=Publish::pubformat(x=100*delta.AUC,digits=digits))]
                }
            }else{tab.deltaAUC <- NULL}
            if (length(object$Brier$contrasts)>0){
                if ("upper"%in%names(object$Brier$contrasts)){
                    if ("p"%in%names(object$AUC$contrasts))
                        tab.deltaBrier <- object$Brier$contrasts[(times%in%ttt)&(model%in%models) & reference!="Null model",data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta Brier (%)"=Publish::formatCI(x=100*delta.Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits),"p-value"=format.pval(p,eps=0.001,digits=pvalue.digits))]
                    else
                        tab.deltaBrier <- object$Brier$contrasts[(times%in%ttt)&(model%in%models) & reference!="Null model",data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta Brier (%)"=Publish::formatCI(x=100*delta.Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
                }else{
                    tab.deltaBrier <- object$Brier$contrasts[(times%in%ttt)&(model%in%models) & reference!="Null model",data.table::data.table(times=times,Model=factor(model,levels=fitted.models)[,drop=TRUE],Reference=reference,"delta Brier (%)"=Publish::pubformat(x=100*delta.Brier,digits=digits))]
                }
            }else{tab.deltaBrier <- NULL}
            ## merging the two tables if there are two
            if (length(tab.deltaBrier)>0) {
                setkey(tab.deltaBrier,times,Model,Reference)
                if (length(tab.deltaAUC)>0) {
                    setkey(tab.deltaAUC,times,Model,Reference)
                    tab.delta <- tab.deltaAUC[tab.deltaBrier]
                    if ("i.p-value"%in%names(tab.delta))
                        setnames(tab.delta,"i.p-value","p-value")
                }else{
                    tab.delta <- tab.deltaBrier
                }
            }else{
                if (length(tab.deltaAUC)>0){
                    setkey(tab.deltaAUC,Model,Reference)
                }
                tab.delta <- tab.deltaAUC
            }
            if (length(tab.delta)>0)
                out <- c(out,list(contrasts=tab.delta))
        }
    }
    # }}}
    # }}}
    out
}


######################################################################
### summary.Score.R ends here
