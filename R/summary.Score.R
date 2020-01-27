### summary.Score.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Dec 26 2019 (08:58) 
## Version: 
## Last-Updated: Jan 27 2020 (07:37) 
##           By: Thomas Alexander Gerds
##     Update #: 27
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
##' @param digits For rounding everything but p-values
##' @param pvalue.digits For rounding p-values
##' @param ... not used
##' @return List of tables
##' @seealso Score
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
summary.Score <- function(object,times,what=c("score","contrasts"),digits=1,pvalue.digits=4,...){
    ## score
    out=AUC=model=lower=upper=Brier=Model=x=reference=delta.AUC=p=delta.Brier=Reference=NULL
    if (object$response.type=="binary"){
        if ("score"%in% tolower(what)){
            if ("upper"%in%names(object$AUC$score)){
                tab.AUC <- object$AUC$score[,data.table::data.table(Model=c("Null model",as.character(model)),
                                                                    "AUC (%)"=c(Publish::pubformat(50.00,digits=digits),
                                                                                Publish::formatCI(x=100*AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits)))]
            }else{
                tab.AUC <- object$AUC$score[,data.table::data.table(Model=c("Null model",as.character(model)),
                                                                    "AUC (%)"=c(Publish::pubformat(50.00,digits=digits),
                                                                                Publish::pubformat(x=100*AUC,digits=digits)))]
            }
            if ("upper"%in%names(object$Brier$score)){
                if ("IPA"%in% names(object$Brier$score))
                    tab.Brier <- object$Brier$score[,data.table::data.table(Model=model,"Brier (%)"=Publish::formatCI(x=100*Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits),IPA=sprintf(paste0("%1.",digits,"f"),100*IPA))]
                else
                    tab.Brier <- object$Brier$score[,data.table::data.table(Model=model,"Brier (%)"=Publish::formatCI(x=100*Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
            }else{
                if ("IPA"%in% names(object$Brier$score))
                    tab.Brier <- object$Brier$score[,data.table::data.table(Model=model,"Brier (%)"=Publish::pubformat(x=100*Brier,digits=digits),IPA=sprintf(paste0("%1.",digits,"f"),100*IPA))]
                else
                    tab.Brier <- object$Brier$score[,data.table::data.table(Model=model,"Brier (%)"=Publish::pubformat(x=100*Brier,digits=digits))]
            }
            setkey(tab.Brier,Model)
            setkey(tab.AUC,Model)
            tab <- tab.AUC[tab.Brier]
            out <- list(score=tab)
            ## print(tab[])
        }
        if ("contrasts"%in% tolower(what)){
            if ("upper"%in%names(object$AUC$contrasts)){
                if ("p"%in%names(object$AUC$contrasts))
                    tab.deltaAUC <- x$AUC$contrasts[,data.table::data.table(Model=model,Reference=reference,"delta AUC (%)"=Publish::formatCI(x=100*delta.AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits),"p-value"=format.pval(p,eps=0.001,digits=pvalue.digits))]
                else
                    tab.deltaAUC <- x$AUC$contrasts[,data.table::data.table(Model=model,Reference=reference,"delta AUC (%)"=Publish::formatCI(x=100*delta.AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
            }else{
                tab.deltaAUC <- x$AUC$contrasts[,data.table::data.table(Model=model,Reference=reference,"delta AUC (%)"=Publish::pubformat(x=100*delta.AUC,digits=pvalue.digits))]
            }
            if ("upper"%in%names(object$Brier$contrasts)){
                if ("p"%in%names(object$Brier$contrasts))
                    tab.deltaBrier <- x$Brier$contrast[reference!="Null model",data.table::data.table(Model=model,Reference=reference,"delta Brier (%)"=Publish::formatCI(x=100*delta.Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits),"p-value"=format.pval(p,eps=0.001,digits=pvalue.digits))]
                else
                    tab.deltaBrier <- x$Brier$contrast[reference!="Null model",data.table::data.table(Model=model,Reference=reference,"delta Brier (%)"=Publish::formatCI(x=100*delta.Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
            }else{
                tab.deltaBrier <- x$Brier$contrast[reference!="Null model",data.table::data.table(Model=model,Reference=reference,"delta Brier (%)"=Publish::pubformat(x=100*delta.Brier,digits=digits))]
            }
            setkey(tab.deltaBrier,Model,Reference)
            setkey(tab.deltaAUC,Model,Reference)
            tab.delta <- tab.deltaAUC[tab.deltaBrier]
            if ("i.p-value"%in%names(tab.delta))
                setnames(tab.delta,"i.p-value","p-value")
            ## tab.delta[1:2,3:4] <- ""
            out <- c(out,list(contrast=tab.delta))
            ## print(tab.delta[])
        }
    }else{
        if (missing(times)) ttt = object$times else ttt <- times
        if ("score"%in% tolower(what)){
            setkey(object$AUC$score,times)
            tab.NULLAUC <- data.table(times=ttt,
                                      Model=rep("Null model",length(ttt)),
                                      "AUC (%)"=rep(Publish::pubformat(50,digits=digits),length(ttt)))
            if ("upper"%in%names(object$AUC$score)){
                tab.AUC <- rbindlist(list(tab.NULLAUC,object$AUC$score[times%in%ttt,data.table::data.table(times=times,Model=as.character(model),"AUC (%)"=Publish::formatCI(x=100*AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]))
            }else{
                tab.AUC <- rbindlist(list(tab.NULLAUC,object$AUC$score[times%in%ttt,data.table::data.table(times=times,Model=as.character(model),"AUC (%)"=Publish::pubformat(x=100*AUC,digits=digits))]))
            }
            if ("upper"%in%names(object$Brier$score)){
                if ("IPA"%in% names(object$Brier$score))
                    tab.Brier <- object$Brier$score[times%in%ttt,data.table::data.table(times=times,Model=model,Brier=Publish::formatCI(x=100*Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits),IPA=sprintf(paste0("%1.",digits,"f"),100*IPA))]
                else
                    tab.Brier <- object$Brier$score[times%in%ttt,data.table::data.table(times=times,Model=model,Brier=Publish::formatCI(x=100*Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
            }else{
                if ("IPA"%in% names(object$Brier$score))
                    tab.Brier <- object$Brier$score[times%in%ttt,data.table::data.table(times=times,Model=model,Brier=Publish::pubformat(x=100*Brier,digits=digits),IPA=sprintf(paste0("%1.",digits,"f"),100*IPA))]
                else
                    tab.Brier <- object$Brier$score[times%in%ttt,data.table::data.table(times=times,Model=model,Brier=Publish::pubformat(x=100*Brier,digits=digits))]
            }
            setnames(tab.Brier,"Brier","Brier (%)")
            setkey(tab.Brier,times,Model)
            setkey(tab.AUC,times,Model)
            tab <- tab.AUC[tab.Brier]
            out <- list(score=tab)
        }
        if ("contrasts"%in% tolower(what)){
            if ("upper"%in%names(object$AUC$contrasts)){
                if ("p"%in%names(object$AUC$contrasts))
                    tab.deltaAUC <- x$AUC$contrasts[times%in%ttt,data.table::data.table(times=times,Model=model,Reference=reference,"delta AUC (%)"=Publish::formatCI(x=100*delta.AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits),"p-value"=format.pval(p,eps=0.001,digits=pvalue.digits))]
                else
                    tab.deltaAUC <- x$AUC$contrasts[times%in%ttt,data.table::data.table(times=times,Model=model,Reference=reference,"delta AUC (%)"=Publish::formatCI(x=100*delta.AUC,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
            }else{
                tab.deltaAUC <- x$AUC$contrasts[times%in%ttt,data.table::data.table(times=times,Model=model,Reference=reference,"delta AUC (%)"=Publish::pubformat(x=100*delta.AUC,digits=digits))]
            }
            if ("upper"%in%names(object$Brier$contrasts)){
                if ("p"%in%names(object$AUC$contrasts))
                    tab.deltaBrier <- x$Brier$contrast[times%in%ttt & reference!="Null model",data.table::data.table(times=times,Model=model,Reference=reference,"delta Brier (%)"=Publish::formatCI(x=100*delta.Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits),"p-value"=format.pval(p,eps=0.001,digits=pvalue.digits))]
                else
                    tab.deltaBrier <- x$Brier$contrast[times%in%ttt & reference!="Null model",data.table::data.table(times=times,Model=model,Reference=reference,"delta Brier (%)"=Publish::formatCI(x=100*delta.Brier,lower=100*lower,upper=100*upper,show.x=1,digits=digits))]
            }else{
                tab.deltaBrier <- x$Brier$contrast[times%in%ttt & reference!="Null model",data.table::data.table(times=times,Model=model,Reference=reference,"delta Brier (%)"=Publish::pubformat(x=100*delta.Brier,digits=digits))]
            }
            setkey(tab.deltaBrier,times,Model,Reference)
            setkey(tab.deltaAUC,times,Model,Reference)
            tab.delta <- tab.deltaAUC[tab.deltaBrier]
            if ("i.p-value"%in%names(tab.delta))
                setnames(tab.delta,"i.p-value","p-value")
            out <- c(out,list(contrast=tab.delta))
        }
    }
    out
}


######################################################################
### summary.Score.R ends here
