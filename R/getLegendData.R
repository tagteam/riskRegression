### getLegendData.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Oct 13 2017 (10:50) 
## Version: 
## Last-Updated: Oct 13 2017 (13:04) 
##           By: Thomas Alexander Gerds
##     Update #: 19
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getLegendData <- function(object,
                          models,
                          brier.in.legend=TRUE,
                          format.brier,
                          auc.in.legend=TRUE,
                          format.auc,
                          drop.null.model=TRUE,
                          scale=100,
                          digits=1,
                          ...){
    model=AUC=lower=upper=Brier=NULL
    if (missing(models)) {
        models <- names(object$models)
    }
    if (drop.null.model==TRUE) models <- models[models!=object$null.model]
    maxlen <- max(nchar(models))
    legend.text.models <- sprintf(paste0("%",maxlen,"s"),models)
    if (missing(format.auc))
        format.auc <- paste0("%1.",digits,"f [%1.",digits,"f;%1.",digits,"f]")
    if (missing(format.brier))
        format.brier <- paste0("%1.",digits,"f [%1.",digits,"f;%1.",digits,"f]")
    if (is.null(object$null.model) || drop.null.model==TRUE){
        keep.null.model <- FALSE
    }else{
        if (brier.in.legend==TRUE){
            keep.null.model=TRUE
        }else{
            keep.null.model=match(object$null.model,models,nomatch=FALSE)
        }
    }
    if (auc.in.legend==TRUE){
        if (is.null(object$AUC)){
            warning("Cannot show AUC as it is not stored in object. Set metrics='auc' in the call of Score.")
            legend.text.auc <- NULL
        }else{
            auc.data <- object$AUC$score[(model%in%models)]
            if (keep.null.model==FALSE){
                auc.data <- auc.data[model!=object$null.model]
            }
            ## user's order
            auc.data[,model:=factor(model,levels=models)]
            setkey(auc.data,model)
            legend.text.auc <- auc.data[,sprintf(fmt=format.auc,scale*AUC,scale*lower,scale*upper)]
        }
    }else{
        legend.text.auc <- NULL
    }
    if (brier.in.legend==TRUE){
        if (is.null(object$Brier)){
            warning("Cannot show Brier score as it is not stored in object. Set metrics='brier' in the call of Score.")
            legend.text.brier <- NULL
        }else{
            brier.data <- object$Brier$score[(model%in%models)]
            if (keep.null.model==FALSE){
                brier.data <- brier.data[model!=object$null.model]
            }
            brier.data[,model:=factor(model,levels=models)]
            setkey(brier.data,model)
            legend.text.brier <- brier.data[,sprintf(fmt=format.brier,scale*Brier,scale*lower,scale*upper)]
        }
    }else{
        legend.text.brier <- NULL
    }
    out <- data.table(cbind(model=legend.text.models,AUC=legend.text.auc,Brier=legend.text.brier))
    cnames <- names(out)
    cnames[1] <- ""
    clen <- lapply(out[1,],nchar)
    fnames <- sprintf(paste0("%",clen,"s"),cnames)
    attr(out,"format.names") <- fnames
    out
}


######################################################################
### getLegendData.R ends here
