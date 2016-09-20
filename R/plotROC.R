### plotROC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun 23 2016 (10:27) 
## Version: 
## last-updated: Aug 30 2016 (14:16) 
##           By: Thomas Alexander Gerds
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Plot ROC curve 
##'
##' @title Plot ROC curve
#' @export
##' @param x Object obtained with \code{Score.list}
##' @param models Choice of models to plot
##' @param timepoint Time point specifying the prediction horizon
##' @param lwd line width
##' @param legend Logical. If \code{TRUE} draw legend.
##' @param legend.title Legend title
##' @param ... Not yet used
plotROC <- function(x,models,timepoint,lwd=3,legend=TRUE,legend.title,...){
    model=FPR=TPR=times=NULL
    pframe <- x$ROC$plotframe
    if (!missing(models)){
        pframe <- pframe[model%in%models]
    }
    browser()
    setkey(pframe,model)
    if (missing(timepoint))
        timepoint <- max(pframe[["times"]])
    else ## can only do one timepoint
        timepoint <- timepoint[[1]]
    pframe <- pframe[times==timepoint]
    plot(0,0,type="n",ylim = 0:1,xlim = 0:1,axes=FALSE,xlab = "1-Specificity",ylab = "Sensitivity")
    prodlim::PercentAxis(1,at=seq(0,1,.25))
    prodlim::PercentAxis(2,at=seq(0,1,.25))
    if (!missing(models)) pframe <- pframe[model %in% models]
    pframe[,col:=as.numeric(as.factor(model))]
    pframe[,lwd:=lwd]
    pframe[,lines(c(0,FPR,1),c(0,TPR,1),type="l",lwd=lwd,col=col),by=model]
    if (legend!=FALSE){
        if (!is.character(legend)){
            if (!missing(models))
                auc <- x$AUC$score[(times==timepoint) & (model%in%models)]
            else
                auc <- x$AUC$score[(times==timepoint)]
            setkey(auc,model)
            ## legend <- paste0(auc$model,": ",Publish::formatCI(100*auc$AUC,100*auc$lower.AUC,100*auc$upper.AUC,showX=TRUE,digits=1))
            legend <- paste0(auc$model,": ",sprintf(fmt="%s [%s;%s]",
                                                    round(100*auc$AUC,digits=1),
                                                    round(100*auc$lower.AUC,digits=1),
                                                    round(100*auc$upper.AUC,digits=1)))
        }
        if (missing(legend.title)) legend.title <- "AUC"
        pframe[,legend(x="bottomright",legend=legend,lwd=lwd[[1]],title=legend.title,col=unique(col))]
    }
    abline(a=0,b=1,col="gray77",lwd=3)
    invisible(x)
}

#----------------------------------------------------------------------
### plotROC.R ends here
