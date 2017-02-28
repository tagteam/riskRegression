### plotROC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun 23 2016 (10:27) 
## Version: 
## last-updated: feb 28 2017 (13:42) 
##           By: Brice Ozenne
##     Update #: 24
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
##' @param x Object obtained with \code{Score.list}
##' @param models Choice of models to plot
##' @param times Time point specifying the prediction horizon
##' @param xlab Label for x-axis
##' @param ylab Label for y-axis
##' @param lwd line width
##' @param legend Logical. If \code{TRUE} draw legend.
##' @param legend.title Legend title
##' @param ... Not yet used
##' @examples
##' ## binary
##' set.seed(18)
##' library(randomForest)
##' library(survival)
##' bdl <- sampleData(40,outcome="binary")
##' bdt <- sampleData(58,outcome="binary")
##' bdl[,y:=factor(Y)]
##' bdt[,y:=factor(Y)]
##' fb1 <- glm(y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=bdl,family="binomial")
##' fb2 <- randomForest(y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=bdl)
##' xb <- Score(list("glm"=fb1,"rf"=fb2),y~1,data=bdt,
##'             plots="roc",metrics="auc")
##' plotROC(xb)
##' ## survival
##' set.seed(18)
##' sdl <- sampleData(40,outcome="survival")
##' sdt <- sampleData(58,outcome="survival")
##' fs1 <- coxph(Surv(time,event)~X3+X5+X6+X7+X8+X10,data=sdl,x=TRUE)
##' fs2 <- coxph(Surv(time,event)~X1+X2+X9,data=sdl,x=TRUE)
##' xs <- Score(list(model1=fs1,model2=fs2),Hist(time,event)~1,data=sdt,
##'             times=5,plots="roc",metrics="auc")
##' plotROC(xs)
##' ## competing risks
##' data(Melanoma)
##' f1 <- CSC(Hist(time,status)~age+sex+epicel+ulcer,data=Melanoma)
##' f2 <- CSC(Hist(time,status)~age+sex+logthick+epicel+ulcer,data=Melanoma)
##' x <- Score(list(model1=f1,model2=f2),Hist(time,status)~1,data=Melanoma,
##'             cause=1,times=5*365.25,plots="roc",metrics="auc")
##' plotROC(x)
#' @export
plotROC <- function(x,
                    models,
                    times,
                    xlab="1-Specificity",
                    ylab="Sensitivity",
                    lwd=3,
                    legend=TRUE,
                    legend.title,
                    ...){
    model=FPR=TPR=times=NULL
    pframe <- x$ROC$plotframe
    if (!missing(models)){
        pframe <- pframe[model%in%models]
    }else{
        if (length(x$nullModel)>0){
            pframe <- pframe[model!=x$nullModel]
        }
    }
    setkey(pframe,model)
    if (x$responseType!="binary"){
        if (missing(times))
            times <- max(pframe[["times"]])
        else ## can only do one times
            times <- times[[1]]
        pframe <- pframe[times==times]
    }
    plot(0,0,type="n",ylim = 0:1,xlim = 0:1,axes=FALSE,xlab = xlab,ylab = ylab)
    prodlim::PercentAxis(1,at=seq(0,1,.25))
    prodlim::PercentAxis(2,at=seq(0,1,.25))
    if (!missing(models)) pframe <- pframe[model %in% models]
    pframe[,col:=as.numeric(as.factor(model))]
    pframe[,lwd:=lwd]
    pframe[,lines(c(0,FPR,1),c(0,TPR,1),type="l",lwd=lwd,col=col),by=model]
    if (legend!=FALSE){
        if (!is.character(legend)){
            if (x$responseType=="binary"){
                if (!missing(models)){
                    auc <- x$AUC$score[(model%in%models)]
                } else{
                    auc <- x$AUC$score
                    if (length(x$nullModel)>0){
                        auc <- auc[model!=x$nullModel]
                    }
                }
            }else{
                if (!missing(models)){
                    auc <- x$AUC$score[(times==times) & (model%in%models)]
                } else{
                    auc <- x$AUC$score[(times==times)]
                    if (length(x$nullModel)>0){
                        auc <- auc[model!=x$nullModel]
                    }
                }
            }
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
