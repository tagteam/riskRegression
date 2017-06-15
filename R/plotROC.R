### plotROC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun 23 2016 (10:27) 
## Version: 
## last-updated: Jun 10 2017 (17:45) 
##           By: Thomas Alexander Gerds
##     Update #: 60
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
##' @title Plot ROC curves
#' @param x Object obtained with function \code{Score}
#' @param models  Choice of models to plot
#' @param times A single time point specifying the prediction horizon
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param col line color
#' @param lwd line width
#' @param lty line style
#' @param cex point size
#' @param pch point style
#' @param legend logical. If \code{1L} draw a legend with the values of AUC.
#' @param add logical. If \code{1L} add lines to an existing plot.
#' @param ... Used for additional control of the subroutines: plot,
#'     axis, lines, legend. See \code{\link{SmartControl}}.
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
                    col,
                    lwd=3,
                    lty=1,
                    cex=1,
                    pch=1,
                    legend=TRUE,
                    add=FALSE,
                    ...){
    if (is.null(x$ROC))
        stop("Object has no information for ROC curves.\nYou should call the function \"riskRegression::Score\" with plots=\"ROC\".")
    model=FPR=TPR=times=NULL
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    pframe <- x$ROC$plotframe
    if (!missing(models)){
        pframe <- pframe[model%in%models]
    }else{
        if (length(x$nullModel)>0){
            pframe <- pframe[model!=x$nullModel]
        }
    }
    setkey(pframe,model)
    mm <- unique(pframe$model)
    lenmm <- length(mm)
    if(missing(col)) col <- rep(cbbPalette,length.out=lenmm)
    names(col) <- mm
    if(missing(lwd)) lwd <- 2
    lwd <- rep(lwd,length.out=lenmm)
    names(lwd) <- mm
    pch <- rep(pch,length.out=lenmm)
    names(pch) <- mm
    if(missing(lwd)) lty <- 1
    lty <- rep(lty,length.out=lenmm)
    names(lty) <- mm
    lines.DefaultArgs <- list(pch=pch,cex=cex,lwd=lwd,type="l",col=col,lty=lty)
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
        }
        auc.legend <- paste0(auc$model,": ",sprintf(fmt="%s [%s;%s]",
                                                    round(100*auc$AUC,digits=1),
                                                    round(100*auc$lower.AUC,digits=1),
                                                    round(100*auc$upper.AUC,digits=1)))
    }else auc.legend <- FALSE
    legend.DefaultArgs <- list(legend=auc.legend,lwd=lwd,col=col,lty=lty,cex=cex,bty="n",y.intersp=1.3,x="bottomright",title="AUC")
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = c(0,1),xlim = c(0,1),ylab=ylab,xlab=xlab)
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(0,1,.25))
    axis2.DefaultArgs <- list(side=2,las=1,at=seq(0,1,.25))
    control <- prodlim::SmartControl(call= list(...),
                                     keys=c("plot","lines","legend","axis1","axis2"),
                                     ignore=NULL,
                                     ignore.case=TRUE,
                                     defaults=list("plot"=plot.DefaultArgs,
                                                   "lines"=lines.DefaultArgs,
                                                   "legend"=legend.DefaultArgs,
                                                   "axis1"=axis1.DefaultArgs,
                                                   "axis2"=axis2.DefaultArgs),
                                     forced=list("plot"=list(axes=FALSE),
                                                 "axis1"=list(side=1),"axis2"=list(side=2)),
                                     verbose=TRUE)
    if (x$responseType!="binary"){
        if (missing(times))
            times <- max(pframe[["times"]])
        else ## can only do one times
            times <- times[[1]]
        pframe <- pframe[times==times]
    }
    if (add==0L) do.call("plot",control$plot)
    ## plot(0,0,type="n",ylim = 0:1,xlim = 0:1,axes=FALSE,xlab = xlab,ylab = ylab)
    control$axis1$labels <- paste(100*control$axis1$at,"%")
    control$axis2$labels <- paste(100*control$axis2$at,"%")
    do.call("axis",control$axis1)
    do.call("axis",control$axis2)
    if (!(is.logical(legend[1]) && legend[1]==FALSE)){
        do.call("legend",control$legend)
    }
    ## pframe[,col:=as.numeric(as.factor(model))]
    ## pframe[,lwd:=lwd]
    ## pframe[,lines(c(0,FPR,1),c(0,TPR,1),type="l",lwd=lwd,col=col),by=model]
    abline(a=0,b=1,col="gray77",lwd=3)
    pframe[,{thisline <- control$line;
        thisline$col=thisline$col[[as.character(model[1])]];
        thisline$lwd=thisline$lwd[[as.character(model[1])]];
        thisline$lty=thisline$lty[[as.character(model[1])]];
        thisline$pch=thisline$pch[[as.character(model[1])]];
        ## thisline$type=thisline$type[[as.character(model[1])]];
        thisline$x=c(0,FPR,1);
        thisline$y=c(0,TPR,1);
        do.call("lines",thisline)},by=model]
    invisible(x)
}

#----------------------------------------------------------------------
### plotROC.R ends here
