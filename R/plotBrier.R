### plotBrier.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 23 2017 (11:07) 
## Version: 
## last-updated: Feb 28 2017 (20:21) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Plot Brier score curves
##'
##' Plot Brier score curves
##' @title Plot Brier curve
##' @param x Object obtained with \code{Score.list}
##' @param models Choice of models to plot
##' @param lwd Line width
##' @param xlim Limits for x-axis 
##' @param ylim Limits for y-axis 
##' @param axes Logical. If \code{TRUE} draw axes.
##' @param ... Not yet used
##' @examples
##' survival
##' ds1=sampleData(40,outcome="survival")
##' ds2=sampleData(40,outcome="survival")
##' f1 <- coxph(Surv(time,event)~X1+X3+X5+X7+X9,data=ds1,x=TRUE)
##' f2 <- coxph(Surv(time,event)~X2+X4+6+X8+X10,data=ds1,x=TRUE)
##' xscore <- Score(list(f1,f2),formula=Hist(time,event)~1,data=ds2,metrics="brier")
#' @export 
plotBrier <- function(x,models,lwd=3,xlim,ylim,axes=TRUE,...){
    times=model=Brier=dimcol=lower.Brier=upper.Brier=NULL
    pframe <- x$Brier$score
    if (missing(xlim)) xlim <- pframe[,range(times)]
    if (missing(ylim)) ylim <- c(0,.25)
    plot(0,0,type="n",ylim = ylim,
         xlim = xlim,
         axes=FALSE,
         xlab = "Time",
         ylab = "Brier score")
    if (axes){
        axis(1)
        prodlim::PercentAxis(2,at=seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5))
    }
    if (!missing(models)) pframe <- pframe[model %in% models]
    pframe[,col:=as.numeric(as.factor(model))]
    pframe[,lwd:=lwd]
    pframe[,lines(times,Brier,type="l",lwd=lwd,col=col),by=model]
    pframe[,dimcol:=prodlim::dimColor(col[[1]],density=55),by=model]
    pframe[,polygon(x=c(times,rev(times)),y=c(lower.Brier,rev(upper.Brier)),col=dimcol,border=NA),by=model]
}


#----------------------------------------------------------------------
### plotBrier.R ends here
