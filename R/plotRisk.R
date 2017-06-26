### plotRisk.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar 13 2017 (16:53) 
## Version: 
## Last-Updated: Jun 25 2017 (12:31) 
##           By: Thomas Alexander Gerds
##     Update #: 33
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' plot predicted risks 
##'
##' Two rival prediction models are applied to the same data.
##' @title plot predicted risks
##' @param x Object obtained with function \code{Score}
##' @param models Choice of two models to plot. The predicted risks of
##'     the first (second) are shown along the x-axis (y-axis).
##' @param times Time point specifying the prediction horizon.
##' @param xlim x-axis limits
##' @param ylim y-axis limits
##' @param xlab x-axis labels
##' @param ylab y-axis labels
##' @param col colour
##' @param pch point type
##' @param cex point size
##' @param ... Used to control the subroutines: plot, axis, lines,
##'     barplot, legend. See \code{\link{SmartControl}}.
##' @return a nice graph
##' @examples
##' ## uncensored
##' learndat = sampleData(40,outcome="binary")
##' testdat = sampleData(40,outcome="binary")
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family="binomial")
##' lr2 = glm(Y~X3+X5+X6,data=learndat,family="binomial")
##' xb=Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,
##'          data=testdat,summary="risks",nullModel=0L)
##' plotRisk(xb)
##' ## survival
##' library(survival)
##' learndat = sampleData(40,outcome="survival")
##' testdat = sampleData(40,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=learndat,x=TRUE)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=learndat,x=TRUE)
##' xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),formula=Surv(time,event)~1,
##'          data=testdat,summary="risks",nullModel=0L)
##' plotRisk(xs,times=5)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
plotRisk <- function(x,
                     models,
                     times,
                     xlim=c(0,1),
                     ylim=c(0,1),
                     xlab,
                     ylab,
                     col,
                     pch=1,
                     cex=1,
                     ...){
    model=ReSpOnSe=risk=status=NULL
    if (!is.null(x$nullModel))
        pframe <- x$risks$score[model!=x$nullModel]
    else
        pframe <- x$risks$score
    if (missing(models)){
        models <- pframe[,unique(model)[1:2]]
    }
    if (is.na(models[2])) stop("Need two models to scatterplot risks")
    pframe <- pframe[model%in%models]
    modelnames <- pframe[,unique(model)]
    if (x$responseType=="binary")
        R <- pframe[model==modelnames[1],ReSpOnSe]
    else{
        warning("Under construction")
        R <- pframe[model==modelnames[1],status]
    }
    if (missing(col)) col <- R+1
    if (missing(xlab)) xlab <- paste0("Risk (%, ",modelnames[1],")")
    if (missing(ylab)) ylab <- paste0("Risk (%, ",modelnames[2],")")
    m1 <- pframe[model==modelnames[1],risk]
    m2 <- pframe[model==modelnames[2],risk]
    # {{{ smart argument control
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,ylab=ylab,xlab=xlab)
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4))
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),mgp=c(4,1,0))
    legend.DefaultArgs <- list(legend=names(x$models),pch=pch,col=col,cex=cex,bty="n",y.intersp=1.3,x="topleft")
    points.DefaultArgs <- list(x=m1,y=m2,pch=pch,cex=cex,col=col)
    abline.DefaultArgs <- list(a=0,b=1,lwd=1,col="gray66")
    control <- prodlim::SmartControl(call= list(...),
                                     keys=c("plot","points","legend","axis1","axis2","abline"),
                                     ignore=NULL,
                                     ignore.case=TRUE,
                                     defaults=list("plot"=plot.DefaultArgs,
                                                   "points"=points.DefaultArgs,
                                                   "abline"=abline.DefaultArgs,
                                                   "legend"=legend.DefaultArgs,
                                                   "axis1"=axis1.DefaultArgs,
                                                   "axis2"=axis2.DefaultArgs),
                                     forced=list("plot"=list(axes=FALSE)),
                                     verbose=TRUE)
    # }}}
    do.call("plot",control$plot)
    do.call("abline",control$abline)
    do.call("points",control$points)
    control$axis2$labels <- paste(100*control$axis2$at,"%")
    control$axis1$labels <- paste(100*control$axis1$at,"%")
    do.call("axis",control$axis1)
    do.call("axis",control$axis2)
}


######################################################################
### plotRisk.R ends here
