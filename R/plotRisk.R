### plotRisk.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar 13 2017 (16:53) 
## Version: 
## Last-Updated: Apr 17 2019 (17:56) 
##           By: Thomas Alexander Gerds
##     Update #: 73
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
##' @param preclipse Value between 0 and 1 defining the preclipse area
##' @param ... Used to control the subroutines: plot, axis, lines,
##'     barplot, legend. See \code{\link{SmartControl}}.
##' @return a nice graph
##' @examples
##' library(prodlim)
##' ## uncensored
##' learndat = sampleData(40,outcome="binary")
##' testdat = sampleData(40,outcome="binary")
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family="binomial")
##' lr2 = glm(Y~X3+X5+X6,data=learndat,family="binomial")
##' xb=Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,
##'          data=testdat,summary="risks",null.model=0L)
##' plotRisk(xb)
##' ## survival
##' library(survival)
##' learndat = sampleData(40,outcome="survival")
##' testdat = sampleData(40,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=learndat,x=TRUE)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=learndat,x=TRUE)
##' xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),formula=Surv(time,event)~1,
##'          data=testdat,summary="risks",null.model=0L)
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
                     pch=3,
                     cex=1,
                     preclipse=0,
                     ...){
    model=ReSpOnSe=risk=status=NULL
    if (is.null(x$risks$score)) stop("No predicted risks in object. You should set summary='risks' when calling Score.")
    if (!is.null(x$null.model)){
        pframe <- x$risks$score[model!=x$null.model]
        pframe[,model:=factor(model)]
    } else{
        pframe <- x$risks$score
    }
    if (x$response.type!="binary"){
        if (missing(times)){
            tp <- max(pframe[["times"]])
            if (length(unique(pframe$times))>1)
                warning("Time point not specified, use max of the available times: ",tp)
        } else{ ## can only do one time point
            tp <- times[[1]]
            if (!(tp%in%unique(pframe$times)))
                stop(paste0("Requested time ",times[[1]]," is not in object"))
        }
        pframe <- pframe[times==tp]
    }else tp <- NULL
    if (missing(models)){
        models <- pframe[,levels(model)[1:2]]
    }
    if (is.na(models[2])) stop("Need two models for a scatterplot of predicted risks")
    pframe <- pframe[model%in%models]
    modelnames <- pframe[,levels(model)]
    if (x$response.type=="binary"){
        R <- pframe[model==modelnames[1],ReSpOnSe]
    }
    else{
        R <- pframe[model==modelnames[1],
        {
            r <- status # 0,1,2?
            r[time>times] <- 2
            r
        }]
    }
    ## m1 <- pframe[model==modelnames[1],.(risk,ID)]
    ## m2 <- pframe[model==modelnames[2],.(risk,ID)]
    pframe[,model:=factor(model)]
    m1 <- levels(pframe$model)[[1]]
    m2 <- levels(pframe$model)[[2]]
    if (missing(xlab)) xlab <- paste0("Risk (%, ",modelnames[1],")")
    if (missing(ylab)) ylab <- paste0("Risk (%, ",modelnames[2],")")
    if (x$response.type=="binary"){
        ppframe <- dcast(pframe,ID~model,value.var="risk")
    }else{
        ppframe <- dcast(pframe,times+ID~model,value.var="risk")
    }
    pred.m1 <- ppframe[[m1]]
    pred.m2 <- ppframe[[m2]]
    if (preclipse>0){
        diff <- abs(pred.m1-pred.m2)
        which <- diff>quantile(diff,preclipse)
        message(paste(sprintf("\nShowing absolute differences > %1.2f",
                              100*quantile(diff,preclipse)),"%"))
        pred.m1 <- pred.m1[which]
        pred.m2 <- pred.m2[which]
        R <- R[which]
    }
    if (x$response.type=="binary"){
        if (missing(col)) col <- factor(R,levels=c(0,1),labels=c("darkgreen","red"))
    }
    else{
        if (missing(col)) col <- as.character(factor(R,levels=c(0,1,2),labels=c("darkorange","red","darkgreen")))
    }
    # {{{ smart argument control
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,ylab=ylab,xlab=xlab)
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4))
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),mgp=c(4,1,0))
    if(x$response.type=="binary") {
        this.legend <- paste0(c("No event","Event")," (n=",c(sum(R==0),sum(R==1)),")")
    }else{
        this.legend <- paste0(c("Censored",
                                "Event",
                                "No event"),
                              " (n=",
                              c(sum(R==0),
                                sum(R==1),
                                sum(R==2)),")")
    }
    legend.DefaultArgs <- list(legend=this.legend,
                               pch=pch,
                               col=if(x$response.type=="binary") c("darkgreen","red") else c("darkorange","red","darkgreen"),
                               cex=cex,
                               bty="n",
                               y.intersp=1.3,
                               x="topleft")
    points.DefaultArgs <- list(x=pred.m1,y=pred.m2,pch=pch,cex=cex,col=col)
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
    do.call("legend",control$legend)
}


######################################################################
### plotRisk.R ends here
