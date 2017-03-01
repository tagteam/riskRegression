### plotAUC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun 23 2016 (09:19) 
## Version: 
## last-updated: Mar  1 2017 (07:13) 
##           By: Thomas Alexander Gerds
##     Update #: 76
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Plot AUC curve
##'
##' @title Plot AUC curve
#' @param x Object obtained with \code{Score.list}
#' @param models Choice of models to plot
#' @param which Character. Either \code{"score"} to show AUC or
#'     \code{"contrasts"} to show differences between AUC.
#' @param xlim Limits for x-axis
#' @param ylim Limits for y-axis
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param col line color
#' @param lwd line width
#' @param lty line style
#' @param cex point size
#' @param pch point style
#' @param type line type
#' @param axes Logical. If \code{TRUE} draw axes.
#' @param percent Logical. If \code{TRUE} scale y-axis in %.
#' @param confint Logical. If \code{TRUE} draw confidence shadows.
#' @param legend Logical. If \code{TRUE} draw legend.
#' @param ... Used for additional control of the subroutines: plot,
#'     axis, lines, legend. See \code{\link{SmartControl}}.
##' @examples
##' library(survival)
##' d=sampleData(100,outcome="survival")
##' nd=sampleData(100,outcome="survival")
##' 
##' f1=coxph(Surv(time,event)~X1+X6+X8,data=d,x=TRUE,y=TRUE)
##' f2=coxph(Surv(time,event)~X2+X5+X9,data=d,x=TRUE,y=TRUE)
##' 
##' xx=Score(list("Cox X1+X6+X8"=f1,"Cox X2+X5+X9"=f2), formula=Surv(time,event)~1,
##'          data=nd, metrics="auc", nullModel=FALSE, times=seq(3:10))
##' plotAUC(xx)
##' plotAUC(xx,confint=TRUE)
##' plotAUC(xx,type="contrasts")
##' a=plotAUC(xx,type="contrasts",confint=TRUE)
##' a+theme_bw()
##' 
#' @export
plotAUC <- function(x,models,which="score",xlim,ylim,xlab,ylab,col,lwd,lty=1,cex=1,pch=1,type="l",axes=1L,percent=1L,confint=0L,legend=1L,...){
    times=contrast=model=AUC=lower.AUC=upper.AUC=lower=upper=delta=reference=NULL
    ## cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    pframe <- switch(which,"score"={x$AUC$score},"contrasts"={x$AUC$contrasts},{stop("argument 'which' has to be either 'score' for AUC or 'contrasts' for differences in AUC.")})
    if (length(pframe$times)<2) stop(paste("Need at least two time points for plotting time-dependent AUC. Object has only ",length(pframe$times),"times"))
    mm <- unique(pframe$model)
    lenmm <- length(mm)
    if(missing(xlab)) xlab <- "Time"
    if(missing(ylab)) ylab <- "AUC"
    if(missing(col)) col <- rep(cbbPalette,length.out=lenmm)
    names(col) <- mm
    if(missing(lwd)) lwd <- 2
    lwd <- rep(lwd,length.out=lenmm)
    names(lwd) <- mm
    if(missing(lwd)) lty <- 1
    lty <- rep(lty,length.out=lenmm)
    names(lty) <- mm
    if (missing(xlim)) xlim <- pframe[,range(times)]
    if (missing(ylim))
        if (which=="score") ylim <- c(0.5,1)
        else ylim <- c(min(pframe$lower),max(pframe$upper))
    lines.DefaultArgs <- list(pch=pch,type=type,cex=cex,lwd=lwd,col=col,lty=lty)
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(0,xlim[2],xlim[2]/4))
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(0,ylim[2],ylim[2]/4),mgp=c(4,1,0))
    legend.DefaultArgs <- list(legend=mm,lwd=lwd,col=col,lty=lty,cex=cex,bty="n",y.intersp=1.3,x="topleft")
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,ylab=ylab,xlab=xlab)
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
                                                 "axis1"=list(side=1)),
                                     verbose=TRUE)
    
    if (which=="score"){
        ## AUC
        if (!missing(models)) pframe <- pframe[model %in% models]
        yticks <- seq(0,1,0.05)
        yticks <- yticks[yticks>=ylim[1] & yticks<=ylim[2]]
        do.call("plot",control$plot)
        pframe[,lines(times,AUC,col=col[model],lwd=lwd[model],lty=lty[model]),by=model]
    }else{
        ## delta AUC
        do.call("plot",control$plot)
        pframe[,contrast:=paste(model,reference,sep=" - ")]
        pframe[,lines(times,AUC,col=col[model],lwd[model],lty[model]),by=contrast]
    }
    ## legend
    if (!(is.logical(legend[1]) && legend[1]==FALSE)){
        do.call("legend",control$legend)
    }
    ## x-axis
    if (confint==TRUE){
        dimcol <- sapply(col,function(cc){prodlim::dimColor(cc)})
        names(dimcol) <- names(col)
        if (which=="score"){
            pframe[,polygon(x=c(times,rev(times)),y=c(lower.AUC,rev(upper.AUC)),col=dimcol[model],border=NA),by=model]
        }else{
            pframe[,polygon(x=c(times,rev(times)),y=c(lower,rev(upper)),col=dimcol[model],border=NA),by=model]
        }
    }
    if (axes){
        control$axis2$labels <- paste(100*control$axis2$at,"%")
        do.call("axis",control$axis1)
        do.call("axis",control$axis2)
    }
    invisible(pframe)
}

#----------------------------------------------------------------------
### plotAUC.R ends here
