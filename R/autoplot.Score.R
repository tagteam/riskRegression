### plotAUC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun 23 2016 (09:19) 
## Version: 
## last-updated: Jan 13 2020 (12:06) 
##           By: Thomas Alexander Gerds
##     Update #: 77
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' ggplot AUC curves
##'
##' @title ggplot AUC curve
##' @param object Object obtained with \code{Score.list}
##' @param models Choice of models to plot
##' @param type Character. Either \code{"score"} to show AUC or \code{"contrasts"} to show differences between AUC.
##' @param lwd Line width
##' @param xlim Limits for x-axis
##' @param ylim Limits for y-axis
##' @param axes Logical. If \code{TRUE} draw axes.
##' @param conf.int Logical. If \code{TRUE} draw confidence shadows.
##' @param ... Not yet used
##' @examples
##' library(survival)
##' d=sampleData(100,outcome="survival")
##' nd=sampleData(100,outcome="survival")
##' f1=coxph(Surv(time,event)~X1+X6+X8,data=d,x=TRUE,y=TRUE)
##' f2=coxph(Surv(time,event)~X2+X5+X9,data=d,x=TRUE,y=TRUE)
##' xx=Score(list(f1,f2), formula=Surv(time,event)~1,
##' data=nd, metrics="auc", null.model=FALSE, times=seq(3:10))
##' g <- autoplot(xx)
##' print(g)
##' aucgraph <- plotAUC(xx)
##' plotAUC(xx,conf.int=TRUE)
##' plotAUC(xx,which="contrasts")
##' plotAUC(xx,which="contrasts",conf.int=TRUE)
##' 
#'
#' @export
autoplot.Score <- function(object,
                    models,
                    type="score",
                    lwd=2,
                    xlim,
                    ylim,
                    axes=TRUE,
                    conf.int=FALSE,
                    ...){
    times=contrast=model=AUC=lower=upper=lower=upper=delta.AUC=reference=NULL
    ## cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    pframe <- switch(type,"score"={object$AUC$score},"contrasts"={object$AUC$contrasts},{stop("Type has to be either 'score' for AUC or 'contrasts' for differences in AUC.")})
    if (length(pframe$times)<2) stop(paste("Need at least two time points for plotting time-dependent AUC. Object has only ",length(pframe$times),"times"))
    if (type=="score"){
        ## AUC
        if (!missing(models)) pframe <- pframe[model %in% models]
        pframe[,lwd:=lwd]
        if (missing(xlim)) xlim <- pframe[,range(times)]
        if (missing(ylim)) ylim <- c(0.5,1)
        yticks <- seq(0,1,0.05)
        yticks <- yticks[yticks>=ylim[1] & yticks<=ylim[2]]
        pp <- ggplot2::ggplot(data=pframe,aes(times,AUC,fill=model,colour=model))
        pp + geom_line(size=lwd)
    }else{
        ## delta AUC
        pframe[,contrast:=paste(model,reference,sep=" - ")]
        pframe[,lwd:=lwd]
        if (missing(xlim)) xlim <- pframe[,range(times)]
        if (missing(ylim)) ylim <- c(min(pframe$lower),max(pframe$upper))
        yticks <- seq(-1,1,0.05)
        yticks <- yticks[yticks>=ylim[1] & yticks<=ylim[2]]
        pp <- ggplot(data=pframe,
                     aes(times,delta.AUC,fill=contrast,colour=contrast)) 
    }
    ## x-axis
    ## pp <- pp+ geom_segment(aes(x=xlim[1],xend=xlim[2],y=ylim[1],yend=ylim[1]))
    pp <- pp+theme_bw() %+replace% theme(axis.line = element_line(colour = "black"), 
                                         panel.grid.major = element_line(), panel.grid.major.x = element_blank(), 
                                         panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
                                         panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
                                         strip.background = element_rect(colour = "black", 
                                                                         size = 0.5), legend.key = element_blank())
    pp <- pp+ scale_fill_manual(values=cbbPalette)+ scale_colour_manual(values=cbbPalette)
    ## add the lines
    pp <- pp + geom_line(size=lwd) + xlim(xlim) + theme(legend.key = element_blank())
    ## y-axis
    pp <- pp + scale_y_continuous(expand=c(0,0),
                                  limits=ylim,
                                  breaks=yticks,
                                  labels=paste(round(100*yticks,1),"%"))
    if (conf.int==TRUE){
        if (type=="score"){
            ## pframe[,polygon(x=c(times,rev(times)),y=c(lower,rev(upper)),col=dimcol,border=NA),by=model]
            pp <- pp + geom_ribbon(aes(ymin=lower,ymax=upper,fill=model,linetype=NA),alpha=0.2)
        }else{
            pp <- pp + geom_ribbon(aes(ymin=lower,ymax=upper,fill=contrast,linetype=NA),alpha=0.2)
        }
    }
    pp
}

#----------------------------------------------------------------------
### plotAUC.R ends here
