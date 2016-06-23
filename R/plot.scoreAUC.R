### plot.scoreAUC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun 23 2016 (09:19) 
## Version: 
## last-updated: Jun 23 2016 (13:27) 
##           By: Thomas Alexander Gerds
##     Update #: 14
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
##' @param x Object obtained with \code{Score.list}
##' @param models Choice of models to plot
##' @param type Character. Either \code{"score"} to show AUC or \code{"test"} to show differences between AUC.
##' @param lwd Line width
##' @param xlim Limits for x-axis
##' @param ylim Limits for y-axis
##' @param col Color
##' @param axes Logical. If \code{TRUE} draw axes.
##' @param confint Logical. If \code{TRUE} draw confidence shadows.
##' @param ... Not yet used
#' @export
plot.scoreAUC <- function(x,models,type="score",lwd=2,xlim,ylim,col,axes=TRUE,confint=FALSE,...){
    times=model=AUC=lower.AUC=upper.AUC=NULL
    pframe <- switch(type,"score"={x$AUC$score},"test"={x$AUC$test},{stop("Type has to be either 'score' for AUC or 'test' for differences in AUC.")})
    ## plot(0,0,type="n",ylim = ylim,xlim = xlim,axes=FALSE,xlab = "Time",ylab = "AUC")
    ## if (axes){
    ## axis(1)
    ## prodlim::PercentAxis(2,at=seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/5))
    ## }
    ## pframe[,lines(times,AUC,type="l",lwd=lwd,col=col),by=model]
    ## pframe[,dimcol:=prodlim::dimColor(col[[1]],density=55),by=model]
    if (type=="score"){
        if (!missing(models)) pframe <- pframe[model %in% models]
        pframe[,col:=as.numeric(as.factor(model))]
        pframe[,lwd:=lwd]
        if (missing(xlim)) xlim <- pframe[,range(times)]
        if (missing(ylim)) ylim <- c(0.5,1)
        if (missing(col)) col <- 1:length(unique(pframe$model))
        pp <- ggplot(data=pframe,aes(times,AUC,fill=model,colour=model))
    }else{
        ## test
        pframe[,contrast:=paste(model,reference,sep=" - ")]
        pframe[,col:=as.numeric(as.factor(contrast))]
        pframe[,lwd:=lwd]
        if (missing(xlim)) xlim <- pframe[,range(times)]
        if (missing(ylim)) ylim <- c(min(pframe$lower),max(pframe$upper))
        if (missing(col)) col <- 1:length(unique(pframe$contrast))
        pp <- ggplot(data=pframe,aes(times,delta,fill=contrast,colour=contrast)) 
    }
    ## pp <- pp + theme_bw() + theme(axis.line = element_line(colour = "black"),
    ## panel.grid.major = element_blank(),
    ## panel.grid.minor = element_blank(),
    ## panel.border = element_blank(),
    ## panel.background = element_blank())
    ## geom_segment(aes(x=xlim[1],xend=xlim[1],y=ylim[1],yend=ylim[2]))
    ## x-axis
    ## pp <- pp+ geom_segment(aes(x=xlim[1],xend=xlim[2],y=ylim[1],yend=ylim[1]))
    pp <- pp+ theme_light()
    pp <- pp + geom_line(size=lwd) + xlim(xlim) + theme(legend.key = element_blank())
    ## + theme_classic()
    ## pp <- theme(axis.text.x= element_text(family, face, colour, size))
    ## pp <- pp+theme(axis.line= element_line(colour="black", size=1,linetype="solid"))
    ## pp <- pp+ theme(axis.line.x = element_line(color="black", size = 1),
    ## axis.line.y = element_line(color="black", size = 1,limits=c(.6,.8)))
    ## y-axis
    pp <- pp + scale_y_continuous(expand=c(0,0),limits=ylim,breaks=seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4),
                                  labels=paste(round(100*seq(ylim[1],ylim[2],(ylim[2]-ylim[1])/4),1),"%"))
    pp
    if (confint==TRUE){
        if (type=="score"){
            ## pframe[,polygon(x=c(times,rev(times)),y=c(lower.AUC,rev(upper.AUC)),col=dimcol,border=NA),by=model]
            pp <- pp + geom_ribbon(aes(ymin=lower.AUC,ymax=upper.AUC,fill=model,linetype=NA),alpha=0.2)
        }else{
            pp <- pp + geom_ribbon(aes(ymin=lower,ymax=upper,fill=contrast,linetype=NA),alpha=0.2)
        }
    }
    pp
}

#----------------------------------------------------------------------
### plot.scoreAUC.R ends here
