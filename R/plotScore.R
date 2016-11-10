### plotScore.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  7 2016 (06:39) 
## Version: 
## last-updated: Jun  5 2016 (10:44) 
##           By: Thomas Alexander Gerds
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' plot score object
#' 
#' plot score object
#' @param x Risk prediction score
#' @param models Choice of models to plot
#' @param xlab label for x-axis
#' @param ylab  label for y-axis
#' @param xlim limits for x-axis
#' @param ylim limits for x-axis
#' @param col line color 
#' @param lwd line width
#' @param lty line type
#' @param type Either \code{"AUC"} or \code{"Brier"}
#' @param add Logical. If \code{TRUE} add to existing plot
#' @param ... passed to plot
#'
#' @method plot highscore
#' @export
plot.highscore <- function(x,models,xlab,ylab,xlim,ylim,col,lwd=3,lty,type,add=FALSE,...){
    model=times=NULL
    X <- x$metric
    if (!is.null(x$score$times)){
        if (missing(models)) models <- unique(x$score$model)
        ttt <- x$score[model==models[[1]],]$times
        if (missing(xlab)) xlab <- "Time"
        if (missing(ylab)) ylab <- switch(X,"AUC"="AUC","Brier"="Brier score",X)
        if (missing(xlim)) xlim <- range(ttt)
        if (missing(ylim)) ylim <- switch(X,"AUC"={c(.5,1)},"Brier"={c(0,.25)},range(x$score[[X]]))
        if (missing(col)) col <- 1:length(models)
        if (missing(lty)) lty <- rep(1,length(models))
        if (missing(type))  type <- switch(X,"AUC"="l","Brier"="s","b")
        ## empty plot
        if (add==FALSE){
            plot(0,0,type="n",axes=FALSE,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
            axis(1)
            axis(2)
        }
        for (m in 1:length(models))
            lines(ttt,x$score[model==models[[m]]][[X]],col=col[[m]],lwd=lwd,type=type)
    }
}

#----------------------------------------------------------------------
### plotScore.R ends here
