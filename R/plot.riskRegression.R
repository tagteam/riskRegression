#' Plotting predicted risk
#' 
#' Show predicted risk obtained by a risk prediction model as a function of
#' time.
#' 
#' 
#' @aliases plot.riskRegression plot.predictedRisk plot.CauseSpecificCox
#' @param x Fitted object obtained with one of \code{ARR}, \code{LRR},
#' \code{riskRegression}.
#' @param cause For CauseSpecificCox models the cause of interest.
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted risk.
#' @param xlim See \code{plot}
#' @param ylim See \code{plot}
#' @param xlab See \code{plot}
#' @param ylab See \code{plot}
#' @param lwd A vector of line thicknesses for the regression coefficients.
#' @param col A vector of colors for the regression coefficients.
#' @param lty A vector of line types for the regression coefficients.
#' @param axes Logical. If \code{FALSE} then do not draw axes.
#' @param percent If true the y-axis is labeled in percent.
#' @param legend If true draw a legend.
#' @param add Logical. If \code{TRUE} then add lines to an existing plot.
#' @param \dots Used for transclusion of smart arguments for \code{plot},
#' \code{lines}, \code{axis} and \code{background}. See function
#' \code{\link{SmartControl}} from prodlim.
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @keywords survival
#' @examples
#' 
#' 
#' library(pec)
#' data(Melanoma)
#' 
#' fit.arr <- ARR(Hist(time,status)~invasion+age+strata(sex),data=Melanoma,cause=1)
#' plot(fit.arr)
#' 
#' fit.csc <- CSC(Hist(time,status)~invasion+age+sex,data=Melanoma,cause=1)
#' plot(fit.csc)
#' 
#'
#' @S3method plot CauseSpecificCox
plot.riskRegression <- function(x,
                                cause,
                                newdata,
                                xlab,
                                ylab,
                                xlim,
                                ylim,
                                lwd,
                                col,
                                lty,
                                axes=TRUE,
                                percent=TRUE,
                                legend=TRUE,
                                add=FALSE,
                                ...){
  # {{{ getting predicted risk
  if (class(x)=="CauseSpecificCox")
    plot.times <- x$eventTimes
  else
    plot.times <- x$time
  if (class(x)=="predictedRisk")
    Y <- split(x$risk,1:NROW(x$risk))
  else{
    if (missing(newdata)){
      ff <- eval(x$call$formula)
      xdat <- unique(eval(x$call$data)[all.vars(rhs(ff))])
      if (NROW(xdat)<5){
        if (class(x)=="CauseSpecificCox"){
          p1 <- predictEventProb(x,newdata=xdat,times=plot.times,cause=cause)
        }
        else{
          p1 <- predict(x,newdata=xdat,times=plot.times)$risk}
        rownames(p1) <- paste("id",1:NROW(xdat))
      }
      else{
        if (class(x)=="CauseSpecificCox"){
          P1 <- predictEventProb(x,
                                 newdata=eval(x$call$data),
                                 times=plot.times,
                                 cause=cause)
        }
        else{
          P1 <- predict(x,newdata=eval(x$call$data),times=plot.times)$risk
        }
        medianP1 <- P1[,sindex(plot.times,median(plot.times))]
        P1 <- P1[order(medianP1),]
        p1 <- P1[round(quantile(1:NROW(P1))),]
        rownames(p1) <- paste("Predicted risk",c("Min","q25","Median","q75","Max"),sep="=")
        warning("Argument newdata is missing.\n",
                "Shown are the cumulative incidence curves from the original data set.\nSelected are curves based on individual risk (min,q25,median,q75,max) at the median time:",
                median(plot.times))
      }
    }
    else{
      p1 <- predict(x,newdata=newdata,time=plot.times)$risk
    }
    Y <- lapply(1:NROW(p1),function(i){p1[i,]})
    if (!is.null(rownames(p1)))
      names(Y) <- rownames(p1)
  }
  nlines <- NROW(Y)
  # }}}
  # {{{ labels, limits etc
  if (missing(ylab)) ylab <- "Cumulative incidence"
  if (missing(xlab)) xlab <- "Time"
  if (missing(xlim)) xlim <- c(0, max(plot.times))
  if (missing(ylim)) ylim <- c(0, 1)
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- 1:nlines
  if (missing(lty)) lty <- rep(1, nlines)
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  # }}}
  # {{{ defaults
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
  lines.DefaultArgs <- list(type="s")
  axis1.DefaultArgs <- list()
  axis2.DefaultArgs <- list(at=seq(0,1,.25),side=2)
  legend.DefaultArgs <- list(legend=names(Y),
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             x="topleft",
                             trimnames=TRUE)
  # }}}
  # {{{ smart control

  smartA <- SmartControl(call=  list(...),
                                   keys=c("plot","lines","legend","confint","marktime","axis1","axis2"),
                                   ignore=c("x","type","cause","newdata","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","marktime","confint","automar","atrisk","timeOrigin","percent","axes","atrisk.args","confint.args","legend.args"),
                                   ignore.case=TRUE,
                                   defaults=list("plot"=plot.DefaultArgs,
                                     "axis1"=axis1.DefaultArgs,
                                     "axis2"=axis2.DefaultArgs,
                                     "legend"=legend.DefaultArgs,
                                     "lines"=lines.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),
                                     "axis1"=list(side=1)),
                                   verbose=TRUE)

  # }}}
  # {{{ empty plot

  if (!add) {
    do.call("plot",smartA$plot)
    if (axes){
      do.call("axis",smartA$axis1)
      if (percent & is.null(smartA$axis1$labels))
        smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
      do.call("axis",smartA$axis2)
    }
  }

  # }}}
  # {{{ adding the lines 
  lines.type <- smartA$lines$type
  nix <- lapply(1:nlines, function(s) {
    lines(x = plot.times,
          y = Y[[s]],
          type = lines.type,
          col = col[s],
          lty = lty[s],
          lwd = lwd[s])
  })
  # }}}
  # {{{ legend
  if(legend==TRUE && !add && !is.null(names(Y))){
    if (smartA$legend$trimnames==TRUE && sapply((nlist <- strsplit(names(Y),"=")),function(x)length(x))==2){
      smartA$legend$legend <- sapply(nlist,function(x)x[[2]])
      smartA$legend$title <- unique(sapply(nlist,function(x)x[[1]]))
    }
    smartA$legend <- smartA$legend[-match("trimnames",names(smartA$legend))]
    save.xpd <- par()$xpd
    par(xpd=TRUE)
    do.call("legend",smartA$legend)
    par(xpd=save.xpd)
  }
  # }}}
}
