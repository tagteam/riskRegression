plot.CauseSpecificCox <- function(x,
                                  cause,
                                  newdata,
                                  add=FALSE,
                                  xlab,
                                  ylab,
                                  xlim,
                                  ylim,
                                  lwd,
                                  col,
                                  lty,
                                  grid=FALSE,
                                  axes=TRUE,
                                  percent=TRUE,
                                  ...){
  # {{{ getting predicted risk
  if (missing(cause)) cause <- x$theCause
  plot.times <- x$response[,"time"]
  if (missing(newdata)){
    P1 <- predictEventProb(x,newdata=eval(x$call$data),time=plot.times)
    medianP1 <- P1[,sindex(plot.times,median(plot.times))]
    P1 <- P1[order(medianP1),]
    p1 <- P1[round(quantile(1:NROW(P1))),]
    warning("Argument newdata is missing.\n",
            "Shown are the cumulative incidence curves from the original data set.\nSelected are curves based on individual risk (min,q25,median,q75,max) at the median time:",
            median(plot.times))
  }
  else{
    p1 <- predict(x,newdata=newdata,time=plot.times)$cuminc
  }
  Y <- lapply(1:NROW(p1),function(i){p1[i,]})
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
  grid.DefaultArgs <- list(h=seq(0,1,0.05),col="gray77",lwd=2)
  axis1.DefaultArgs <- list()
  axis2.DefaultArgs <- list(at=seq(0,1,.25),side=2)
  legend.DefaultArgs <- list(legend=names(Y),lwd=lwd,col=col,lty=lty,cex=1.5,bty="n",y.intersp=1.3,x="topright",trimnames=TRUE)
  # }}}
  # {{{ smart control
  smartA <- prodlim:::SmartControl(call=  list(...),
                                   keys=c("plot","lines","legend","confint","grid","marktime","axis1","axis2"),
                                   ignore=c("x","type","cause","newdata","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","marktime","confint","automar","atrisk","timeOrigin","percent","axes","atrisk.args","confint.args","legend.args"),
                                   ignore.case=TRUE,
                                   defaults=list("plot"=plot.DefaultArgs,
                                     "axis1"=axis1.DefaultArgs,
                                     "axis2"=axis2.DefaultArgs,
                                     "grid"=grid.DefaultArgs,
                                     "legend"=legend.DefaultArgs,
                                     "lines"=lines.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),
                                     "axis1"=list(side=1)),
                                   verbose=TRUE)
  # }}}
  # {{{ empty plot
  if (!add) {
    do.call("plot",smartA$plot)
    if (grid==TRUE){
      do.call("abline",smartA$grid)
    }
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
  ## if(legend==TRUE && !add && !is.null(names(Y))){
  ## if (smartA$legend$trimnames==TRUE){
  ## smartA$legend$legend <- sapply(strsplit(names(Y),"="),function(x)x[[2]])
  ## smartA$legend$title <- unique(sapply(strsplit(names(Y),"="),function(x)x[[1]]))
  ## }
  ## smartA$legend <- smartA$legend[-match("trimnames",names(smartA$legend))]
  ## save.xpd <- par()$xpd
  ## par(xpd=TRUE)
  ## do.call("legend",smartA$legend)
  ## par(xpd=save.xpd)
  # }}}
}
