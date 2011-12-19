plotEffects <- function(x,
                        formula,
                        level,
                        refLine=TRUE,
                        confint=.95,
                        xlim,
                        ylim,
                        xlab="Time",
                        ylab="Cumulative coefficient",
                        col,
                        lty,
                        lwd,
                        add=FALSE,
                        legend,
                        axes=TRUE,
                        ...){
  # {{{ find variables and coefficients with confidence limits

  timevars <- all.vars(x$design$timevar$formula)
  if (length(timevars)==0) stop("No variable with time-varying effect in model design.")
  if (missing(formula)) formula <- x$design$timevar$formula
  var <- all.vars(formula)[1]
  time <- x$time
  matchVar <- grep(var,colnames(x$timeVaryingEffects$coef))
  matchVarnames <- grep(var,colnames(x$timeVaryingEffects$coef),val=TRUE)
  coef <- x$timeVaryingEffects$coef[,matchVar,drop=FALSE]
  se <- sqrt(x$timeVaryingEffects$var[,matchVar,drop=FALSE])
  zval <- qnorm(1- (1-confint)/2, 0,1)
  lower <- coef-zval*se
  upper <- coef + zval*se

  # }}}
  # {{{  plotting limits, colors, etc
  if (missing(ylim))
    ylim <- c(floor(min(lower)),ceiling(max(upper)))
  if (missing(xlim))
    xlim=c(0,max(time))
  if (missing(col))
    col <- 1:12
  col <- rep(col,length.out=NCOL(coef))
  if (missing(lty))
    lty <- 1
  lty <- rep(lty,length.out=NCOL(coef))
  if (missing(lwd))
    lwd <- 2
  lwd <- rep(lwd,length.out=NCOL(coef))
  # }}}
  # {{{  setting default arguments
  if (missing(legend)) legend <- length(matchVar)>1
  background.DefaultArgs <- list(xlim=xlim,ylim=ylim,horizontal=seq(0,1,.25),vertical=NULL,bg="white",fg="gray88")
  axis1.DefaultArgs <- list()
  axis2.DefaultArgs <- list(side=2)
  lines.DefaultArgs <- list(type="s")
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
  
  legend.DefaultArgs <- list(legend=matchVarnames,
                             trimnames=FALSE,
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             x="topright")
  ## confint.DefaultArgs <- list(x=x,newdata=newdata,type=type,citype="shadow",times=plot.times,cause=cause,density=55,col=col[1:nlines],lwd=rep(2,nlines),lty=rep(3,nlines))
  # }}}
  # {{{ smart control

  smartA <- SmartControl(call=  list(...),
                         keys=c("plot","lines","legend","confint","background","axis1","axis2"),
                         ignore=c("x","formula","refLine","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","confint","percent","axes"),
                         defaults=list("plot"=plot.DefaultArgs,
                           "lines"=lines.DefaultArgs,
                           "legend"=legend.DefaultArgs,
                           ## "confint"=confint.DefaultArgs,
                           "background"=background.DefaultArgs,
                           "axis1"=axis1.DefaultArgs,
                           "axis2"=axis2.DefaultArgs),
                         forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                         ignore.case=TRUE,
                         replaceDefaults=FALSE,
                         verbose=TRUE)

  # }}}
  # {{{  plot and backGround
  if (!add) {
    do.call("plot",smartA$plot)
  }
  # }}}
  # {{{  axes

  if (!add) {
    if (axes){
      do.call("axis",smartA$axis1)
      do.call("axis",smartA$axis2)
    }
  }

  # }}}
  # {{{  legend

  if(legend==TRUE && !add && !is.null(matchVarnames)){

    if (smartA$legend$trimnames==TRUE){
      smartA$legend$legend <- sapply(strsplit(names(Y),"="),function(x)x[[2]])
    }
    smartA$legend <- smartA$legend[-match("trimnames",names(smartA$legend))]
    save.xpd <- par()$xpd
    par(xpd=TRUE)
    do.call("legend",smartA$legend)
    par(xpd=save.xpd)
  }

  # }}}
  # {{{ adding the lines

  if (length(matchVar)>1){
    ref <- x$refLevels[var]
    if (refLine==TRUE)
      abline(h=0,col="gray55",lwd=2)
    levs <- colnames(coef)
    if (!missing(level)) {
      levs <- levs[match(level,levs,nomatch=0)]
    }
    if (length(levs)==0) stop(paste("Could not find level(s): ",paste(level,collapse=", ")),"\nAvailable levels: ",paste(colnames(coef),collapse=", "))
    nix <- lapply(1:length(levs),function(l){
      i <- match(levs[l],colnames(coef),nomatch=0)
      lines(time,coef[,i],col=col[l],lwd=lwd[l],lty=lty[l],type="s")
      ## confidence shadows
      ccrgb=as.list(col2rgb(col[l],alpha=TRUE))
      names(ccrgb) <- c("red","green","blue","alpha")
      ccrgb$alpha=55
      cc=do.call("rgb",c(ccrgb,list(max=255)))
      polygon(x=c(time,rev(time)),y=c(lower[,i],rev(upper[,i])),col=cc,border=NA)
    })
  }
  else{
    if (refLine==TRUE)
      abline(h=0,col="gray55",lwd=2)
    lines(time,coef,lwd=lwd[1],lty=lty[1],col=col[1],type="s")
    ## confidence shadows
    ccrgb=as.list(col2rgb(col[1],alpha=TRUE))
    names(ccrgb) <- c("red","green","blue","alpha")
    ccrgb$alpha=55
    cc=do.call("rgb",c(ccrgb,list(max=255)))
    polygon(x=c(time,rev(time)),y=c(lower,rev(upper)),col=cc,border=NA)
  }

  # }}}
}
