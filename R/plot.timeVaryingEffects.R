plot.timeVaryingEffects <- function(x,
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
                                    ...){
  if (missing(formula)){
    formula <- x$formula
  }
  time <- x$coef[,1,drop=TRUE]
  var <- all.vars(formula)[1]
  matchVar <- grep(var,colnames(x$coef))
  coef <- x$coef[,matchVar,drop=FALSE]
  se <- sqrt(x$var[,matchVar,drop=FALSE])
  zval <- qnorm(1- (1-confint)/2, 0,1)
  lower <- coef-zval*se
  upper <- coef + zval*se
  if (missing(ylim))
    ylim <- c(floor(min(lower)),ceiling(max(upper)))
  if (missing(xlim))
    xlim=c(0,max(time))
  if (!add){
    plot(0,0,type="n",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,...)
  }
  if (missing(col))
    col <- 1:12
  col <- rep(col,length.out=NCOL(coef))
  if (missing(lty))
    lty <- 1
  lty <- rep(lty,length.out=NCOL(coef))
  if (missing(lwd))
    lwd <- 2
  lwd <- rep(lwd,length.out=NCOL(coef))
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
      lines(time,coef[,i],col=col[l],lwd=lwd[l],lty=lty[l])
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
    lines(time,coef,lwd=lwd[l],lty=lty[l],col=col[1])
    ## confidence shadows
    ccrgb=as.list(col2rgb(col[1],alpha=TRUE))
    names(ccrgb) <- c("red","green","blue","alpha")
    ccrgb$alpha=55
    cc=do.call("rgb",c(ccrgb,list(max=255)))
    polygon(x=c(time,rev(time)),y=c(lower,rev(upper)),col=cc,border=NA)
  }
}
