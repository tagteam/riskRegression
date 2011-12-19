plotCoef <- function(x,formula,confint=.95,ylim,xlab="Time",ylab="Cumulative coefficient",add=FALSE,...){
  if (missing(formula)){
    formula <- x$design$timevar$formula
  }
  var <- all.vars(formula)[1]
  matchVar <- grep(var,colnames(x$timeVarCoef))
  coef <- x$timeVarCoef[,matchVar,drop=FALSE]
  ## var <- apply(x$timeVarVar[,matchVar,drop=FALSE],1,function(x){
  ## xx <- x
  ## xx[xx<0] <- 0
  ## xx
  ## })
  ## se <- t(sqrt(var))
  se <- sqrt(x$timeVarVar[,matchVar,drop=FALSE])
  zval <- qnorm(1- (1-confint)/2, 0,1)
  lower <- coef-zval*se
  upper <- coef + zval*se
  browser()
  if (missing(ylim))
      ylim <- c(floor(min(lower)),ceiling(max(upper)))
  if (!add)
  plot(0,0,type="n",xlab=xlab,ylab=ylab,xlim=c(0,max(time)),ylim=ylim)
  col <- 1:12
  if (length(matchVar)>1){
    ref <- x$refLevels[var]
    abline(h=0,col=1,lwd=2)
    nix <- lapply(1:NCOL(coef),function(i){
      lines(time,coef[,i],col=col[1+i],lwd=2)
      ## confidence shadows
      ccrgb=as.list(col2rgb(col[1+i],alpha=TRUE))
      names(ccrgb) <- c("red","green","blue","alpha")
      ccrgb$alpha=55
      cc=do.call("rgb",c(ccrgb,list(max=255)))
      polygon(x=c(time,rev(time)),y=c(lower[,i],rev(upper[,i])),col=cc,border=NA)
    })
  }
  else{
    abline(h=0,col="gray55",lwd=2)
    lines(time,coef,col=col[1],lwd=2)
    ## confidence shadows
    ccrgb=as.list(col2rgb(col[1],alpha=TRUE))
    names(ccrgb) <- c("red","green","blue","alpha")
    ccrgb$alpha=55
    cc=do.call("rgb",c(ccrgb,list(max=255)))
    polygon(x=c(time,rev(time)),y=c(lower,rev(upper)),col=cc,border=NA)
  }
}
