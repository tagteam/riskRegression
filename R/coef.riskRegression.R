#' @S3method coef riskRegression
coef.riskRegression <- function(object,digits=3,eps=10^-4,...){

  cvars <- object$design$const
  Flevels <- object$factorLevels
  cat("\nTime constant regression coefficients:\n\n")
  if (is.null(object$timeConstantEffects)){
    cat("\nNone.\n")
    coefMat <- NULL
  }
  else{
    const.coef <- object$timeConstantEffects$coef
    const.se <- sqrt(diag(object$timeConstantEffects$var))
    wald <- const.coef/const.se
    waldp <- (1 - pnorm(abs(wald))) * 2
    format.waldp <- format.pval(waldp,digits=digits,eps=eps)
    names(format.waldp) <- names(waldp)
    format.waldp[const.se==0] <- NA
    if (any(const.se==0))
      warning("Some standard errors are zero. It seems that the model did not converge")
    coefMat <- do.call("rbind",lapply(cvars,function(v){
      covname <- strsplit(v,":")[[1]][[1]]
      if (is.null(Flevels[[covname]])){
        out <- c(v,signif(c(const.coef[v],exp(const.coef[v]),const.se[v],wald[v]),digits),format.waldp[v])
      }
      else{
        rlev <- object$refLevels[[covname]]
        out <- do.call("rbind",lapply(Flevels[[covname]],function(l){
          V <- paste(covname,l,sep=":")
          if (match(V,paste(covname,rlev,sep=":"),nomatch=FALSE))
            c(paste(covname,rlev,sep=":"),"--","--","--","--","--")
          else
            c(V,signif(c(const.coef[V],exp(const.coef[V]),const.se[V],wald[V]),digits),format.waldp[V])
        }))
      }
      out
    }))
    colnames(coefMat) <- c("Factor","Coef","exp(Coef)","StandardError","z","Pvalue")
    rownames(coefMat) <- rep("",NROW(coefMat))
    print(coefMat,quote=FALSE,right=TRUE)
    tp <- object$design$timepower[object$design$timepower>0]
    if (any(tp))
      cat(paste("\n\nNote:The coeffient(s) for the following variable(s)\n",
                sapply(names(tp),function(x){
                  paste("\t",x," (power=",as.character(tp[x]),")\n",sep="")}),
                "are interpreted as per factor unit  multiplied by time^power.\n",sep=""))
  }
  invisible(coefMat)
}
