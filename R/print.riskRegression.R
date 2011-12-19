print.riskRegression <- function(object,
                                 times,
                                 digits=3,
                                 eps=10^-4,
                                 verbose=TRUE,
                                 ...) {
  # {{{ echo model type, IPCW and link function
  cat("\nriskRegression: Competing risks regression model \n")
  cat("\nIPCW estimation. The weights are based on\n",
      switch(object$censModel,
             "KM"={"the Kaplan-Meier estimate" },
             "cox"={"a Cox model" },                             
             "aalen"={"a non-parametric additive Aalen model"}),
      "for the censoring distribution.\n",sep="")
  ##   cat("\n",rep("_",options()$width/2),"\n",sep="")
  summary(object$response)
  ##   cat("\n",rep("_",options()$width/2),"\n",sep="")
  cat("\nLink function: \'",
      switch(object$link,"prop"="proportional","logistic"="logistic","additive"="additive","relative"="relative"),
      "\' yielding ",
      switch(object$link,
             "prop"="sub-hazard ratios (Fine & Gray 1999)",
             "logistic"="odds ratios",
             "additive"="absolute risk differences",
             "relative"="absolute risk ratios"),
      ", see help(riskRegression).\n",
      sep="")

  # }}}
  # {{{ find covariates and factor levels 

  cvars <- all.vars(object$design$const$formula)
  tvars <- all.vars(object$design$timevar$formula)
  Flevels <- object$factorLevels

  # }}}
  # {{{ time varying coefs

  if (!is.null(tvars)){
    cat("\nCovariates with time-varying effects:\n\n")
    nix <- lapply(tvars,function(v){
      if (is.null(flevs <- Flevels[[v]])){
        cat(" ",v," (numeric)\n",sep="")
      }
      else{
        cat(" ",v," (factor with levels: ",paste(flevs,collapse=", "),")\n",sep="")
      }
    })
  }
  if (length(tvars)==0){
    cat("No covariates with time-varying coefficient specified.\n")
  }
  else{
    cat("\nThe effects of these variables depend on time.")
  }
  cat("The column 'Intercept' is the baseline risk")
  cat(" where all the covariates have value zero\n\n")
  if (missing(times)) times <- quantile(object$time)
  showTimes <- sindex(eval.times=times,jump.times=object$time)
  showMat <- signif(exp(object$timeVaryingEffects$coef[showTimes,-1,drop=FALSE]),digits)
  rownames(showMat) <- signif(object$timeVaryingEffects$coef[showTimes,1],2)
  print(showMat)
  cat("\nShown are selected time points, use\n\nplot.riskRegression\n\nto investigate the full shape.\n\n")

  # }}}
  # {{{ time constant coefs
  if (!is.null(cvars)){
    cat("\nCovariates with time-constant effects:\n\n")
    nix <- lapply(cvars,function(v){
      if (is.null(flevs <- Flevels[[v]])){
        cat(" ",v," (numeric)\n",sep="")
      }
      else{
        cat(" ",v," (factor with levels: ",paste(flevs,collapse=", "),")\n",sep="")
      }
    })
  }
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
    if (!is.null(coefMat)){
      colnames(coefMat) <- c("Factor","Coef","exp(Coef)","StandardError","z","Pvalue")
      rownames(coefMat) <- rep("",NROW(coefMat))
      print(coefMat,quote=FALSE,right=TRUE)
      cat(paste("\n\nNote: The values exp(Coef) are",switch(object$link,
                                                            "prop"="sub-hazard ratios (Fine & Gray 1999)",
                                                            "logistic"="odds ratios",
                                                            "additive"="absolute risk differences",
                                                            "relative"="absolute risk ratios")),"\n")
      tp <- sapply(object$design$const$specialArguments,function(x)!is.null(x$power))
      if (any(tp))
        cat(paste("\n\nNote:The coeffient(s) for the following variable(s)\n",
                  sapply(names(object$design$const$specialArguments[tp]),function(x){
                    paste("\t",x," (power=",as.character(object$design$const$specialArguments[[x]]),")\n",sep="")}),
                  "are interpreted as per factor unit  multiplied by time^power.\n",sep=""))
    }
  }

  # }}}
  invisible(coefMat)
}
