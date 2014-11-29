#' @S3method summary riskRegression
summary.riskRegression <- function(object,
                                   times,
                                   digits=3,
                                   pvalue.digits=4,
                                   eps=10^-4,
                                   verbose=TRUE,
                                   ...) {
    # {{{ echo model type, IPCW and link function
    if (verbose){
        cat("\nriskRegression: Competing risks regression model \n")
        cat("\nIPCW estimation. The weights are based on\n",
            switch(object$censModel,
                   "KM"={"the Kaplan-Meier estimate" },
                   "cox"={"a Cox model" },                             
                   "aalen"={"a non-parametric additive Aalen model"}),
            " for the censoring distribution.\n",sep="")
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
    }
    # }}}
    # {{{ find covariates and factor levels
    cvars <- object$design$const
    tvars <- object$design$timevar
    Flevels <- object$factorLevels
    # }}}
    # {{{ time varying coefs
    if (verbose){
        if (!is.null(tvars)){
            cat("\nCovariates with time-varying effects:\n\n")
            nix <- lapply(tvars,function(v){
                if (is.null(flevs <- Flevels[[v]])){
                    cat(" ",v," (numeric)\n",sep="")
                }
                else{
                    cat(" ",v," (factor with levels: ",paste(flevs,collapse=", "),")\n",sep="")
                }
            })}
        if (length(tvars)==0){
            cat("No covariates with time-varying coefficient specified.\n")
        }
        else{
            cat("\nThe effects of these variables depend on time.")
        }
        cat("The column 'Intercept' is the baseline risk")
        cat(" where all the covariates have value zero\n\n")
    }
    if (missing(times)) times <- quantile(object$time)
    showTimes <- prodlim::sindex(eval.times=times,jump.times=object$time)
    showMat <- format(exp(object$timeVaryingEffects$coef[showTimes,-1,drop=FALSE]),digits=digits,nsmall=digits)
    rownames(showMat) <- signif(object$timeVaryingEffects$coef[showTimes,1],2)
    if (verbose){
        print(showMat)
        cat("\nShown are selected time points, use\n\nplot.riskRegression\n\nto investigate the full shape.\n\n")
    }
    # }}}
    # {{{ time constant coefs
    if (!is.null(cvars)){
        if (verbose){
            cat("\nCovariates with time-constant effects:\n\n")
            nix <- lapply(cvars,function(v){
                if (is.null(flevs <- Flevels[[v]])){
                    cat(" ",v," (numeric)\n",sep="")
                }
                else{
                    cat(" ",v," (factor with levels: ",paste(flevs,collapse=", "),")\n",sep="")
                }
            })}
    }
    if (verbose)
        cat("\nTime constant regression coefficients:\n\n")
    if (is.null(object$timeConstantEffects)){
        if (verbose)
            cat("\nNone.\n")
        coefMat <- NULL
    }
    else{
        const.coef <- object$timeConstantEffects$coef
        const.se <- sqrt(diag(object$timeConstantEffects$var))
        wald <- const.coef/const.se
        waldp <- (1 - pnorm(abs(wald))) * 2
        lower <- format(exp(const.coef-qnorm(0.975)*const.se),digits=digits,nsmall=digits)
        upper <- format(exp(const.coef+qnorm(0.975)*const.se),digits=digits,nsmall=digits)
        format.waldp <- format.pval(waldp,digits=pvalue.digits,eps=eps)
        names(format.waldp) <- names(waldp)
        format.waldp[const.se==0] <- NA
        if (any(const.se==0))
            warning("Some standard errors are zero. It seems that the model did not converge")
        coefMat <- do.call("rbind",lapply(cvars,function(v){
            covname <- strsplit(v,":")[[1]][[1]]
            if (is.null(Flevels[[covname]])){
                if (const.se[v]==0){
                    lower[v]=NA
                    upper[v]=NA
                }
                out <- c(v,
                         format(c(const.coef[v],exp(const.coef[v]),const.se[v],wald[v]),digits=digits,nsmall=digits),
                         paste("[",lower[v],";",upper[v],"]",sep=""),
                         format.waldp[v])
            } else{
                rlev <- object$refLevels[[covname]]
                out <- do.call("rbind",lapply(Flevels[[covname]],function(l){
                    V <- paste(covname,l,sep=":")
                    if (match(V,paste(covname,rlev,sep=":"),nomatch=FALSE))
                        c(paste(covname,rlev,sep=":"),"0","1","--","--","--","--")
                    else{
                        if (const.se[V]==0){
                            lower[V]=NA
                            upper[V]=NA
                        }
                        c(V,
                          format(c(const.coef[V],exp(const.coef[V]),const.se[V],wald[V]),digits=digits,nsmall=digits),
                          paste("[",lower[V],";",upper[V],"]",sep=""),
                          format.waldp[V])
                    }
                }))
            }
            out
        }))
        colnames(coefMat) <- c("Factor","Coef","exp(Coef)","StandardError","z","CI_95","Pvalue")
        rownames(coefMat) <- rep("",NROW(coefMat))
        if (verbose){
            print(coefMat,quote=FALSE,right=TRUE)
            cat(paste("\n\nNote: The values exp(Coef) are",switch(object$link,
                                                                  "prop"="sub-hazard ratios (Fine & Gray 1999)",
                                                                  "logistic"="odds ratios",
                                                                  "additive"="absolute risk differences",
                                                                  "relative"="absolute risk ratios")),"\n")
            tp <- x$design$timepower!=0
            if (any(tp))
                cat(paste("\n\nNote:The coeffient(s) for the following variable(s)\n",
                          paste(names(x$design$timepower)[tp],collapse=", "),
                          "are interpreted as per factor unit multiplied by time^power.\n",sep=""))
        }
    }
    # }}}
    invisible(coefMat)
}
