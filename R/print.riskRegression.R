#' Print function for riskRegression models
#'
#' Print function for riskRegression models
#' @param x Object obtained with ARR, LRR or riskRegression
#' @param times Time points at which to show time-dependent
#'     coefficients
#' @param digits Number of digits for all numbers but p-values
#' @param eps p-values smaller than this number are shown as such
#' @param verbose Level of verbosity
#' @param conf.int level of confidence. default is 0.95
#' @param ... not used
#'
#' @method print riskRegression
#' @export
print.riskRegression <- function(x,
                                 times,
                                 digits=3,
                                 eps=10^-4,
                                 verbose=TRUE,
                                 conf.int=.95,
                                 ...) {
    alpha <- 1-conf.int
    # {{{ echo model type, IPCW and link function
    cat("Competing risks regression model \n")
    cat("\nIPCW weights: ",
        switch(tolower(x$cens.model),
               "km"={"marginal Kaplan-Meier" },
               "cox"={"Cox regression model" },                             
               "aalen"={"non-parametric additive Aalen model"}),
        " for the censoring distribution.",sep="")
    summary(x$response)
    cat("\nLink: \'",
        switch(x$link,"prop"="cloglog","logistic"="logistic","additive"="linear","relative"="log"),
        "\' yielding ",
        switch(x$link,
               "prop"="sub-hazard ratios (Fine & Gray 1999)",
               "logistic"="odds ratios",
               "additive"="absolute risk differences",
               "relative"="absolute risk ratios"),
        ## ", see help(riskRegression).\n",
        sep="")

    # }}}
    # {{{ find covariates and factor levels 
    cvars <- x$design$const
    if (Ipos <- match("Intercept",x$design$timevar,nomatch=0))
        tvars <- x$design$timevar[-Ipos]
    else
        tvars <- x$design$timevar
    Flevels <- x$factorLevels
    
    # }}}
    # {{{ time varying coefs
    if (length(tvars)>0){
        cat("\nCovariates with time-varying effects:\n\n")
        nix <- lapply(tvars,function(v){
            if (is.null(flevs <- Flevels[[v]])){
                cat(" ",v," (numeric)\n",sep="")
            }
            else{
                cat(" ",v," (factor with levels: ",paste(flevs,collapse=", "),")\n",sep="")
            }
        })
    } else{
        cat("\nNo covariates with time-varying coefficient specified.\n")
    }
    if (is.null(x$timeConstantEffects$coef)){
        cat("\nNo time constant regression coefficients in model.\n")
        coefMat <- NULL
    }
    else{
        cat("\nTime constant regression coefficients:\n")
        const.coef <- x$timeConstantEffects$coef
        const.se <- sqrt(diag(x$timeConstantEffects$var))
        wald <- const.coef/const.se
        waldp <- (1 - pnorm(abs(wald))) * 2
        format.waldp <- format.pval(waldp,digits=digits,eps=eps)
        names(format.waldp) <- names(waldp)
        format.waldp[const.se==0] <- NA
        if (any(const.se==0))
            warning("Some standard errors are zero. It seems that the model did not converge")
        lower <- const.coef + qnorm(alpha/2) * const.se
        upper <- const.coef + qnorm(1- alpha/2) * const.se
        levs <- rep("",length(cvars))
        ## levs <- sapply(cvars,function(v)Flevels[[v]])
        Flevels <- Flevels[attr(x$timeConstantEffects,"variables")]
        if (length(Flevels)>0){
            fvars <- lapply(1:length(Flevels),function(i){paste0(names(Flevels)[[i]],Flevels[[i]])})
            names(fvars) <- names(Flevels)
            for (f in names(fvars)){
                pos.active <- match(fvars[[f]][[2]],cvars,nomatch=0)
                if (pos.active==0) {
                    warning("Problems with factor variable names or levels.")
                }else{
                    cvars[pos.active:(pos.active+length(fvars[[f]])-2)] <- c(f,rep("",length(fvars[[f]])-2))
                    levs[pos.active:(pos.active+length(fvars[[f]])-2)] <- Flevels[[f]][-1]
                }
            }
        }
        coefMat <- data.frame(Variable=cvars,
                              Levels=levs,
                              Coef=exp(const.coef),
                              Lower=exp(lower),
                              Upper=exp(upper),
                              Pvalue=waldp)
        rownames(coefMat) <- NULL
        if (!is.null(coefMat)){
            pcoefMat <- coefMat
            pcoefMat$Pvalue <- format.waldp
            print(pcoefMat,row.names=FALSE,quote=FALSE,right=TRUE,digits=digits)
            cat(paste("\n\nNote: The coefficients (Coef) are",switch(x$link,"prop"="sub-hazard ratios (Fine & Gray 1999)","logistic"="odds ratios","additive"="absolute risk differences","relative"="absolute risk ratios")),"\n")
            tp <- x$design$timepower!=0
            if (any(tp))
                cat(paste("\n\nNote:The coeffient(s) for the variable(s)\n",
                          paste(names(x$design$timepower)[tp],collapse=", "),
                          " are to be interpreted as effect per unit multiplied by time^power.\n",sep=""))
        }
    }
    # }}}
    invisible(coefMat)
}
