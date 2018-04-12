### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: apr 12 2018 (12:18) 
##           By: Brice Ozenne
##     Update #: 111
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Print average treatment effects
#'
#' Print average treatment effects
#' @param x object obtained with function \code{ate}
#' @param bootci.method Character. Method for constructing bootstrap confidence intervals.
#' Either "perc" (the default), "norm", "basic", "stud", or "bca".
#' Argument passed to \code{boot::boot.ci}.
#' @param digits Number of digits
#' @param ... passed to print
#'
#' @details When using bootstrap resampling the p-values are computing using a test-inversion method,
#' i.e. find the critical confidence level such that one side of the confidence interval overlap the null hypothesis.
#' The p-value is 1 minus the critical confidence level.
#' 
#' @method print ate
#' @export
print.ate <- function(x, bootci.method = x$bootci.method, digits = 3, ...){

                                        # {{{ compute confidence intervals and p-values when using the bootstrap
    if(!is.null(x$boot) && bootci.method != x$bootci.method){
        out <- calcCIboot(boot = x$boot,
                          meanRisk = x$meanRisk,
                          riskComparison = x$riskComparison,
                          type = bootci.method,
                          conf = x$conf.level,
                          TD = x$TD)
        x$bootci.method <- bootci.method
        x$meanRisk <- out$meanRisk
        x$riskComparison <- out$riskComparison

    }
                                        # }}}

                                        # {{{ display
    cat("The treatment variable ",x$treatment," has the following options:\n",sep="")
    cat(paste(x$contrasts,collapse=", "),"\n")
    cat("\nMean risks on probability scale [0,1] in hypothetical worlds\nin which all subjects are treated with one of the treatment options:\n\n")
    print(x$meanRisk,digits=digits,...)
    cat("\nComparison of risks on probability scale [0,1] between\nhypothetical worlds are interpretated as if the treatment was randomized:\n\n")    
    print(x$riskComparison,digits=digits,...)
    ##
    if(x$se && (x$conf.level > 0 && (x$conf.level < 1))){
        if(x$B==0){
            cat("\nWald confidence intervals are based on asymptotic standard errors.",sep="")
        }else {
            bootci.method <- switch(bootci.method,
                                    "norm" = "Normal",
                                    "basic" = "Basic",
                                    "stud" = "Studentized",
                                    "perc" = "Percentile",
                                    "bca" = "BCa")
            cat("\n",bootci.method," bootstrap confidence intervals based on ",x$B," bootstrap samples\nthat were drawn with replacement from the original data.",sep="")
        }
        if(x$nsim.band>0){
            cat("\nConfidence bands are based on ",x$nsim.band," simulations",sep="")
        }
        cat("\nConfidence level:",x$conf.level,"\n")
    }
                                        # }}}
    
    ## export
    return(invisible(x))
}


#----------------------------------------------------------------------
### print.ate.R ends here
