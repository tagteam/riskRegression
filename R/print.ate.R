### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: Mar 26 2018 (08:42) 
##           By: Thomas Alexander Gerds
##     Update #: 22
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
#' @param digits Number of digits
#' @param ... passed to print
#'
#' @method print ate
#' @export
print.ate <- function(x,digits=3,...){
    cat("The treatment variable ",x$treatment," has the following options:\n",sep="")
    cat(paste(x$contrasts,collapse=", "),"\n")
    cat("\nMean risks on probability scale [0,1] in hypothetical worlds\nin which all subjects are treated with one of the treatment options:\n\n")
    print(x$meanRisk,digits=digits,...)
    cat("\nComparison of risks on probability scale [0,1] between\nhypothetical worlds are interpretated as if the treatment was randomized:\n\n")    
    print(x$riskComparison,digits=digits,...)
    ##
    if(x$se && (x$conf.level > 0 && (x$conf.level < 1))){
        if(x$n.bootstrap==0){
            cat("\nWald confidence intervals are based on asymptotic standard errors.",sep="")
        }else {
            cat("\nBootstrap confidence intervals are based on ",x$n.bootstrap," bootstrap samples\nthat were drawn with replacement from the original data.",sep="")
        }
        if(x$nsim.band>0){
            cat("\nConfidence bands are based on ",x$nsim.band," simulations",sep="")
        }
        cat("\nConfidence level:",x$conf.level,"\n")
    }
    invisible(x)
}


#----------------------------------------------------------------------
### print.ate.R ends here
