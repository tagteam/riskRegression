### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: Apr 20 2017 (18:21) 
##           By: Thomas Alexander Gerds
##     Update #: 13
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
    if( (x$conf.level > 0 && (x$conf.level < 1)) ){

        if(x$n.bootstrap==0){
            cat("\n",100*x$conf.level,"% Wald confidence intervals.\n",sep="")
        }else {
            cat("\n",100*x$conf.level,"% bootstrap confidence intervals are based on ",x$n.bootstrap," bootstrap samples\nthat were drawn with replacement from the original data.\n",sep="")
        }
    }
    invisible(x)
    }


#----------------------------------------------------------------------
### print.ate.R ends here
