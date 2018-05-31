### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: maj 31 2018 (11:53) 
##           By: Brice Ozenne
##     Update #: 133
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
#' @details When using bootstrap resampling the p-values are computing using a test-inversion method,
#' i.e. find the critical confidence level such that one side of the confidence interval overlap the null hypothesis.
#' The p-value is 1 minus the critical confidence level.
#' 
#' @method print ate
#' @export
print.ate <- function(x, bootci.method = x$bootci.method, digits = 3, ...){

                                        # {{{ display
    if(!is.null(x$treatment)){
        cat("The treatment variable ",x$treatment," has the following options:\n",sep="")
        cat(paste(x$contrasts,collapse=", "),"\n")
    }


    if(!is.null(x$treatment)){
        cat("\nMean risks on probability scale [0,1] in hypothetical worlds\nin which all subjects are treated with one of the treatment options:\n\n")
    }else{
        cat("Average risks within strata on probability scale [0,1]:\n\n")
    }
    print(x$meanRisk,digits=digits,...)
    if(!is.null(x$treatment)){
        id.name <- names(x$riskComparison)[1:3]

        cat("\nComparison of risks on probability scale [0,1] between\nhypothetical worlds are interpretated as if the treatment was randomized:\n\n")

        cat("     > risk difference \n\n")
        ## simplify names
        name.diff <- grep("^diff",names(x$riskComparison),value = TRUE)
        dt.tempo <- x$riskComparison[,.SD,.SDcols = c(id.name,name.diff)]
        names(dt.tempo) <- gsub("diff.","",names(dt.tempo),fixed = TRUE)
        ## print
        print(dt.tempo,digits=digits,...)
        
        cat("\n     > risk ratio \n\n")
        ## simplify names
        name.ratio <- grep("^ratio",names(x$riskComparison),value = TRUE)
        dt.tempo <- x$riskComparison[,.SD,.SDcols = c(id.name,name.ratio)]
        names(dt.tempo) <- gsub("ratio.","",names(dt.tempo),fixed = TRUE)
        ## print
        print(dt.tempo,digits=digits,...)
    }
    ##
    if(x$se && !is.null(x$conf.level)){
        if(!is.null(x$boot)){            
            type.boot <- switch(x$type.boot,
                                "norm" = "Normal",
                                "basic" = "Basic",
                                "stud" = "Studentized",
                                "perc" = "Percentile",
                                "bca" = "BCa")
            cat("\n",type.boot," bootstrap confidence intervals based on ",x$B," bootstrap samples\nthat were drawn with replacement from the original data.",sep="")
        }else {
            cat("\nWald confidence intervals and p-values are based on asymptotic standard errors.",sep="")
            cat("\nConfidence bands are based on ",x$nsim.band," simulations",sep="")
            cat("\nTransformation used to compute the confidence intervals/bands/p-values:",sep="")
            cat("\nmean risk      : ",object$ate.transform)
            cat("\nrisk difference: ",object$diffAte.transform)
            cat("\nrisk ratio     : ",object$ratioAte.transform)
        }
        cat("\nConfidence level:",x$conf.level,"\n")
    }
                                        # }}}
    
    ## export
    return(invisible(x))
}


#----------------------------------------------------------------------
### print.ate.R ends here
