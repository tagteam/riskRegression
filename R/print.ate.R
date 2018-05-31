### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: maj 31 2018 (18:09) 
##           By: Brice Ozenne
##     Update #: 150
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * print.ate (documentation)
#' @title Print Average Treatment Effects
#' @description Print average treatment effects.
#' @name print.ate
#' 
#' @param x object obtained with function \code{ate}
#' @param digits [integer, >0] Number of digits.
#' @param type [character vector] what to displayed.
#' Can be any combination of \code{"meanRisk"}, \code{"diffRisk"}, and \code{"ratioRisk"}.
#' @param ... passed to print
#'
#' @details to display confidence intervals/bands and p.value,
#' the \code{confint} method needs to be applied on the object.
#'
#' @seealso
#' \code{\link{confint.ate}} to compute confidence intervals/bands.
#' \code{\link{ate}} to compute the average treatment effects.

## * print.ate (code)
#' @rdname print.ate
#' @method print ate
#' @export
print.ate <- function(x, digits = 3, type = c("meanRisk","diffRisk","ratioRisk"), ...){

    ## check arguments
    type <- match.arg(type,c("meanRisk","diffRisk","ratioRisk"), several.ok = TRUE)
    
                                        # {{{ display
    if(!is.null(x$treatment)){
        cat("The treatment variable ",x$treatment," has the following options:\n",sep="")
        cat(paste(x$contrasts,collapse=", "),"\n")
    }


    if("meanRisk" %in% type){
        if(!is.null(x$treatment)){
            cat("\nMean risks on probability scale [0,1] in hypothetical worlds\nin which all subjects are treated with one of the treatment options:\n\n")
        }else{
            cat("Average risks within strata on probability scale [0,1]:\n\n")
        }
        print(x$meanRisk,digits=digits,...)
    }

    if(!is.null(x$treatment) && ("diffRisk" %in% type || "ratioRisk" %in% type)){
        id.name <- names(x$riskComparison)[1:3]

        cat("\nComparison of risks on probability scale [0,1] between\nhypothetical worlds are interpretated as if the treatment was randomized:\n\n")

        if("diffRisk" %in% type){
            cat("     > risk difference \n\n")
            ## simplify names
            name.diff <- grep("^diff",names(x$riskComparison),value = TRUE)
            dt.tempo <- x$riskComparison[,.SD,.SDcols = c(id.name,name.diff)]
            names(dt.tempo) <- gsub("diff.","",names(dt.tempo),fixed = TRUE)
            ## print
            print(dt.tempo,digits=digits,...)
        }

        if("ratioRisk" %in% type){
            cat("\n     > risk ratio \n\n")
            ## simplify names
            name.ratio <- grep("^ratio",names(x$riskComparison),value = TRUE)
            dt.tempo <- x$riskComparison[,.SD,.SDcols = c(id.name,name.ratio)]
            names(dt.tempo) <- gsub("ratio.","",names(dt.tempo),fixed = TRUE)
            ## print
            print(dt.tempo,digits=digits,...)
        }
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
            if("meanRisk" %in% type){
                cat("\nmean risk      : ",x$meanRisk.transform)
            }
            if("diffRisk" %in% type){
                cat("\nrisk difference: ",x$diffRisk.transform)
            }
            if("ratioRisk" %in% type){
                cat("\nrisk ratio     : ",x$ratioRisk.transform)
            }
        }
        cat("\nConfidence level:",x$conf.level,"\n")
    }
                                        # }}}
    
    ## export
    return(invisible(x))
}


#----------------------------------------------------------------------
### print.ate.R ends here
