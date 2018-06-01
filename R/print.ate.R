### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: jun  1 2018 (17:16) 
##           By: Brice Ozenne
##     Update #: 199
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
        dt.tempo <- data.table::copy(x$meanRisk)
        order.col <- c(names(dt.tempo)[1:2],"meanRisk")
        if(!is.null(x$boot) && !is.null(x$conf.level)){
            order.col <- c(order.col,"meanRiskBoot")
        }

        ## round and merge column containing CI and CB
        numeric.col <- names(dt.tempo)[-(1:2)]
        dt.tempo[, c(numeric.col) := round(.SD, digits = digits) , .SDcols = numeric.col]        
        if(x$se && !is.null(x$conf.level)){
            dt.tempo[, c("conf.interval") := paste0("[",lower," ; ",upper,"]")]
            dt.tempo[,c("lower","upper") := NULL]
            order.col <- c(order.col,"se","conf.interval")
        }else if(x$se && is.null(x$boot)){
            order.col <- c(order.col,"se")
        }
        if(x$band && !is.null(x$conf.level)){
            dt.tempo[, c("conf.band") := paste0("[",lowerBand," ; ",upperBand,"]")]
            dt.tempo[,c("lowerBand","upperBand") := NULL]
            order.col <- c(order.col,"quantileBand","conf.band")
        }
        
        ## print
        data.table::setcolorder(dt.tempo, neworder = order.col)
        print(dt.tempo,...)
    }

    if(!is.null(x$treatment) && ("diffRisk" %in% type || "ratioRisk" %in% type)){
        id.name <- names(x$riskComparison)[1:3]

        cat("\nComparison of risks on probability scale [0,1] between\nhypothetical worlds are interpretated as if the treatment was randomized:\n\n")

        if("diffRisk" %in% type){
            cat("     > risk difference \n\n")
            ## only pick diff
            name.diff <- grep("^diff",names(x$riskComparison),value = TRUE)
            dt.tempo <- x$riskComparison[,.SD,.SDcols = c(id.name,name.diff)]

            ## simplify names
            names(dt.tempo) <- gsub("diff.","",names(dt.tempo),fixed = TRUE)
            order.col <- c(names(dt.tempo)[1:3],"diff")
            if(!is.null(x$boot) && !is.null(x$conf.level)){
                order.col <- c(order.col,"diffMeanBoot")
            }
             
            ## round and merge column containing CI and CB
            numeric.col <- setdiff(names(dt.tempo)[-(1:2)],"p.value")
            dt.tempo[, c(numeric.col) := round(.SD, digits = digits) , .SDcols = numeric.col]
            if(x$se && !is.null(x$conf.level)){
                dt.tempo[, c("conf.interval") := paste0("[",lower," ; ",upper,"]")]
                dt.tempo[,c("lower","upper") := NULL]
                order.col <- c(order.col,"se","conf.interval","p.value")
            }else if(x$se && is.null(x$boot)){
                order.col <- c(order.col,"se")
            }
            if(x$band && !is.null(x$conf.level)){
                dt.tempo[, c("conf.band") := paste0("[",lowerBand," ; ",upperBand,"]")]
                dt.tempo[,c("lowerBand","upperBand") := NULL]
                order.col <- c(order.col,"quantileBand","conf.band")
            }

            ## round p.value
            if(x$se && !is.null(x$conf.level)){
                dt.tempo$p.value <- sapply(dt.tempo$p.value, function(iP){
                    if(is.na(iP)){
                        as.numeric(NA)
                    }else if(iP<10^(-digits)){
                        return(paste0("<",10^(-digits)))
                    }else{
                        round(iP, digits = digits)
                    }
                })
            }
            
            ## print
            data.table::setcolorder(dt.tempo, neworder = order.col)
            print(dt.tempo,...)
        }

        if("ratioRisk" %in% type){
            cat("\n     > risk ratio \n\n")
            ## only pick ratio
            name.ratio <- grep("^ratio",names(x$riskComparison),value = TRUE)
            dt.tempo <- x$riskComparison[,.SD,.SDcols = c(id.name,name.ratio)]

            ## simplify names
            names(dt.tempo) <- gsub("ratio.","",names(dt.tempo),fixed = TRUE)
            order.col <- c(names(dt.tempo)[1:3],"ratio")
            if(!is.null(x$boot) && !is.null(x$conf.level)){
                order.col <- c(order.col,"ratioMeanBoot")
            }
            
            ## round and merge column containing CI and CB
            numeric.col <- setdiff(names(dt.tempo)[-(1:2)],"p.value")
            dt.tempo[, c(numeric.col) := round(.SD, digits = digits) , .SDcols = numeric.col]
            if(x$se && !is.null(x$conf.level)){
                dt.tempo[, c("conf.interval") := paste0("[",lower," ; ",upper,"]")]
                dt.tempo[,c("lower","upper") := NULL]
                order.col <- c(order.col,"se","conf.interval","p.value")
            }else if(x$se && is.null(x$boot)){
                order.col <- c(order.col,"se")
            }
            if(x$band && !is.null(x$conf.level)){
                dt.tempo[, c("conf.band") := paste0("[",lowerBand," ; ",upperBand,"]")]
                dt.tempo[,c("lowerBand","upperBand") := NULL]
                order.col <- c(order.col,"quantileBand","conf.band")
            }

            ## round p.value
            if(x$se && !is.null(x$conf.level)){
                dt.tempo$p.value <- sapply(dt.tempo$p.value, function(iP){
                    if(iP<10^(-digits)){return(paste0("<",10^(-digits)))}else{round(iP, digits = digits)}
                })
            }
            
            ## print
            data.table::setcolorder(dt.tempo, neworder = order.col)
            print(dt.tempo,...)
        }
    }
    ##
    if(x$se && !is.null(x$conf.level)){
        if(!is.null(x$boot)){
            bootci.method <- switch(x$bootci.method,
                                    "norm" = "Normal",
                                    "basic" = "Basic",
                                    "stud" = "Studentized",
                                    "perc" = "Percentile",
                                    "bca" = "BCa",
                                    "wald" = "Wald",
                                    "quantile" = "Percentile")
            cat("\n",bootci.method," bootstrap confidence intervals based on ",x$B," bootstrap samples\nthat were drawn with replacement from the original data.",sep="")
        }else {
            txt.band <- if(x$band){"/bands"}else{""}
            cat("\n",x$conf.level*100,"% Wald confidence intervals",txt.band," and p-values are based on asymptotic standard errors.",sep="")
            if(x$band){
                cat("\nQuantile for the confidence bands has been computed using ",x$nsim.band," simulations.",sep="")
            }
            cat("\nTransformation used to compute the confidence intervals/bands/p-values:",sep="")
            if("meanRisk" %in% type){
                cat("\n  ",x$meanRisk.transform," (mean risk)",sep="")
            }
            if("diffRisk" %in% type){
                cat("   ",x$diffRisk.transform, " (risk difference)",sep="")
            }
            if("ratioRisk" %in% type){
                cat("   ",x$ratioRisk.transform," (risk ratio)",sep="")
            }
        }
    }
    cat("\n\n")
                                        # }}}
    
    ## export
    return(invisible(x))
}


#----------------------------------------------------------------------
### print.ate.R ends here
