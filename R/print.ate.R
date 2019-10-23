### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: okt 23 2019 (16:25) 
##           By: Brice Ozenne
##     Update #: 334
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
#' @param estimator [character] The type of estimator relative to which the estimates should be displayed. 
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
print.ate <- function(x, digits = 3, type = c("meanRisk","diffRisk","ratioRisk"), estimator = x$estimator[1], ...){

    lower <- upper <- lowerBand <- upperBand <- NULL ## [:CRANcheck:] data.table
    
    ## ** check arguments
    type <- match.arg(type,c("meanRisk","diffRisk","ratioRisk"), several.ok = TRUE)
    estimator <- match.arg(estimator, choices =  x$estimator, several.ok = FALSE)

    ## ** Display: specification of ate
    cat("    Estimation of the Average Treatment Effect \n\n")

    cat("- Event of interest               : ",x$event,"=",x$cause,"\n", sep = "")
    if(!is.null(x$treatment)){
        cat("- Levels of the treatment variable: ",x$treatment,"=", paste(x$contrasts,collapse=", "),"\n", sep = "")
    }
    cat("- Type of estimator               : ",switch(estimator,
                                                      "Gformula" = "G-formula",
                                                      "IPTW"= "inverse probability of treatment weighting",
                                                      "AIPTW" = "double robust"),"\n", sep = "")
    if(x$se[[1]] && !is.null(x$conf.level)){
        cat("- Statistical inference           : ")
        if(!is.null(x$boot)){
            bootci.method <- switch(x$bootci.method,
                                    "norm" = "Normal",
                                    "basic" = "Basic",
                                    "stud" = "Studentized",
                                    "perc" = "Percentile",
                                    "bca" = "BCa",
                                    "wald" = "Wald",
                                    "quantile" = "Percentile")
            cat(bootci.method," bootstrap based on ",x$B," bootstrap samples\n",
                "that were drawn with replacement from the original data.",sep="")
        }else {
            cat("using iid decomposition of the statistic (asymptotic normality, robust standard error) \n")
            if(x$meanRisk.transform!="none" && "meanRisk" %in% type){
                cat("                                  : on the ",x$meanRisk.transform," scale for the mean risk \n",sep="")
            }
            if(x$diffRisk.transform!="none" && "diffRisk" %in% type){
                cat("                                  : on the ",x$diffRisk.transform," scale for the risk difference \n",sep ="")
            }
            if(x$ratioRisk.transform!="none" && "ratioRisk" %in% type){
                cat("                                  : on the ",x$ratioRisk.transform," scale for the risk ratio \n",sep="")
            }
        }
        cat("- Confidence level                : ", x$conf.level,"\n", sep ="")
        if(x$band){
            cat("- Confidence bands                : computed using ", x$nsim.band," simulations", sep ="")
        }
    
    }
    cat("\n\n")

    ## ** Display: meanRisk
    if("meanRisk" %in% type){
        cat("Average risk: between time zero and 'time',\n")
        if(!is.null(x$treatment)){
            cat("              in hypothetical worlds in which all subjects are treated with one of the treatment options,\n")
        }else{ 
            cat("             within strata,\n")
        }
        cat("              reported on the scale [0,1] (probability scale)\n\n")
        keep.cols <- c(names(x$meanRisk)[!grepl("meanRisk",names(x$meanRisk))],
                       grep(estimator,names(x$meanRisk), value = TRUE))
                       
        dt.tempo <- data.table::copy(x$meanRisk[,.SD,.SDcols = keep.cols])
        ## order.col <- c(names(dt.tempo)[1:2],"meanRisk")
        ## if(!is.null(x$boot) && !is.null(x$conf.level)){
            ## order.col <- c(order.col,"bootstrap")
        ## }

        ## simplify names (needs to be done in two steps)
        names(dt.tempo) <- gsub("meanRisk\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
        names(dt.tempo)[names(dt.tempo)=="meanRisk"] <- "Average risk"

        ## merge into CI and CB
        if(!is.null(x$conf.level)){
            if(x$se){
                dt.tempo[, c("conf.interval") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lower),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upper),
                                                        "]")]
                dt.tempo[,c("lower","upper") := NULL]
                ## order.col <- c(order.col,"se","conf.interval")
            }
            if(x$band){
                dt.tempo[, c("conf.band") := paste0("[",
                                                    sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                    " ; ",
                                                    sprintf(paste0("%1.",digits,"f"),upperBand),
                                                    "]")]
                dt.tempo[,c("lowerBand","upperBand") := NULL]
                ## order.col <- c(order.col,"quantileBand","conf.band")
            }
        }else if(is.null(x$boot)){
            ## if(x$se){
                ## order.col <- c(order.col,"se")
            ## }
        }

        ## print
        ## data.table::setcolorder(dt.tempo, neworder = order.col)
        print(dt.tempo,digits=digits,...)
        cat("\n")
    }

    if(!is.null(x$treatment) && ("diffRisk" %in% type || "ratioRisk" %in% type)){
        if("diffRisk" %in% type){
            cat("\n")
            cat("Difference in risks: (B-A) between time zero and 'time',\n")
            if(!is.null(x$treatment)){
                cat("                      comparing an hypothetical world in which all subjects are treated with one treatment option (A),\n",
                    "                            to an hypothetical world in which all subjects are treated with the other treatment options (B),\n")
            }else{  
                cat("                      between two strata,\n")
            }
            cat("                      reported on the scale [-1,1] (difference of two probabilities)\n\n")

            ## only pick diff
            keep.cols <- c(names(x$riskComparison)[!grepl("diff|ratio",names(x$riskComparison))],
                           grep(paste0("diff\\.",estimator),names(x$riskComparison), value = TRUE))
            dt.tempo <- x$riskComparison[,.SD,.SDcols = keep.cols]
            ## order.col <- c(names(dt.tempo)[1:3],"diff")
            ## if(!is.null(x$boot) && !is.null(x$conf.level)){
            ## order.col <- c(order.col,"bootstrap")
            ## }

            ## simplify names (needs to be done in two steps)
            names(dt.tempo) <- gsub("diff\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
            names(dt.tempo)[names(dt.tempo)=="diff"] <- "risk difference"

            ## merge into CI and CB
            if(!is.null(x$conf.level)){
                if(x$se){
                    dt.tempo[, c("conf.interval") := paste0("[",
                                                            sprintf(paste0("%1.",digits,"f"),lower),
                                                            " ; ",
                                                            sprintf(paste0("%1.",digits,"f"),upper),
                                                            "]")]
                    dt.tempo[,c("lower","upper") := NULL]
                    ## order.col <- c(order.col,"se","conf.interval","p.value")
                }
                if(x$band){
                    dt.tempo[, c("conf.band") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upperBand),
                                                        "]")]
                    dt.tempo[,c("lowerBand","upperBand") := NULL]
                    ## order.col <- c(order.col,"quantileBand","conf.band")
                }
            }else if(is.null(x$boot)){
                ## if(x$se){
                    ## order.col <- c(order.col,"se")
                ## }
            }

            if(x$se[[1]] && !is.null(x$conf.level)){
                dt.tempo$p.value <- format.pval(dt.tempo$p.value,digits=digits,eps=10^{-digits})
            }

            ## print
            ## data.table::setcolorder(dt.tempo, neworder = order.col)
            print(dt.tempo,digits=digits,...)
            cat("\n")
        }
        if("ratioRisk" %in% type){
            cat("\n")
            cat("Ratio of risks: (B/A) between time zero and 'time',\n")
            if(!is.null(x$treatment)){
                cat("                comparing an hypothetical world in which all subjects are treated with one treatment option (A),\n",
                    "                       to an hypothetical world in which all subjects are treated with the other treatment options (B),\n")
            }else{
                cat("                between two strata,\n")
            }
            cat("                reported on the scale ]0,+oo[ (ratio of two probabilities):\n\n")

            ## only pick ratio
            keep.cols <- c(names(x$riskComparison)[!grepl("diff|ratio",names(x$riskComparison))],
                           grep(paste0("ratio\\.",estimator),names(x$riskComparison), value = TRUE))
            dt.tempo <- x$riskComparison[,.SD,.SDcols = keep.cols]
            ## order.col <- c(names(dt.tempo)[1:3],"ratio")
            ## if(!is.null(x$boot) && !is.null(x$conf.level)){
            ## order.col <- c(order.col,"bootstrap")
            ## }
            
            ## simplify names (needs to be done in two steps)
            names(dt.tempo) <- gsub("ratio\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
            names(dt.tempo)[names(dt.tempo)=="ratio"] <- "risk ratio"

            ## merge into CI and CB
            if(!is.null(x$conf.level)){
                if(x$se){
                    dt.tempo[, c("conf.interval") := paste0("[",
                                                            sprintf(paste0("%1.",digits,"f"),lower),
                                                            " ; ",
                                                            sprintf(paste0("%1.",digits,"f"),upper),
                                                            "]")]
                    dt.tempo[,c("lower","upper") := NULL]
                    ## order.col <- c(order.col,"se","conf.interval","p.value")
                }
                if(x$band){
                    dt.tempo[, c("conf.band") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upperBand),
                                                        "]")]
                    dt.tempo[,c("lowerBand","upperBand") := NULL]
                    ## order.col <- c(order.col,"quantileBand","conf.band")
                }
            }else if(is.null(x$boot)){
                ## if(x$se){
                    ## order.col <- c(order.col,"se")
                ## }
            }

            if(x$se[[1]] && !is.null(x$conf.level)){
                dt.tempo$p.value <- format.pval(dt.tempo$p.value,digits=digits,eps=10^{-digits})
            }
            
            ## print
            ## data.table::setcolorder(dt.tempo, neworder = order.col)
            print(dt.tempo,digits=digits,...)
        }
        cat("\n")
    }

    # }}}
    ## export
    return(invisible(x))
}


#----------------------------------------------------------------------
### print.ate.R ends here
