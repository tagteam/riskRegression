### summary.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 29 2019 (13:18) 
## Version: 
## Last-Updated: jun 24 2020 (09:07) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.ate (documentation)
#' @title Summary Average Treatment Effects
#' @description Summary average treatment effects.
#' @name summary.ate
#' 
#' @param object object obtained with function \code{ate}
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
#' @rdname summary.ate
#' @method summary ate
#' @export
summary.ate <- function(object, digits = 3, type = c("meanRisk","diffRisk","ratioRisk"), estimator = object$estimator[1],
                        ...){
    short <- list(...)$short ## called by the print function
    lower <- upper <- lowerBand <- upperBand <- NULL ## [:CRANcheck:] data.table
    
    ## ** check arguments
    type <- match.arg(type,c("meanRisk","diffRisk","ratioRisk"), several.ok = TRUE)
    estimator <- match.arg(estimator, choices =  object$estimator, several.ok = FALSE)

    ## ** Display: specification of ate
    cat("    Estimation of the Average Treatment Effect \n\n")

    if(!identical(short,TRUE)){
        cat("- Event of interest               : ",object$event,"=",object$cause,"\n", sep = "")
        if(!is.null(object$treatment)){
            cat("- Levels of the treatment variable: ",object$treatment,"=", paste(object$contrasts,collapse=", "),"\n", sep = "")
        }
        cat("- Type of estimator               : ",switch(estimator,
                                                          "Gformula" = "G-formula",
                                                          "IPTW"= "inverse probability of treatment weighting",
                                                          "AIPTW" = "double robust"),"\n", sep = "")
        if(object$se[[1]] && !is.null(object$conf.level)){
            cat("- Statistical inference           : ")
            if(!is.null(object$boot)){
                bootci.method <- switch(object$bootci.method,
                                        "norm" = "Normal",
                                        "basic" = "Basic",
                                        "stud" = "Studentized",
                                        "perc" = "Percentile",
                                        "bca" = "BCa",
                                        "wald" = "Wald",
                                        "quantile" = "Percentile")
                cat(bootci.method," bootstrap based on ",object$B," bootstrap samples\n",
                    "                                    that were drawn with replacement from the original data.\n",sep="")
            }else {
                cat("using iid decomposition of the statistic (asymptotic normality, robust standard error) \n")
                if(object$meanRisk.transform!="none" && "meanRisk" %in% type){
                    cat("                                  : on the ",object$meanRisk.transform," scale for the mean risk \n",sep="")
                }
                if(object$diffRisk.transform!="none" && "diffRisk" %in% type){
                    cat("                                  : on the ",object$diffRisk.transform," scale for the risk difference \n",sep ="")
                }
                if(object$ratioRisk.transform!="none" && "ratioRisk" %in% type){
                    cat("                                  : on the ",object$ratioRisk.transform," scale for the risk ratio \n",sep="")
                }
            }
            cat("- Confidence level                : ", object$conf.level,"\n", sep ="")
            if(object$band){
                cat("- Confidence bands                : computed using ", object$nsim.band," simulations\n", sep ="")
            }
    
        }
    }
    
    ## ** Display: meanRisk
    if("meanRisk" %in% type){
        if(!identical(short,TRUE)){cat("\n")}
        cat("Average risk: between time zero and 'time',\n")

        if(!identical(short,TRUE)){
            if(!is.null(object$treatment)){
                cat("              in hypothetical worlds in which all subjects are treated with one of the treatment options,\n")
            }else{ 
                cat("             within strata,\n")
            }
            cat("              reported on the scale [0,1] (probability scale)\n\n")
        }
        
        keep.cols <- c(names(object$meanRisk)[!grepl("meanRisk",names(object$meanRisk))],
                       grep(estimator,names(object$meanRisk), value = TRUE))
        if(all(is.na(object$meanRisk$time))){
            keep.cols <- setdiff(keep.cols, "time")
        }
        dt.tempo <- data.table::copy(object$meanRisk[,.SD,.SDcols = keep.cols])
        ## order.col <- c(names(dt.tempo)[1:2],"meanRisk")
        ## if(!is.null(object$boot) && !is.null(object$conf.level)){
        ## order.col <- c(order.col,"bootstrap")
        ## }

        ## simplify names (needs to be done in two steps)
        names(dt.tempo) <- gsub("meanRisk\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
        names(dt.tempo)[names(dt.tempo)=="meanRisk"] <- "average risk"

        ## merge into CI and CB
        if(!is.null(object$conf.level)){
            if(object$se){
                dt.tempo[, c("conf.interval") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lower),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upper),
                                                        "]")]
                dt.tempo[,c("lower","upper") := NULL]
                ## order.col <- c(order.col,"se","conf.interval")
            }
            if(object$band){
                dt.tempo[, c("conf.band") := paste0("[",
                                                    sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                    " ; ",
                                                    sprintf(paste0("%1.",digits,"f"),upperBand),
                                                    "]")]
                dt.tempo[,c("lowerBand","upperBand") := NULL]
                ## order.col <- c(order.col,"quantileBand","conf.band")
            }
        }else if(is.null(object$boot)){
            ## if(object$se){
            ## order.col <- c(order.col,"se")
            ## }
        }

        ## print
        ## data.table::setcolorder(dt.tempo, neworder = order.col)
        print(dt.tempo,digits=digits,...)
        cat("\n")
    }

    if(!is.null(object$treatment) && ("diffRisk" %in% type || "ratioRisk" %in% type)){
        if("diffRisk" %in% type){

            if(!identical(short,TRUE)){cat("\n")}
            cat("Difference in risks: (B-A) between time zero and 'time',\n")

            if(!identical(short,TRUE)){
                if(!is.null(object$treatment)){
                    cat("                      comparing an hypothetical world in which all subjects are treated with one treatment option (A),\n",
                        "                            to an hypothetical world in which all subjects are treated with the other treatment options (B),\n")
                }else{  
                    cat("                      between two strata,\n")
                }
                cat("                      reported on the scale [-1,1] (difference between two probabilities)\n\n")
            }
            
            ## only pick diff
            keep.cols <- c(names(object$riskComparison)[!grepl("diff|ratio",names(object$riskComparison))],
                           grep(paste0("diff\\.",estimator),names(object$riskComparison), value = TRUE))
            if(all(is.na(object$riskComparison$time))){
                keep.cols <- setdiff(keep.cols, "time")
            }
            dt.tempo <- object$riskComparison[,.SD,.SDcols = keep.cols]
            ## order.col <- c(names(dt.tempo)[1:3],"diff")
            ## if(!is.null(object$boot) && !is.null(object$conf.level)){
            ## order.col <- c(order.col,"bootstrap")
            ## }

            ## simplify names (needs to be done in two steps)
            names(dt.tempo) <- gsub("diff\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
            names(dt.tempo)[names(dt.tempo)=="diff"] <- "risk difference"

            ## merge into CI and CB
            if(!is.null(object$conf.level)){
                if(object$se){
                    dt.tempo[, c("conf.interval") := paste0("[",
                                                            sprintf(paste0("%1.",digits,"f"),lower),
                                                            " ; ",
                                                            sprintf(paste0("%1.",digits,"f"),upper),
                                                            "]")]
                    dt.tempo[,c("lower","upper") := NULL]
                    ## order.col <- c(order.col,"se","conf.interval","p.value")
                }
                if(object$band){
                    dt.tempo[, c("conf.band") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upperBand),
                                                        "]")]
                    dt.tempo[,c("lowerBand","upperBand") := NULL]
                    ## order.col <- c(order.col,"quantileBand","conf.band")
                }
            }else if(is.null(object$boot)){
                ## if(object$se){
                    ## order.col <- c(order.col,"se")
                ## }
            }

            if(object$se[[1]] && !is.null(object$conf.level)){
                dt.tempo$p.value <- format.pval(dt.tempo$p.value,digits=digits,eps=10^{-digits})
            }

            ## print
            ## data.table::setcolorder(dt.tempo, neworder = order.col)
            print(dt.tempo,digits=digits,...)
            cat("\n")
        }
        if("ratioRisk" %in% type){

            if(!identical(short,TRUE)){cat("\n")}
            cat("Ratio of risks: (B/A) between time zero and 'time',\n")

            if(!identical(short,TRUE)){
                if(!is.null(object$treatment)){
                    cat("                comparing an hypothetical world in which all subjects are treated with one treatment option (A),\n",
                        "                      to an hypothetical world in which all subjects are treated with the other treatment options (B),\n")
                }else{
                    cat("                between two strata,\n")
                }
                cat("                reported on the scale ]0,+oo[ (ratio of two probabilities):\n\n")
            }
            
            ## only pick ratio
            keep.cols <- c(names(object$riskComparison)[!grepl("diff|ratio",names(object$riskComparison))],
                           grep(paste0("ratio\\.",estimator),names(object$riskComparison), value = TRUE))
            if(all(is.na(object$riskComparison$time))){
                keep.cols <- setdiff(keep.cols, "time")
            }
            dt.tempo <- object$riskComparison[,.SD,.SDcols = keep.cols]
            ## order.col <- c(names(dt.tempo)[1:3],"ratio")
            ## if(!is.null(object$boot) && !is.null(object$conf.level)){
            ## order.col <- c(order.col,"bootstrap")
            ## }
            
            ## simplify names (needs to be done in two steps)
            names(dt.tempo) <- gsub("ratio\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
            names(dt.tempo)[names(dt.tempo)=="ratio"] <- "risk ratio"

            ## merge into CI and CB
            if(!is.null(object$conf.level)){
                if(object$se){
                    dt.tempo[, c("conf.interval") := paste0("[",
                                                            sprintf(paste0("%1.",digits,"f"),lower),
                                                            " ; ",
                                                            sprintf(paste0("%1.",digits,"f"),upper),
                                                            "]")]
                    dt.tempo[,c("lower","upper") := NULL]
                    ## order.col <- c(order.col,"se","conf.interval","p.value")
                }
                if(object$band){
                    dt.tempo[, c("conf.band") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upperBand),
                                                        "]")]
                    dt.tempo[,c("lowerBand","upperBand") := NULL]
                    ## order.col <- c(order.col,"quantileBand","conf.band")
                }
            }else if(is.null(object$boot)){
                ## if(object$se){
                    ## order.col <- c(order.col,"se")
                ## }
            }

            if(object$se[[1]] && !is.null(object$conf.level)){
                dt.tempo$p.value <- format.pval(dt.tempo$p.value,digits=digits,eps=10^{-digits})
            }
            
            ## print
            ## data.table::setcolorder(dt.tempo, neworder = order.col)
            print(dt.tempo,digits=digits,...)
        }
        cat("\n")
    }

    ## export
    return(invisible(object))
}


######################################################################
### summary.ate.R ends here
