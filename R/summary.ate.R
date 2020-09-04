### summary.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 29 2019 (13:18) 
## Version: 
## Last-Updated: sep  4 2020 (10:37) 
##           By: Brice Ozenne
##     Update #: 79
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
#' @param se [logical] should the standard error of the risks be displayed?
#' @param quantile [logical] should the quantile of the confidence bands be displayed?
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
summary.ate <- function(object,  estimator = object$estimator[1],
                        type = c("meanRisk","diffRisk","ratioRisk"), se = object$se, quantile = FALSE,
                        digits = 3,
                        ...){
    short <- list(...)$short ## called by the print function
    lower <- upper <- lowerBand <- upperBand <- NULL ## [:CRANcheck:] data.table
    allContrasts <- object$allContrasts
    contrasts <- attr(allContrasts,"contrasts")

    test.se <- !is.null(object$se) && object$se
    test.ci <- !is.null(object$ci) && object$ci
    test.band <- !is.null(object$band) && object$band
    test.p.value <- !is.null(object$p.value) && object$p.value
    
    ## ** check arguments
    type <- match.arg(type,c("meanRisk","diffRisk","ratioRisk"), several.ok = TRUE)
    estimator <- match.arg(estimator, choices =  object$estimator, several.ok = FALSE)

    ## ** keep columns
    keep.cols <- NULL
    if(test.se){
        keep.cols <- c(keep.cols,"se")
    }
    if(test.ci){
        keep.cols <- c(keep.cols,"lower","upper")
        if(test.p.value){keep.cols <- c(keep.cols,"p.value")}
    }
    if(test.band){
        if(object$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
            if(quantile){
                keep.cols <- c(keep.cols,"quantileBand")
            }        
            keep.cols <- c(keep.cols,"lowerBand","upperBand")
        }
        if(test.p.value){keep.cols <- c(keep.cols,"adj.p.value")}
    }
    
    ## ** Display: specification of ate
    if(length(object$causes)>1){        
        cat("    Estimation of the Average Treatment Effect for cause ",object$theCause," \n\n",sep="")
    }else{
        cat("    Estimation of the Average Treatment Effect \n\n",sep="")
    }
    
    if(!identical(short,TRUE)){
        cat("- Event of interest               : ",object$event,"=",object$cause,"\n", sep = "")
        if(!is.null(object$treatment)){
            cat("- Levels of the treatment variable: ",object$treatment,"=", paste(object$contrasts,collapse=", "),"\n", sep = "")
        }
        cat("- Type of estimator               : ",switch(estimator,
                                                          "Gformula" = "G-formula",
                                                          "IPTW"= "inverse probability of treatment weighting",
                                                          "AIPTW" = "double robust"),"\n", sep = "")
        if(test.ci || test.band){
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
                cat("using iid decomposition of the statistic (robust standard errors) \n")
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

            txt.alternative <- switch(object$alternative,
                                      "two.sided" = "two-sided tests, alternative: unequal risk",
                                      "greater" = "one-sided tests, alternative: greater risk",
                                      "less" = "one-sided tests, alternative: less risk")
            cat("- Confidence level                : ", object$conf.level," (",txt.alternative,")\n", sep ="")
            if(test.band){
                if(object$method.band=="maxT-simulation"){
                    cat("- Confidence bands/adj. p-values  : single-step maxT computed using ", object$n.sim," simulations\n", sep ="")
                }else if(object$method.band=="maxT-integration"){
                    cat("- Confidence bands/adj. p-values  : single-step maxT computed using numerical integration\n", sep ="")
                }else{
                    cat("- Confidence bands/adj. p-values  : ",object$method.band," method\n", sep ="")
                }
            }
    
        }
        if("meanRisk" %in% type || !is.null(object$treatment) && ("diffRisk" %in% type || "ratioRisk" %in% type)){cat("\n")}        
    }

    ## ** Display: meanRisk
    if("meanRisk" %in% type){
        cat("Average risk: between time zero and 'time',\n")

        if(!identical(short,TRUE)){
            if(!is.null(object$treatment)){
                cat("              in hypothetical worlds in which all subjects are treated with one of the treatment options,\n")
            }else{ 
                cat("             within strata,\n")
            }
            cat("              reported on the scale [0,1] (probability scale)\n\n")
        }
        keep.cols.mR <- c(names(object$meanRisk)[!grepl("meanRisk",names(object$meanRisk))],
                          paste0("meanRisk.",estimator),
                          if(length(keep.cols)>0){paste0("meanRisk.",estimator,".",setdiff(keep.cols,c("p.value","adj.p.value")))}
                          )
        if(all(is.na(object$meanRisk$time))){
            keep.cols.mR <- setdiff(keep.cols.mR, "time")
        }
        dt.tempo <- data.table::copy(object$meanRisk[,.SD,.SDcols = keep.cols.mR])
        if(!is.null(contrasts)){
            dt.tempo <- dt.tempo[dt.tempo[[1]] %in% contrasts]
        }
        ## order.col <- c(names(dt.tempo)[1:2],"meanRisk")
        ## if(!is.null(object$boot) && !is.null(object$conf.level)){
        ## order.col <- c(order.col,"bootstrap")
        ## }

        ## simplify names (needs to be done in two steps)
        names(dt.tempo) <- gsub("meanRisk\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
        names(dt.tempo)[names(dt.tempo)=="meanRisk"] <- "average risk"

        ## merge into CI and CB
        if(test.ci){
            dt.tempo[, c("lower") := paste0("[",
                                                    sprintf(paste0("%1.",digits,"f"),lower),
                                                    " ; ",
                                                    sprintf(paste0("%1.",digits,"f"),upper),
                                                    "]")]
            dt.tempo[,c("upper") := NULL]
            data.table::setnames(dt.tempo, old = "lower", new = "conf.interval")
            ## order.col <- c(order.col,"se","conf.interval")
        }
        if(test.band && object$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
            dt.tempo[, c("lowerBand") := paste0("[",
                                                sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                " ; ",
                                                sprintf(paste0("%1.",digits,"f"),upperBand),
                                                "]")]
            dt.tempo[,c("upperBand") := NULL]
            data.table::setnames(dt.tempo, old = "lowerBand", new = "conf.band")
            ## order.col <- c(order.col,"quantileBand","conf.band")
        }
        ## print
        ## data.table::setcolorder(dt.tempo, neworder = order.col)
        print(dt.tempo,digits=digits,...)
        if(!is.null(object$treatment) && ("diffRisk" %in% type || "ratioRisk" %in% type)){cat("\n")}
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
            keep.cols.dR <- c(names(object$riskComparison)[!grepl("diff|ratio",names(object$riskComparison))],
                              paste0("diff.",estimator),
                              if(length(keep.cols)>0){paste0("diff.",estimator,".",keep.cols)})
            if(all(is.na(object$riskComparison$time))){
                keep.cols.dR <- setdiff(keep.cols.dR, "time")
            }
            dt.tempo <- object$riskComparison[,.SD,.SDcols = keep.cols.dR]
            if(!is.null(allContrasts)){
                dt.tempo <- dt.tempo[paste0(dt.tempo[[1]],".",dt.tempo[[2]]) %in% paste0(allContrasts[1,],".",allContrasts[2,])]
            }
            ## order.col <- c(names(dt.tempo)[1:3],"diff")
            ## if(!is.null(object$boot) && !is.null(object$conf.level)){
            ## order.col <- c(order.col,"bootstrap")
            ## }

            ## simplify names (needs to be done in two steps)
            names(dt.tempo) <- gsub("diff\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
            names(dt.tempo)[names(dt.tempo)=="diff"] <- "risk difference"

            ## merge into CI and CB
            if(test.ci){
                dt.tempo[, c("lower") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lower),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upper),
                                                        "]")]
                dt.tempo[,c("upper") := NULL]
                data.table::setnames(dt.tempo, old = "lower", new = "conf.interval")
                                       
                ## order.col <- c(order.col,"se","conf.interval","p.value")
                if(test.p.value){
                    dt.tempo$p.value <- format.pval(dt.tempo$p.value,digits=digits,eps=10^{-digits})
                }
            }
            if(test.band){
                if(object$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
                    dt.tempo[, c("lowerBand") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upperBand),
                                                        "]")]
                    dt.tempo[,c("upperBand") := NULL]
                    data.table::setnames(dt.tempo, old = "lowerBand", new = "conf.band")
                }
                ## order.col <- c(order.col,"quantileBand","conf.band")
                if(test.p.value){
                    dt.tempo$adj.p.value <- format.pval(dt.tempo$adj.p.value,digits=digits,eps=10^{-digits})
                }
            }


            ## print
            print(dt.tempo,digits=digits,...)
            if("ratioRisk" %in% type){cat("\n")}
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
            keep.cols.rR <- c(names(object$riskComparison)[!grepl("diff|ratio",names(object$riskComparison))],
                              paste0("ratio.",estimator),
                              if(length(keep.cols)>0){paste0("ratio.",estimator,".",keep.cols)})
            if(all(is.na(object$riskComparison$time))){
                keep.cols.rR <- setdiff(keep.cols.rR, "time")
            }
            dt.tempo <- object$riskComparison[,.SD,.SDcols = keep.cols.rR]
            if(!is.null(allContrasts)){
                dt.tempo <- dt.tempo[paste0(dt.tempo[[1]],".",dt.tempo[[2]]) %in% paste0(allContrasts[1,],".",allContrasts[2,])]
            }
            ## order.col <- c(names(dt.tempo)[1:3],"ratio")
            ## if(!is.null(object$boot) && !is.null(object$conf.level)){
            ## order.col <- c(order.col,"bootstrap")
            ## }
            
            ## simplify names (needs to be done in two steps)
            names(dt.tempo) <- gsub("ratio\\.","",gsub(paste0("\\.",estimator),"",names(dt.tempo)))
            names(dt.tempo)[names(dt.tempo)=="ratio"] <- "risk ratio"

            ## merge into CI and CB
            if(test.ci){
                dt.tempo[, c("lower") := paste0("[",
                                                sprintf(paste0("%1.",digits,"f"),lower),
                                                " ; ",
                                                sprintf(paste0("%1.",digits,"f"),upper),
                                                "]")]
                dt.tempo[,c("upper") := NULL]
                data.table::setnames(dt.tempo, old = "lower", new = "conf.interval")
                ## order.col <- c(order.col,"se","conf.interval","p.value")
                if(test.p.value){
                    dt.tempo$p.value <- format.pval(dt.tempo$p.value,digits=digits,eps=10^{-digits})
                }
            }
            if(test.band){
                if(object$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
                    dt.tempo[, c("lowerBand") := paste0("[",
                                                        sprintf(paste0("%1.",digits,"f"),lowerBand),
                                                        " ; ",
                                                        sprintf(paste0("%1.",digits,"f"),upperBand),
                                                        "]")]
                    dt.tempo[,c("upperBand") := NULL]
                    data.table::setnames(dt.tempo, old = "lowerBand", new = "conf.band")
                }
                ## order.col <- c(order.col,"quantileBand","conf.band")
                if(test.p.value){
                    dt.tempo$adj.p.value <- format.pval(dt.tempo$adj.p.value,digits=digits,eps=10^{-digits})
                }
            }

            
            ## print
            ## data.table::setcolorder(dt.tempo, neworder = order.col)
            print(dt.tempo,digits=digits,...)
        }
    }

    ## export
    return(invisible(object))
}


######################################################################
### summary.ate.R ends here
