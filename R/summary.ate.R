### summary.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 29 2019 (13:18) 
## Version: 
## Last-Updated: okt 26 2020 (09:39) 
##           By: Brice Ozenne
##     Update #: 267
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
#' @param type [character vector] what to displayed.
#' Can be \code{"meanRisk"} to display the risks specific to each treatment group,
#' \code{"diffRisk"} to display the difference in risks between treatment groups,
#' or \code{"ratioRisk"} to display the ratio of risks between treatment groups,.
#' @param estimator [character] The type of estimator relative to which the estimates should be displayed. 
#' @param se [logical] should the standard error of the risks be displayed?
#' @param quantile [logical] should the quantile of the confidence bands be displayed?
#' @param estimate.boot [logical] should the average estimate on the bootstrap samples be displayed?
#' @param digits [integer, >0] Number of digits.
#' @param short [logical] If \code{TRUE}, only displays the estimated risks.
#' @param ... passed to confint
#'
#' @details to display confidence intervals/bands and p.value,
#' the \code{confint} method needs to be applied on the object.
#'
#' @seealso
#' \code{\link{as.data.table}} to extract the estimates in a \code{data.table} object.
#' \code{\link{autoplot.ate}} for a graphical representation the standardized risks.
#' \code{\link{confint.ate}} to compute p-values and adjusted p-values
#' or perform statistical inference using a transformation.
#' \code{\link{confint.ate}} to compute (pointwise/simultaneous) confidence intervals and (unadjusted/adjusted) p-values, possibly using a transformation.

## * print.ate (code)
#' @rdname summary.ate
#' @method summary ate
#' @export
summary.ate <- function(object,  estimator = object$estimator[1], short = FALSE,
                        type = c("meanRisk","diffRisk"), 
                        se = FALSE, quantile = FALSE, estimate.boot = TRUE,
                        digits = 3,
                        ...){

    lower <- upper <- lowerBand <- upperBand <- NULL ## [:CRANcheck:] data.table

    ## ** update confint
    args.confint <- c("conf.level","n.sim","contrasts","allContrasts","meanRisk.transform","diffRisk.transform","ratioRisk.transform","seed",
                      "ci","band","p.value","method.band","alternative","bootci.method")
    if(any(args.confint %in% names(list(...)))){
        object[c("meanRisk","diffRisk","ratioRisk","inference","inference.allContrasts","inference.contrasts","transform")] <- confint(object,...)
    }
        
    if(!is.null(object$inference.allContrasts)){
        allContrasts <- object$inference.allContrasts
        contrasts <- object$inference.contrasts
    }else{
        contrasts <- object$contrasts
        allContrasts <- utils::combn(contrasts, m = 2)
    }

    ## ** check arguments
    type <- match.arg(type,c("meanRisk","diffRisk","ratioRisk"), several.ok = TRUE)
    if("diffRisk" %in% type && "ratioRisk" %in% type){
        stop("Cannot simulatenously display the risk difference and the risk ratio \n")
    }
    estimator <- match.arg(estimator, choices =  object$estimator, several.ok = FALSE)

    ## ** Display: specification of ate
    if(short==0){
        print(object, display.results = FALSE)

        if(object$inference$ci || object$inference$band || object$inference$p.value){
            cat("\n Testing procedure\n")
            index.transform <- which(object$transform!="none")
            if(length(index.transform)>0){
                name.tempo <- c("mean risk","risk difference","risk ratio")
                cat("   (on the ",paste(object$transform[index.transform],collapse="/")," scale for the ",paste(name.tempo[index.transform],collapse="/"),") \n",sep="")
            }

            if(object$inference$band == 1){
                if(object$inference$alternative=="two.sided"){
                    cat(" - Null hypothesis     : given two treatments (A,B), equal risks at all timepoints \n", sep ="")            
                }else if(object$inference$alternative=="greater"){
                    cat(" - Null hypothesis     : given two treatments (A,B), \n",
                        "                         risk under B is equal or smaller than the risk under A at all timepoints \n", sep ="")            
                }else if(object$inference$alternative=="less"){
                    cat(" - Null hypothesis     : given two treatments (A,B), \n",
                        "                         risk under B is equal or greater than the risk under A at all timepoints \n", sep ="")            
                }
            }else if(object$inference$band == 2){
                if(object$inference$alternative=="two.sided"){
                    cat(" - Null hypothesis     : equal risks at all timepoints for all treatments \n", sep ="")            
                }else if(object$inference$alternative=="greater"){
                    cat(" - Null hypothesis     : risks under the experimental treatments are all equal or smaller than the risks under the reference treatment(s) \n", sep ="")            
                }else if(object$inference$alternative=="less"){
                    cat(" - Null hypothesis     : risks under the experimental treatments are all equal or greater than the risks under the reference treatment(s) \n", sep ="")            
                }
            }else if(object$inference$ci){
                if(object$inference$alternative=="two.sided"){
                    cat(" - Null hypothesis     : given two treatments (A,B) and a specific timepoint, equal risks \n", sep ="")            
                }else if(object$inference$alternative=="greater"){
                    cat(" - Null hypothesis     : given two treatments (A,B) and a specific timepoint, \n",
                        "                         risk under B is equal or smaller than the risk under A \n", sep ="")            
                }else if(object$inference$alternative=="less"){
                    cat(" - Null hypothesis     : given two treatments (A,B) and a specific timepoint, \n",
                        "                         risk under B is equal or greater than the risk under A \n", sep ="")            
                }
            }
            cat(" - Confidence level    : ", object$inference$conf.level,"\n", sep ="")
            if(object$inference$band){
                if(is.na(object$inference$method.band)){
                    cat(" - Multiple comparisons: no adjustment\n", sep ="")
                }else if(object$inference$method.band=="maxT-simulation"){
                    cat(" - Multiple comparisons: single-step max-T adjustment computed using ", object$inference$n.sim," simulations\n", sep ="")
                }else if(object$inference$method.band=="maxT-integration"){
                    cat(" - Multiple comparisons: single-step max-T adjustment computed using numerical intergration\n", sep ="")
                }else{
                    cat(" - Multiple comparisons: ",object$inference$method.band," adjustment\n", sep ="")
                }
            }
        }
    }

    ##     ## ** prepare
    ## txt.alternative <- switch(object$inference$alternative,
    ##                           "two.sided" = "two-sided tests, alternative: unequal risk",
    ##                           "greater" = "one-sided tests, alternative: greater risk",
    ##                           "less" = "one-sided tests, alternative: less risk")

    if(short==0){
        cat("\n Results: \n")
    }
    if(identical(type,"meanRisk")){
        ## ** only risks
        ## prepare
        iIndexRow <- which((object$meanRisk$estimator == estimator) * (object$meanRisk$treatment %in% contrasts) == 1)

        dt.tempo <- object$meanRisk[iIndexRow,.SD, .SDcols = setdiff(names(object$meanRisk),"estimator")]
        if(!is.na(object$variable["strata"])){
            setnames(dt.tempo, old = "treatment", new = object$variable["strata"])
        }else if(!is.na(object$variable["treatment"])){
            setnames(dt.tempo, old = "treatment", new = object$variable["treatment"])
        }
        if(object$inference$se && !se){
            dt.tempo[,c("se") := NULL]
        }        
        if(all(is.na(dt.tempo$time))){
            dt.tempo[,c("time") := NULL]
        }
        setnames(dt.tempo, old = "estimate", new = "risk")
        if(object$inference$bootstrap){
            if(!estimate.boot){
                dt.tempo[,c("estimate.boot") := NULL]
            }else{
                setnames(dt.tempo, old = "estimate.boot", new = "risk.boot")
            }
        }
        if(object$inference$band && object$inference$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
            dt.tempo[,c("lowerBand") := Publish::formatCI(lower = .SD$lowerBand, upper = .SD$upperBand)]
            dt.tempo[,c("upperBand") := NULL]
            setnames(dt.tempo, old = "lowerBand", new = "simultaneous ci")
            if(quantile==FALSE){
                dt.tempo[,c("quantileBand") := NULL]
            }
        }
        if(object$inference$ci){
            dt.tempo[,c("lower") := Publish::formatCI(lower = .SD$lower, upper = .SD$upper)]
            dt.tempo[,c("upper") := NULL]
            setnames(dt.tempo, old = "lower", new = "ci")
        }

        ## display
        if(short == 1){
            cat(" - Standardized risk between time zero and 'time'\n",sep="")
            cat("\n")        
            print(dt.tempo, digits=digits, row.names = FALSE)
            cat("\n")
        }else if(short == 0){
            cat(" - Standardized risk between time zero and 'time', reported on the scale [0;1] (probability scale)\n",sep="")
            if(!is.na(object$variable["strata"])){
                cat("   (average risk within strata)\n")
            }else{ 
                cat("   (average risk when treating all subjects with one treatment)\n",sep="")
            }
            cat("\n")        
            print(dt.tempo, digits=digits, row.names = FALSE)
            cat("\n")
        }
        
        
    }else{
        ## ** difference or ratio
        ## prepare
        if("diffRisk" %in% type){
            iIndexRow <- which((object$diffRisk$estimator == estimator) * (interaction(object$diffRisk$A,object$diffRisk$B) %in% interaction(allContrasts[1,],allContrasts[2,])) == 1)
            dt.tempo <- object$diffRisk[iIndexRow,.SD, .SDcols = setdiff(names(object$diffRisk),"estimator")]
        }else if("ratioRisk" %in% type){
            iIndexRow <- which((object$ratioRisk$estimator == estimator) * (interaction(object$ratioRisk$A,object$ratioRisk$B) %in% interaction(allContrasts[1,],allContrasts[2,])) == 1)
            dt.tempo <- object$ratioRisk[iIndexRow,.SD, .SDcols = setdiff(names(object$ratioRisk),"estimator")]
        }
        if(!is.na(object$variable["treatment"])){
            setnames(dt.tempo, old = c("A","B"), new = paste(object$variable["treatment"],c("A","B"),sep="="))
        }else if(!is.na(object$variable["strata"])){
            setnames(dt.tempo, old = c("A","B"), new = paste(object$variable["strata"],c("A","B"),sep="="))
        }
        if("meanRisk" %in% type == FALSE){
            dt.tempo[,c("estimate.A","estimate.B") := NULL]
        }else{
            if(!is.na(object$variable["treatment"])){
                setnames(dt.tempo, old = c("estimate.A","estimate.B"), new = c(paste0("risk(",object$variable["treatment"],"=A)"),paste0("risk(",object$variable["treatment"],"=B)")))
            }else if(!is.na(object$variable["strata"])){
                setnames(dt.tempo, old = c("estimate.A","estimate.B"), new = c("risk(",object$variable["strata"],"=A)","risk(",object$variable["strata"],"=B)"))
            }
        }
        if("diffRisk" %in% type){
            setnames(dt.tempo, old = "estimate", new = "difference")
        }else if("ratioRisk" %in% type){
            setnames(dt.tempo, old = "estimate", new = "ratio")
        }
        if(object$inference$se && !se){
            dt.tempo[,c("se") := NULL]
        }
        if(all(is.na(dt.tempo$time))){
            dt.tempo[,c("time") := NULL]
        }
        if(object$inference$bootstrap){
            if(!estimate.boot){
                dt.tempo[,c("estimate.boot") := NULL]
            }else{
                if("diffRisk" %in% type){
                    setnames(dt.tempo, old = "estimate.boot", new = "difference.boot")
                }else if("ratioRisk" %in% type){
                    setnames(dt.tempo, old = "estimate.boot", new = "ratio.boot")
                }
            }
        }
        if(object$inference$band && object$inference$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
            dt.tempo[,c("lowerBand") := Publish::formatCI(lower = .SD$lowerBand, upper = .SD$upperBand)]
            dt.tempo[,c("upperBand") := NULL]
            setnames(dt.tempo, old = "lowerBand", new = "simultaneous ci")
            if(quantile==FALSE){
                dt.tempo[,c("quantileBand") := NULL]
            }
            if("adj.p.value" %in% names(dt.tempo)){
                setnames(dt.tempo, old = "adj.p.value", new = "adjusted p.value")
            }
        }
        if(object$inference$ci){
            dt.tempo[,c("lower") := Publish::formatCI(lower = .SD$lower, upper = .SD$upper)]
            dt.tempo[,c("upper") := NULL]
            setnames(dt.tempo, old = "lower", new = "ci")
        }

        ## display
        if("diffRisk" %in% type){
            
            if(short==0){
                cat(" - Difference in standardized risk (B-A) between time zero and 'time' \n")
                cat("                reported on the scale [-1;1] (difference between two probabilities)\n",sep="")
                if(!is.na(object$variable["strata"])){
                    cat(" (difference between the average risk in two different strata)\n")
                }else{ 
                    cat(" (difference in average risks when treating all subjects with the experimental treatment (B),\n",
                        "                               vs. treating all subjects with the reference treatment (A))\n")
                }
                cat("\n")
                print(dt.tempo,digits=digits,row.names = FALSE)
                cat("\n")
            }else if(short==1){
                cat(" - Difference in standardized risk (B-A) between time zero and 'time' \n")
                cat("\n")
                print(dt.tempo,digits=digits,row.names = FALSE)
                cat("\n")
            }
        }else if("ratioRisk" %in% type){
            if(short==0){
                cat(" - Ratio of standardized risks (B/A) between time zero and 'time' \n")
                cat("                reported on the scale ]0;1] (ratio of two probabilities)\n",sep="")
                if(!is.na(object$variable["strata"])){
                    cat(" (ratio between the average risk in two different strata)\n")
                }else{ 
                    cat(" (ratio of average risks when treating all subjects with the experimental treatment (B),\n",
                        "                          vs. treating all subjects with the reference treatment (A))\n")
                }
                cat("\n")
                print(dt.tempo,digits=digits,row.names = FALSE)
                cat("\n")
            }else if(short == 1){
                cat(" - Ratio of standardized risks (B/A) between time zero and 'time' \n")
                cat("\n")
                print(dt.tempo,digits=digits,row.names = FALSE)
                cat("\n")
            }
        }
    }

    ## CI/Band
    if(short == 0){
        if("diffRisk" %in% type){
            cat(" difference      : estimated difference in standardized risks \n")
            if(object$inference$bootstrap && estimate.boot){
                cat(" difference.boot : average value over the bootstrap samples \n")
            }
        }else if("ratioRisk" %in% type){
            cat(" ratio           : estimated ratio between standardized risks \n")
            if(object$inference$bootstrap && estimate.boot){
                cat(" ratio.boot      : average value over the bootstrap samples \n")
            }
        }else{
            cat(" risk            : estimated standardized risk \n")
            if(object$inference$bootstrap && estimate.boot){
                cat(" risk.boot       : average value over the bootstrap samples \n")
            }
        }
        if(object$inference$ci){
            cat(" ci              : pointwise confidence intervals \n")
            if(object$inference$p.value && any(c("diffRisk","ratioRisk") %in% type)){
                cat(" p.value         : (unadjusted) p-value \n")
            }
        }
    
    if(object$inference$band && object$inference$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
        if(object$inference$band==1){
            cat(" simultaneous ci : simulatenous confidence intervals over time\n")
        }else if(object$inference$band==2){
            cat(" simulatenous ci : simulatenous confidence intervals over time and treatment\n")
        }
    }
    
    if(object$inference$p.value && object$inference$band && any(c("diffRisk","ratioRisk") %in% type)){
        if(object$inference$band==1){
            cat(" adjusted p.value: p-value adjusted for multiple comparisons over time\n")
        }else if(object$inference$band==2){
            cat(" adjusted p.value: p-value adjusted for multiple comparisons over time and treatment\n")
        }
    }
    }
    
    ## export
    return(invisible(object))
}


######################################################################
### summary.ate.R ends here
