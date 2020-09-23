### summary.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 29 2019 (13:18) 
## Version: 
## Last-Updated: sep 23 2020 (17:28) 
##           By: Brice Ozenne
##     Update #: 157
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
#' Can be \code{"risk"} to display the risks specific to each treatment group,
#' \code{"difference"} to display the difference in risks between treatment groups,
#' or \code{"ratio"} to display the ratio of risks between treatment groups,.
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
                        type = c("risk","difference"), se = object$inference$se, quantile = FALSE, estimate.boot = TRUE,
                        digits = 3,
                        ...){
    short <- list(...)$short ## called by the print function
    lower <- upper <- lowerBand <- upperBand <- NULL ## [:CRANcheck:] data.table

    if(!is.null(object$allContrasts)){
        allContrasts <- object$allContrasts
        contrasts <- attr(allContrasts,"contrasts")
    }else{
        contrasts <- object$contrasts
        allContrasts <- utils::combn(contrasts, m = 2)
    }

    test.se <- !is.null(object$se) && object$se
    test.ci <- !is.null(object$ci) && object$ci
    test.band <- !is.null(object$band) && object$band
    test.p.value <- !is.null(object$p.value) && object$p.value
    
    ## ** check arguments
    type <- match.arg(type,c("risk","difference","ratio"), several.ok = TRUE)
    if("difference" %in% type && "ratio" %in% type){
        stop("Cannot simulatenously display the risk difference and the risk ratio \n")
    }
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

    ## ** prepare
    txt.alternative <- switch(object$inference$alternative,
                              "two.sided" = "two-sided tests, alternative: unequal risk",
                              "greater" = "one-sided tests, alternative: greater risk",
                              "less" = "one-sided tests, alternative: less risk")



    ## ** Display: specification of ate
    print(object, display.results = FALSE)
    if(object$inference$ci || object$inference$band || object$inference$p.value){
        index.transform <- which(object$transform!="none")
        if(length(index.transform)>0){
            name.tempo <- c("mean risk","risk difference","risk ratio")
            cat("   (on the ",paste(object$transform[index.transform],collapse="/")," scale for the ",paste(name.tempo[index.transform],collapse="/"),") \n",sep="")
        }
        cat(" - Confidence level     : ", object$inference$conf.level," (",txt.alternative,")\n", sep ="")
    }
    if(object$inference$band){
        if(object$inference$method.band=="maxT-simulation"){
            text.band <- paste("- Confidence bands/adj. p-values  : single-step maxT computed using ", object$inference$n.sim," simulations\n", sep ="")
        }else if(object$inference$method.band=="maxT-integration"){
            text.band <- paste("- Confidence bands/adj. p-values  : single-step maxT computed using numerical integration\n", sep ="")
        }else if(object$inference$method.band=="bonferroni"){
            text.band <- paste("- Confidence bands/adj. p-values  : Bonferroni adjustment\n", sep ="")
        }else{
            text.band <- paste("- adj. p-values  : ",object$inference$method.band," method\n", sep ="")
        }
    }
    
    cat("\n Results: \n")
    if(identical(type,"risk")){
        ## ** only risks
        ## prepare
        iIndexRow <- which((object$meanRisk$estimator == estimator) * (object$meanRisk$treatment %in% contrasts) == 1)

        dt.tempo <- object$meanRisk[iIndexRow,.SD, .SDcols = setdiff(names(object$meanRisk),"estimator")]
        if(!is.na(object$variable["strata"])){
            setnames(dt.tempo, old = "treatment", new = object$variable["strata"])
        }else if(!is.na(object$variable["treatment"])){
            setnames(dt.tempo, old = "treatment", new = object$variable["treatment"])
        }
        if(!se){
            dt.tempo[,c("se") := NULL]
        }
        if(all(is.na(dt.tempo$time))){
            dt.tempo[,c("time") := NULL]
        }
        if(object$inference$B>0 && !estimate.boot){
           dt.tempo[,c("estimate.boot") := NULL]
        }
        if(object$inference$ci){
            dt.tempo[,c("lower") := Publish::formatCI(lower = .SD$lower, upper = .SD$upper)]
            dt.tempo[,c("upper") := NULL]
            setnames(dt.tempo, old = "lower", new = "ci")
        }
        if(object$inference$band && object$inference$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
            dt.tempo[,c("lowerBand") := Publish::formatCI(lower = .SD$lowerBand, upper = .SD$upperBand)]
            dt.tempo[,c("upperBand") := NULL]
            setnames(dt.tempo, old = "lowerBand", new = "band")
            if(quantile==FALSE){
                dt.tempo[,c("quantileBand") := NULL]
            }
        }

        ## display
        cat(" - Risk between time zero and 'time', reported on the scale [0;1] (probability scale)\n",sep="")
        if(!is.na(object$variable["strata"])){
            cat("   (average within strata)\n")
        }else{ 
            cat("   (average in hypothetical worlds in which all subjects are treated with one of the treatment options)\n")
        }
        cat("\n")        
        print(dt.tempo,digits=digits,row.names = FALSE,...)
        cat("\n")
    }else{
        ## ** difference or ratio
        ## prepare
        if("difference" %in% type){
            iIndexRow <- which((object$diffRisk$estimator == estimator) * (interaction(object$diffRisk$A,object$diffRisk$B) %in% interaction(allContrasts[1,],allContrasts[2,])) == 1)
            dt.tempo <- object$diffRisk[iIndexRow,.SD, .SDcols = setdiff(names(object$diffRisk),"estimator")]
        }else if("ratio" %in% type){
            iIndexRow <- which((object$ratioRisk$estimator == estimator) * (interaction(object$ratioRisk$A,object$ratioRisk$B) %in% interaction(allContrasts[1,],allContrasts[2,])) == 1)
            dt.tempo <- object$ratioRisk[iIndexRow,.SD, .SDcols = setdiff(names(object$ratioRisk),"estimator")]
        }
        if("risk" %in% type == FALSE){
            dt.tempo[,c("estimate.A","estimate.B") := NULL]
        }
        if(!se){
            dt.tempo[,c("se") := NULL]
        }
        if(all(is.na(dt.tempo$time))){
            dt.tempo[,c("time") := NULL]
        }
        if(object$inference$B>0 && !estimate.boot){
            dt.tempo[,c("estimate.boot") := NULL]
        }
        if(object$inference$ci){
            dt.tempo[,c("lower") := Publish::formatCI(lower = .SD$lower, upper = .SD$upper)]
            dt.tempo[,c("upper") := NULL]
            setnames(dt.tempo, old = "lower", new = "ci")
        }
        if(object$inference$band && object$inference$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
            dt.tempo[,c("lowerBand") := Publish::formatCI(lower = .SD$lowerBand, upper = .SD$upperBand)]
            dt.tempo[,c("upperBand") := NULL]
            setnames(dt.tempo, old = "lowerBand", new = "band")
            if(quantile==FALSE){
                dt.tempo[,c("quantileBand") := NULL]
            }
        }

        ## display
        if("difference" %in% type){
            cat(" - Difference in risk (B-A) between time zero and 'time' \n",
                "                reported on the scale [-1;1] (difference between two probabilities)\n",sep="")
            if(!is.na(object$variable["strata"])){
                cat(" (difference between two strata average)\n")
            }else{ 
                cat(" (average risks difference between a hypothetical world where all subjects are treated with one treatment option (A),\n",
                    "                         and another hypothetical world where all subjects are treated with another treatment option (B))\n")
            }
        }else if("ratio" %in% type){
            cat(" - Ratio of risks (B/A) between time zero and 'time' \n",
                "                reported on the scale ]0;1] (ratio of two probabilities)\n",sep="")
            if(!is.na(object$variable["strata"])){
                cat(" (ratio between two strata average)\n")
            }else{ 
                cat(" (ratio of the average risks between a hypothetical world where all subjects are treated with one treatment option (A),\n",
                    "                           and another hypothetical world where all subjects are treated with another treatment option (B))\n")
            }
        }
        cat("\n")
        print(dt.tempo,digits=digits,row.names = FALSE,...)
        cat("\n")
    }

    ## CI/Band
    if(object$inference$B>0 && estimate.boot){
        cat(" estimate.boot: average value over the bootstrap samples \n")
    }
    if(object$inference$ci){
        cat(" ci           : pointwise confidence intervals \n")
    }
    if(object$inference$p.value && any(c("difference","ratio") %in% type)){
        cat(" p.value      : unadjusted p-value \n")
    }
    
    if(object$inference$band && object$inference$method.band %in% c("bonferroni","maxT-integration","maxT-simulation")){
        if(object$inference$band==1){
            cat(" band         : simulatenous confidence intervals over time\n")
        }else if(object$inference$band==2){
            cat(" band         : simulatenous confidence intervals over time and treatment\n")
        }
    }
    
    if(object$inference$p.value && object$inference$band && any(c("difference","ratio") %in% type)){
        if(object$inference$band==1){
            cat(" adj.p.value  : p-value adjusted for multiple comparisons over time\n")
        }else if(object$inference$band==2){
            cat(" adj.p.value  : p-value adjusted for multiple comparisons over time and treatment\n")
        }
    }
    
    
    ## export
    return(invisible(object))
}


######################################################################
### summary.ate.R ends here
