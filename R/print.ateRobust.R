### print.ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2018 (15:28) 
## Version: 
## Last-Updated: sep 19 2018 (08:20) 
##           By: Brice Ozenne
##     Update #: 64
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.ateRobust (documentation)
#' @title Print Average Treatment Effect
#' @description Print average treatment effect.
#' @name print.ateRobust
#' 
#' @param x object obtained with the function \code{ateRobust}.
#' @param digits [integer, >0] indicating the number of decimal places.
#' @param augment.cens [logical] should the standard errors account for the augmentation term
#' in direction of the model for the censoring mechanism.
#' accounting for the uncertainty in the outcome and propensity score models be output
#' @param ... Passed to print.
#' 


## * print.ateRobust (code)
#' @rdname print.ateRobust
#' @method print ateRobust
#' @export
print.ateRobust <- function(x, digits = 3, augment.cens = x$augment.cens,...){

    dt.x <- as.data.table(x)
    keep.col <- c("estimator","value")

    ## ** round
    dt.x[,c("value") := format(round(.SD$value, digits = digits), digit = 3, nsmall = 3)]

    ## ** add confidence interval and p.value
    if(!is.null(x$conf.level)){
        dt.x[,c("lower") := format(round(.SD$lower, digits = digits), digit = 3, nsmall = 3)]
        dt.x[,c("upper") := format(round(.SD$upper, digits = digits), digit = 3, nsmall = 3)]

        indexNNA <- which(!is.na(dt.x$se))
        dt.x[,c("[lower ; upper]") := ""]
        dt.x[indexNNA,c("[lower ; upper]") := paste("[",.SD$lower ," ; ",.SD$upper ,"]", sep ="")]

        lowerValue <- 10^(-digits)
        dt.x[,p.value := sapply(p.value, function(iP){
            if(is.na(iP)){""}else if(iP < lowerValue){paste0("<",formatC(lowerValue))}else{format(iP, digits = digits, nsmall = digits)}
        })]
        
        keep.col <- c(keep.col, "[lower ; upper]","p.value")
    }

    ## ** prepare for display
    ## use keep.col instead of "value" because the last column can either be "value" or "value [CI inf ; CI sup]"
    dt.printMeanRisk <- rbind(cbind(dt.x[dt.x$estimand=="risk.0",.SD,.SDcols = setdiff(keep.col, "p.value")],treatment = x$level.treatment[1]),
                              cbind(dt.x[dt.x$estimand=="risk.1",.SD,.SDcols = setdiff(keep.col, "p.value")],treatment = x$level.treatment[2]))
    setkeyv(dt.printMeanRisk, cols = c("estimator","treatment"))
    setcolorder(dt.printMeanRisk, c(keep.col[1],"treatment",setdiff(keep.col[-1], "p.value"))) 
    dt.printMeanRisk[dt.printMeanRisk$treatment==x$level.treatment[2], c("estimator") := ""]

    dt.printDiffRisk <- dt.x[dt.x$estimand=="ate.diff",.SD,.SDcols = keep.col]

    ## ** display
    cat("Mean risks on probability scale [0,1] in hypothetical worlds\nin which all subjects are treated with one of the treatment options:\n\n")
    print(dt.printMeanRisk, row.names = FALSE)

    cat("\nComparison of risks on probability scale [0,1] between\nhypothetical worlds are interpreted as if the treatment was randomized:\n\n")
    cat("     * risk difference \n\n")
    print(dt.printDiffRisk, row.names = FALSE)

    ## output
    return(invisible(x))
}



######################################################################
### print.ateRobust.R ends here
