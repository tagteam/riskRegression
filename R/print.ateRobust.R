### print.ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2018 (15:28) 
## Version: 
## Last-Updated: Jan 29 2019 (10:49) 
##           By: Thomas Alexander Gerds
##     Update #: 104
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
print.ateRobust <- function(x,
                            digits = 3,
                            augment.cens = x$augment.cens,
                            ...){

    p.value=estimand=treatment=lower=upper=NULL
    dt.x <- as.data.table(x)
    keep.col <- c("method","value")

    ## ** add confidence interval and p.value
    if(!is.null(x$se[[1]])&& x$se[[1]]==TRUE){
        ## dt.x[,c("lower") := format(round(.SD$lower, digits = digits), digits = 3, nsmall = 3)]
        ## dt.x[,c("upper") := format(round(.SD$upper, digits = digits), digits = 3, nsmall = 3)]
        indexNNA <- which(!is.na(dt.x$se))
        dt.x[,c("[lower ; upper]") := ""]
        dt.x[indexNNA,sprintf(fmt=paste0("[%1.",digits,"f;%1.",digits,"f]"),lower,upper)]
        dt.x[indexNNA,c("[lower ; upper]") := sprintf(fmt=paste0("[%1.",digits,"f;%1.",digits,"f]"),lower,upper)]
        ## lowerValue <- 10^(-digits)
        dt.x[,p.value := format.pval(p.value,eps=10^{-digits},digits=digits,na.form="")]
        ## sapply(p.value, function(iP){
        ## if(is.na(iP)){""}else if(iP < lowerValue){paste0("<",formatC(lowerValue))}else{format(iP, digits = digits, nsmall = digits)}
        keep.col <- c(keep.col, "[lower ; upper]","p.value")
    }

    ## ** prepare for display
    ## use keep.col instead of "value" because the last
    ## column can either be "value" or "value [CI inf ; CI sup]"
    dt.printMeanRisk <- rbind(cbind(dt.x[estimand=="risk.0",.SD,.SDcols = setdiff(keep.col, "p.value")],treatment = x$level.treatment[1]),
                              cbind(dt.x[estimand=="risk.1",.SD,.SDcols = setdiff(keep.col, "p.value")],treatment = x$level.treatment[2]))
    data.table::setkeyv(dt.printMeanRisk, cols = c("method","treatment"))
    data.table::setcolorder(dt.printMeanRisk, c(keep.col[1],"treatment",setdiff(keep.col[-1], "p.value"))) 
    dt.printMeanRisk[treatment==x$level.treatment[2], c("method") := ""]

    dt.printDiffRisk <- dt.x[estimand=="ate.diff",.SD,.SDcols = keep.col]

    ## ** display
    cat("\nAverage risk:\n")
    print(dt.printMeanRisk, row.names = FALSE,digits=digits)
    cat("Average risks of the event between time zero and 'time'")
    cat("are standardized to covariate distribution of all subjects\nare given on probability scale [0,1] in hypothetical worlds\nin which all subjects are treated with one of the treatment options\nstandardized to all subjects.\n\n")

    cat("Risk difference:\n\n")
    print(dt.printDiffRisk, row.names = FALSE,digits=digits)
    cat("\nComparisons of risks on probability scale [0,1] between hypothetical worlds are given as 
 treatment value ",x$level.treatment[2]," minus ",x$level.treatment[1],"
and are interpreted as what would have been observed had treatment been randomized.\n\n",sep="")
    ## output
    return(invisible(x))
}



######################################################################
### print.ateRobust.R ends here
