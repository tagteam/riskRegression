### print.Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: May 31 2016 (11:32) 
## Version: 
## last-updated: Aug  7 2017 (09:12) 
##           By: Thomas Alexander Gerds
##     Update #: 22
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

##' Print method for risk prediction scores
##'
##' @title Print Score object
##' @param x Object obtained with \code{Score.list}
##' @param digits Number of digits
##' @param ... passed to print
#'
#' @method print Score
#' @export
print.Score <- function(x,digits=3,...){
    for (s in c(x$summary)){
        cat(paste0("\nSummary statistics ",s,":\n"))
        print(x[[s]],digits=digits, ...)
    }
    for (m in c(x$metrics)){
        cat(paste0("\nMetric ",m,":\n"))
        print(x[[m]],digits=digits, ...)
    }
    for (p in c(x$plots)){
        cat(paste0("\nData for ",p," plot are stored in the object as x[[\"",p,"\"]].\n"))
        ## print(x[[m]],digits=digits, ...)
    }
}

#' @method print scoreAUC
#' @export
print.scoreAUC <- function(x,B,digits=3,...){
    cat("\nResults by model:\n\n")
    print(x$score,digits=digits,...)
    if (length(x$contrasts)>0){
        cat("\nResults of model comparisons:\n\n")
        print(x$contrasts,digits=digits,...)
    }
}
#' @method print scoreBrier
#' @export
print.scoreBrier <- function(x,B,digits=3,...){
    cat("\nResults by model:\n\n")
    print(x$score,digits=digits,...)
    if (length(x$contrasts)>0){
        cat("\nResults of model comparisons:\n\n")
        print(x$contrasts,digits=digits,...)
    }
}


#----------------------------------------------------------------------
### print.Score.R ends here
