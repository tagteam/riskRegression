### print.Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: May 31 2016 (11:32) 
## Version: 
## last-updated: Feb 24 2017 (14:19) 
##           By: Thomas Alexander Gerds
##     Update #: 11
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
    B <- x$splitMethod$B
    for (m in c(x$summary,x$metrics,x$plots)){
        cat(paste0("\nMetric ",m,":\n"))
        print(x[[m]],B,digits=digits, ...)
    }
}


#----------------------------------------------------------------------
### print.Score.R ends here
