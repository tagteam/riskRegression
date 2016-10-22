### print.Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: May 31 2016 (11:32) 
## Version: 
## last-updated: Oct  2 2016 (09:10) 
##           By: Thomas Alexander Gerds
##     Update #: 9
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
#' @export
##' @param x Object obtained with \code{Score.list}
##' @param digits Number of digits
##' @param ... passed to print
print.Score <- function(x,digits=3,...){
    B <- x$splitMethod$B
    for (m in c(x$summary,x$metrics)){
        cat(paste0("\nMetric ",m,":\n"))
        print(x[[m]],B,digits=digits, ...)
    }
}

##' Print metric specific element of risk prediction assessment
##'
##' @title Print metric specific results of risk prediction assessment
#' @export
##' @param x Element of the result of \code{Score.list}
##' @param B Number of splits if any.
##' @param digits Number of digits
##' @param ... passed to print
print.highscore <- function(x,B=0,digits=3,...){
    p=model=reference=times=NULL
    if (B>0){
        cat("\nScores\n\n")
        print(x$score,digits=digits)
        cat(paste0("\nCross-validation (average of ",B," steps)\n\n"))
        if (!is.null(x$contrasts)){
            cat("\nTests\n\n")
            if (match("times",colnames(x$contrasts),nomatch=0))
                print(x$contrasts[,list("p-value (median)"=median(p)),by=list(model,reference,times)],digits=digits,...)
            else
                print(x$contrasts[,list("p-value (median)"=median(p)),by=list(model,reference)],digits=digits,...)
            cat(paste0("\nMultisplit test (",B," splits)\n\n"))
        }
    }else{
        cat("\nScores\n\n")
        print(x$score,digits=digits)
        cat("\nTests\n\n")
        if (!is.null(x$contrasts) && NROW(x$contrasts)>0)
            print(x$contrasts,digits=digits,...)
    }
}


#----------------------------------------------------------------------
### print.Score.R ends here
