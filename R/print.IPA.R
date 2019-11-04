### print.IPA.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  4 2019 (09:07) 
## Version: 
## Last-Updated: Nov  4 2019 (12:09) 
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
##' Print method for IPA 
##'
##' @title Print IPA object
##' @param x Object obtained with \code{IPA}
##' @param digits Number of digits
##' @param percent Logical. If \code{TRUE} show percentages.
##' @param ... passed to print
##' 
##' @method print IPA
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
print.IPA <- function(x,percent=TRUE,digits=2,...){
    if (missing(digits)){
        if (percent==TRUE) digits <- 1 else digits <- 3
    }
    if(percent==TRUE){
        X <- copy(x)
        data.table::setDT(X)
        fmt <- paste0("%1.",digits[[1]],"f")
        X[,Brier:=sprintf(fmt=fmt,100*Brier)]
        X[,IPA:=sprintf(fmt=fmt,100*IPA)]
        if (match("IPA.drop",colnames(X),nomatch=0)) X[,IPA.drop:=sprintf(fmt=fmt,100*IPA.drop)]
        print(X,...)
        message("\nNOTE: Values are multiplied by 100 and given in % (use print(...,percent=FALSE) to avoid this.")
    }else{
        print(x,digits=digits,...)
    }
    message("NOTE: IPA.drop = IPA(Full model) - IPA.")
}


######################################################################
### print.IPA.R ends here
