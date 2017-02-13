### print.predictCSC.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 11 2017 (10:01) 
## Version: 
## last-updated: Feb 11 2017 (12:24) 
##           By: Thomas Alexander Gerds
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @method print predictCSC
#' @export
print.predictCSC <- function(x,...){
    absRisk <- x$absRisk
    if (is.null(x$absRisk.se)){
        colnames(absRisk) <- x$times
        print(absRisk)
    }
    else{
        lower <- pmax(0,absRisk + qnorm((1-x$conf.level)/2) * x$absRisk.se)
        upper <- pmin(1,absRisk - qnorm((1-x$conf.level)/2) * x$absRisk.se)
        colnames(lower) <- colnames(upper) <- colnames(absRisk) <- x$times
        print.listof(list("absolute risk"= absRisk,"lower"=lower,"upper"=upper))
    }
}


#----------------------------------------------------------------------
### print.predictCSC.R ends here
