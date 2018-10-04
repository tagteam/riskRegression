### as.data.table.ate.R --- 
##----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: Oct  3 2018 (16:23) 
##           By: Thomas Alexander Gerds
##     Update #: 93
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * as.data.table.ate (documentation)
#' @title Turn ate Object Into a \code{data.table}
#' @description Turn ate object into a \code{data.table}.
#' @name as.data.table.ateRobust
#' 
#' @param x object obtained with function \code{ate}
#' @param keep.rownames Not used.
#' @param ... Not used.
#'

## * as.data.table.ateRobust (code)
#' @rdname as.data.table.ateRobust
#' @export
as.data.table.ateRobust <- function(x, keep.rownames = FALSE, ...){
    name.method <- colnames(x$ate.value)
    ## if(x$se==FALSE){
    ## x$ate.se[,grep("Gformula|^IPTW",name.method)] <- NA
    ## }
    iDT.risk.0 <- data.table(method = name.method,  estimand = "risk.0", value = x$ate.value["risk.0",])
    iDT.risk.1 <- data.table(method = name.method,  estimand = "risk.1", value = x$ate.value["risk.1",])
    iDT.ate.diff <- data.table(method = name.method,  estimand = "ate.diff", value = x$ate.value["ate.diff",])

    if(x$se){
        iDT.risk.0[, c("se") := x$ate.se["risk.0",]]
        iDT.risk.1[, c("se") := x$ate.se["risk.1",]]
        iDT.ate.diff[, c("se") := x$ate.se["ate.diff",]]
    }

    if(!is.null(x$conf.level)){
        iDT.risk.0[, c("lower","upper") := list(x$ate.lower["risk.0",],x$ate.upper["risk.0",])]
        iDT.risk.1[, c("lower","upper") := list(x$ate.lower["risk.1",],x$ate.upper["risk.1",])]
        iDT.ate.diff[, c("lower","upper") := list(x$ate.lower["ate.diff",],x$ate.upper["ate.diff",])]

        iDT.risk.0[, c("p.value") := x$ate.p.value["risk.0",]]
        iDT.risk.1[, c("p.value") := x$ate.p.value["risk.1",]]
        iDT.ate.diff[, c("p.value") := x$ate.p.value["ate.diff",]]
    }

    ## ** reorder method
    out <- rbind(iDT.risk.0, iDT.risk.1, iDT.ate.diff)
    out[, c("method") := factor(.SD$method, levels = name.method)]

    ## ** export
    return(out[])
  
}



######################################################################
### as.data.table.ate.R ends here
