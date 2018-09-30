### as.data.table.ate.R --- 
##----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: sep 19 2018 (08:17) 
##           By: Brice Ozenne
##     Update #: 88
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
#' @param se [logical] Should standard errors be displayed?
#' @param ... Not used.
#'

## * as.data.table.ateRobust (code)
#' @rdname as.data.table.ateRobust
#' @export
as.data.table.ateRobust <- function(x, keep.rownames = FALSE, se = TRUE, ...){

    name.estimator <- colnames(x$ate.value)
    if(x$se==FALSE){
        x$ate.se[,grep("Gformula|^IPTW",name.estimator)] <- NA
    }

    iDT.risk.0 <- data.table(estimator = name.estimator,  estimand = "risk.0", value = x$ate.value["risk.0",])
    iDT.risk.1 <- data.table(estimator = name.estimator,  estimand = "risk.1", value = x$ate.value["risk.1",])
    iDT.ate.diff <- data.table(estimator = name.estimator,  estimand = "ate.diff", value = x$ate.value["ate.diff",])

    if(se){
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

    ## ** reorder estimator
    out <- rbind(iDT.risk.0, iDT.risk.1, iDT.ate.diff)
    out[, c("estimator") := factor(.SD$estimator, levels = name.estimator)]

    ## ** export
    return(out)
  
}



######################################################################
### as.data.table.ate.R ends here
