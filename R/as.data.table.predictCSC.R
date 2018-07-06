### as.data.table.predictCSC.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: jul  6 2018 (13:59) 
##           By: Brice Ozenne
##     Update #: 47
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * as.data.table.predictCSC (documentation)
#' @title Turn predictCSC Object Into a \code{data.table}
#' @description Turn predictCSC object into a \code{data.table}.
#' 
#' @param x object obtained with function \code{predictCSC}
#' @param keep.rownames not used
#' @param se should standard errors/quantile for confidence bands be displayed?
#' @param ... not used

## * as.data.table.predictCSC (code)
#' @rdname as.data.table.predictCSC
#' @export
as.data.table.predictCSC <- function(x, keep.rownames = FALSE, se = TRUE, ...){
    times=NULL
    
    n.obs <- NROW(x[["absRisk"]])
    nd <- data.table(observation = 1:n.obs)
    if (!is.null(x$newdata)){
        nd <- cbind(nd, x$newdata)
    }
    n.times <- NCOL(x$absRisk) 
        
    out <- data.table::rbindlist(lapply(1:n.times,function(tt){
        ndtt=copy(nd)
        if(x$keep.times){
            nd[,times:=x$times[tt]]
        }
        if (!is.null(x$strata)){
            nd[,strata:=x$strata]
        }
        ar <- cbind(absRisk = x[["absRisk"]][,tt])
    
    vec.names <- c("absRisk")

    if (x$se==1L){
        if(se){
            ar <- cbind(ar,
                        absRisk.se=x[["absRisk.se"]][,tt])
            
        }
        if(!is.null(x$conf.level)){
            ar <- cbind(ar,
                        absRisk.lower=x[["absRisk.lower"]][,tt],
                        absRisk.upper=x[["absRisk.upper"]][,tt]
                        )
        }
    }

        if (x$band==1L){
            if(!is.null(x$conf.level)){
                ar <- cbind(ar,
                            absRisk.quantileBand=x[["absRisk.quantileBand"]],
                            absRisk.lowerBand=x[["absRisk.lowerBand"]][,tt],
                            absRisk.upperBand=x[["absRisk.upperBand"]][,tt])
            }
        }
        ## setDT(tyc)
        nd <- cbind(nd,ar)
        nd   
    }))    
  
  return(out)
  
}



######################################################################
### as.data.table.predictCSC.R ends here
