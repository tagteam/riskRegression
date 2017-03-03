### as.data.table.predictCSC.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: Mar  3 2017 (21:06) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Turn predictCSC object into a \code{data.table}
#' @description Turn predictCSC object into a \code{data.table}
#' @param x object obtained with function \code{predictCSC}
#' @param keep.rownames not used
#' @param ... not used
#' @export
as.data.table.predictCSC <- function(x,keep.rownames=FALSE,...){
    times=NULL
    nd=x$newdata
    if (is.null(nd))
        return(nd)
    else{
        data.table::setDT(nd)
        out <- data.table::rbindlist(lapply(1:length(x$times),function(tt){
            ndtt=copy(nd)
            nd[,times:=x$times[[tt]]]
            if (!is.null(x$strata))
                nd[,strata:=x$strata]
            ar <- cbind(x[["absRisk"]][,tt])
            colnames(ar) <- "absRisk"
            if (x$se==1L){
                ar <- cbind(ar,x[[paste0("absRisk",".se")]][,tt],x[[paste0("absRisk",".lower")]][,tt],x[[paste0("absRisk",".upper")]][,tt])
                colnames(ar) <- paste0("absRisk",c("",".se",".lower",".upper"))
            }
            ## setDT(tyc)
            nd <- cbind(nd,ar)
            nd   
        }))    
    }
}



######################################################################
### as.data.table.predictCSC.R ends here
