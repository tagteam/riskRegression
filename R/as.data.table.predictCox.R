### as.data.table.predictCox.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: Mar  3 2017 (09:29) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Turn predictCox object into a \code{data.table}
#' @description Turn predictCox object into a \code{data.table}
#' @param x object obtained with function \code{predictCox}
#' @param keep.rownames not used
#' @param ... not used
#' @export
as.data.table.predictCox <- function(x,keep.rownames=FALSE,...){
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
            for (name in x$type){
                tyc <- cbind(x[[name]][,tt])
                colnames(tyc) <- name
                if (x$se==1L){
                    tyc <- cbind(tyc,x[[paste0(name,".se")]][,tt],x[[paste0(name,".lower")]][,tt],x[[paste0(name,".upper")]][,tt])
                    colnames(tyc) <- paste0(name,c("",".se",".lower",".upper"))
                }
                ## setDT(tyc)
                nd <- cbind(nd,tyc)
            }
            nd   
        }))    
    }
}



######################################################################
### as.data.table.predictCox.R ends here
