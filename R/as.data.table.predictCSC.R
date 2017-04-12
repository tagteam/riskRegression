### as.data.table.predictCSC.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: apr 12 2017 (16:06) 
##           By: Brice Ozenne
##     Update #: 15
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
  n.obs <- NROW(x[["absRisk"]])
  nd <- data.table(observation = 1:n.obs)
  if (!is.null(x$newdata)){
    nd <- cbind(nd, x$newdata)
  }
        
  out <- data.table::rbindlist(lapply(1:length(x$times),function(tt){
    ndtt=copy(nd)
    nd[,times:=x$times[[tt]]]
    if (!is.null(x$strata))
      nd[,strata:=x$strata]
    ar <- cbind(x[["absRisk"]][,tt])
    colnames(ar) <- "absRisk"
    if (x$se==1L){
      ar <- cbind(ar,
                  x[["absRisk.se"]][,tt],
                  x[["absRisk.lower"]][,tt],
                  x[["absRisk.upper"]][,tt])
      colnames(ar) <- paste0("absRisk",c("",".se",".lower",".upper"))
    }
    ## setDT(tyc)
    nd <- cbind(nd,ar)
    nd   
  }))    
  
  return(out)
  
}



######################################################################
### as.data.table.predictCSC.R ends here
