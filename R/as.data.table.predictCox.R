### as.data.table.predictCox.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: apr 12 2017 (15:34) 
##           By: Brice Ozenne
##     Update #: 13
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
  times=NULL

  n.obs <- NROW(x[[x$type[1]]])
  nd <- data.table(observation = 1:n.obs)
  if (!is.null(x$newdata)){
      nd <- cbind(nd, x$newdata)
  }
  out <- data.table::rbindlist(lapply(1:length(x$times),function(tt){
      ndtt=copy(nd)
      nd[,times:=x$times[[tt]]]
      if (!is.null(x$strata))
          nd[,strata:=x$strata]
      for (name in x$type){
          tyc <- cbind(x[[name]][,tt])
          colnames(tyc) <- name
          vec.names <- c("")
          if (x$se==1L){
              tyc <- cbind(tyc,
                           x[[paste0(name,".se")]][,tt],
                           x[[paste0(name,".lower")]][,tt],
                           x[[paste0(name,".upper")]][,tt])
              vec.names <- c(vec.names,".se",".lower",".upper")
          }
          if (x$band==1L){
            tyc <- cbind(tyc,
                         x[[paste0(name,".lowerBand")]][,tt],
                         x[[paste0(name,".upperBand")]][,tt])
            vec.names <- c(vec.names,".lowerBand",".upperBand")
          }
          colnames(tyc) <- paste0(name,vec.names)
          
          ## setDT(tyc)
          nd <- cbind(nd,tyc)
      }
      nd   
  }))    
 
  return(out)
  
}



######################################################################
### as.data.table.predictCox.R ends here
