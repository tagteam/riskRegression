### as.data.table.predictCox.R --- 
##----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: maj 31 2018 (11:54) 
##           By: Brice Ozenne
##     Update #: 82
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * as.data.table.predictCox (documentation)
#' @title Turn predictCox object into a \code{data.table}
#' @description Turn predictCox object into a \code{data.table}
#' @name as.data.table.predictCox
#' 
#' @param x object obtained with function \code{predictCox}
#' @param keep.rownames not used.
#' @param se should standard errors/quantile for confidence bands be displayed?
#' @param ... not used.
#'
#' @details
#' The columns \code{.seBand} corresponds to the column \code{.se} inflated
#' by the ratio between the quantile for the confidence bands and the quantile for the confidence interval.#' 


## * as.data.table.predictCox (code)
#' @rdname as.data.table.predictCox
#' @export
as.data.table.predictCox <- function(x, keep.rownames = FALSE, se = TRUE,...){
    times=NULL

    n.obs <- NROW(x[[x$type[1]]])
    nd <- data.table(observation = 1:n.obs)
    if (!is.null(x$newdata)){
        nd <- cbind(nd, x$newdata)
    }
    if(is.null(x$times)){
        stop("Cannot convert to a data.table object when times is missing in object \n",
             "set the argument \'keep.time\' to TRUE when calling the predict method \n")
    }

    if(!is.matrix(x[[x$type[1]]])){ ## baseline hazard
        out <- as.data.table(x[c("times","strata",x$type)])    
    }else{
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
                    if(se){
                        tyc <- cbind(tyc,
                                     x[[paste0(name,".se")]][,tt]
                                     )
                        vec.names <- c(vec.names,".se")
                    }
                    if(!is.null(x$conf.level)){
                        tyc <- cbind(tyc,
                                     x[[paste0(name,".lower")]][,tt],
                                     x[[paste0(name,".upper")]][,tt])
                        vec.names <- c(vec.names,".lower",".upper")
                    }
                }
                if (x$band==1L){
                    if(se && !is.null(x$conf.level)){
                        tyc <- cbind(tyc,
                                     x[[paste0(name,".quantileBand")]]
                                     )
                        vec.names <- c(vec.names,".quantileBand")                    
                    }
                    if(!is.null(x$conf.level)){
                        tyc <- cbind(tyc,
                                     x[[paste0(name,".lowerBand")]][,tt],
                                     x[[paste0(name,".upperBand")]][,tt])
                        vec.names <- c(vec.names,".lowerBand",".upperBand")
                    }
                }

                colnames(tyc) <- paste0(name,vec.names)
          
                ## setDT(tyc)
                nd <- cbind(nd,tyc)
            }
            nd   
        }))    
    }
    
    return(out)
  
}



######################################################################
### as.data.table.predictCox.R ends here
