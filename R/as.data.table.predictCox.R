### as.data.table.predictCox.R --- 
##----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar  3 2017 (09:28) 
## Version: 
## Last-Updated: Mar  3 2025 (13:00) 
##           By: Thomas Alexander Gerds
##     Update #: 178
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * as.data.table.predictCox (documentation)
#' @title Turn predictCox Object Into a \code{data.table}
#' @description Turn predictCox object into a \code{data.table}.
#' @name as.data.table.predictCox
#' @param x object obtained with function \code{predictCox}
#' @param keep.rownames Not used.
#' @param se [logical] Should standard errors/quantile for confidence bands be displayed?
#' @param ... Not used.
## * as.data.table.predictCox (code)
#' @rdname as.data.table.predictCox
#' @export
as.data.table.predictCox <- function(x, keep.rownames = FALSE, se = TRUE,...){
    times=NULL
    
    n.obs <- NROW(x[[x$type[1]]])
    if (!is.null(x$status)){
        nd <- x$status[,.SD,.SDcols = c("nevent","strata")]
    }else if(x$baseline){
        nd <- data.table::data.table(observation = 1:n.obs, strata = x$strata)
    }else{
        nd <- data.table::data.table(observation = 1:n.obs, x$newdata, strata = x$strata)
    }
    
    if(is.null(x$times)){
        stop("Cannot convert to a data.table object when times is missing in object \n",
             "set the argument \'keep.times\' to TRUE when calling the predict method \n")
    }
    if(!is.matrix(x[[x$type[1]]])){ ## baseline hazard
        out <- as.data.table(x[c("times",x$type)])
        if (!is.null(nd)){
            out <- cbind(nd,out)
        }
    }else{        
        if(x$diag || identical(as.character(x$type),"lp")){
            n.times <- 1
        }else{
            n.times <- length(x$times)
        }        
        out <- data.table::rbindlist(lapply(1:n.times,function(tt){
            ndtt=data.table::copy(nd)
            if(x$diag && length(x$times)>0){
                nd[,times:=x$times]
            }else if(length(x$times)>0){
                nd[,times:=x$times[[tt]]]
            }
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
                    if(!is.null(x[[paste0(name,".transform")]])){
                        tyc <- cbind(tyc,
                                     x[[paste0(name,".lower")]][,tt],
                                     x[[paste0(name,".upper")]][,tt])
                        vec.names <- c(vec.names,".lower",".upper")
                    }
                }
                if (x$band==1L){
                    if(se[[1]] && !is.null(x[[paste0(name,".transform")]])){
                        tyc <- cbind(tyc,
                                     x[[paste0(name,".quantileBand")]]
                                     )
                        vec.names <- c(vec.names,".quantileBand")                    
                    }
                    if(!is.null(x[[paste0(name,".transform")]])){
                        tyc <- cbind(tyc,
                                     x[[paste0(name,".lowerBand")]][,tt],
                                     x[[paste0(name,".upperBand")]][,tt])
                        vec.names <- c(vec.names,".lowerBand",".upperBand")
                    }
                }
                colnames(tyc) <- paste0(name,vec.names)
          
                nd <- cbind(nd,tyc)
            }
            nd   
        }))
    }
    
    return(out)
  
}



######################################################################
### as.data.table.predictCox.R ends here
