### splitStrataVar.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:38) 
## Version: 
## Last-Updated: Apr 27 2025 (07:39) 
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

#' @title Reconstruct each of the strata variables
#' @description Reconstruct each of the strata variables from the strata variable stored in the coxph object.
#' @name splitStrataVar
#' 
#' @param object a coxph object.
#'
#' @author Brice Ozenne broz@@sund.ku.dk and Thomas A. Gerds tag@@biostat.ku.dk
#'
splitStrataVar <- function(object){
    grid.levels <- expand.grid(object$xlevels)
    vec.levels <- apply(grid.levels,1,paste,collapse=", ")
    ## find original name of the strata variables
    xterms <- delete.response(object$terms)
    name.strataVar.original <- all.vars(xterms)[attr(xterms,"specials")[["strata"]]]
    name.strataVar <- attr(xterms,"term.labels")[attr(xterms,"specials")[["strata"]]]
    n.strataVar <- length(name.strataVar)
    ## identify value relative to each strata variable
    df.strata <- NULL
    for(iStrataVar in 1:n.strataVar){ # iStrataVar <- 1
        value2level <- setNames(grid.levels[,name.strataVar[iStrataVar]], vec.levels)
        if(iStrataVar==1){
            df.strata <- data.frame(unname(value2level[object$strata]))
        }else{
            df.strata <- cbind(df.strata,
                               unname(value2level[object$strata]))
        }
    }
    ## export
    names(df.strata) <- name.strataVar.original
    return(df.strata)
}

######################################################################
### splitStrataVar.R ends here
