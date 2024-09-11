### predictCoxPL.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  3 2024 (15:15) 
## Version: 
## Last-Updated: sep 10 2024 (12:20) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * predictCoxPL (code)
#' @title Deprecated Function for Product Limit Estimation of Survival Probabilities .
#' @description Depreciated function for Product Limit Estimation of Survival Probabilities from a Cox model.
#' Use the \code{\link{predictCox}} function instead with argument \code{product.limit=TRUE}.
#'
#' @inheritParams predictCox
#' @param ... additional arguments to be passed to \code{\link{predictCox}}.
#'
#' @export
predictCoxPL <- function(object, ...){

    message("The function predictCoxPL is deprecated. \n",
            "predictCox with argument \'product.limit = TRUE\' should be used instead. \n")
    predictCox(object, product.limit = TRUE, ...)
    
}


##----------------------------------------------------------------------
### predictCoxPL.R ends here
