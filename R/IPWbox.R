### IPWbox.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 16 2026 (15:54) 
## Version: 
## Last-Updated: apr 17 2026 (17:41) 
##           By: Brice Ozenne
##     Update #: 43
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * IPWbox
##' @title Encapsulate Weights
##' @description Object containing user-defined weights that can be used with the \code{ate} function.
##'
##' @param variable [character vector] name of the treatment variable or (time,censoring) variables .
##' @param proba [matrix] value of the weight for each observation (rows) and treatment group or timepoints (columns).
##' @param IF [array] value of the influence function for each observation (rows) relative to each weight (columns,slice).
##' @param id [list of a single vector ] optional id variable.
##'
##' @return An object
##'
##' @details Enable to specified user-defined treatment and censoring weights with the ate function.
##' For instance if one wants to truncate the weights or restrict the population to non-extreme weights.
##' 
##' @examples
##'
#' library(data.table)
##'
#' #############################################
#' #### Survival settings without censoring ####
#' #### ATE with glm                        ####
#' #############################################
#' 
#' #### generate data ####
#' n <- 250
#' dt <- sampleData(n, outcome="binary")
#' dt[, X2 := as.numeric(X2)]
#' 
#' #### compute the ATE ####
#' m.event <- glm(formula = Y ~ X1+X2+X3, data=dt, family = "binomial")
#' m.treatment <- glm(X3~X1+X2+X5+X8, data=dt,family=binomial(link="logit"))
#' ateFit <- ate(m.event, treatment = m.treatment, data = dt, store = c(weights = TRUE))
#' model.tables(ateFit)
#'
#' range(weights(ateFit, type = "IPTW"))
#'
#' manualWeights <- weights(ateFit, type = "IPTWbox")
#' manualWeights$proba
#' ateFit.bis <- ate(m.event, treatment = manualWeights, data = dt)
#' model.tables(ateFit.bis) ## same as when providing the logistic model
#' 
#' #### set weights above 5 to 0
#' manualWeights.trim <- drop1(manualWeights, scope = 5)
#' ateFit.trim <- ate(m.event, treatment = manualWeights.trim, data = dt)
#' model.tables(ateFit.trim)
#' 
##' @export
IPWbox <- function(variable, proba, IF = NULL, id = NULL){

    ## ** check user-input
    if(length(variable)!=1 || !is.character(variable)){
        stop("Argument \'variable\' should be character of length 1. \n")
    }
    if(!is.matrix(proba)){
        stop("Argument \'proba\' should be a matrix. \n")
    }
    if(!is.numeric(proba)){
        stop("Argument \'proba\' should be a numeric matrix. \n")
    }
    if(any(is.na(proba))){
        stop("Argument \'proba\' should not contain any missing values. \n")
    }
    if(is.null(colnames(proba))){
        stop("Argument \'proba\' should have column names. \n")
    }
    if(any(proba<0)){
        stop("Argument \'proba\' should take non-negative values. \n")
    }
    if(any(proba>1)){
        stop("Argument \'proba\' should take values above 1. \n")
    }
    if(!is.null(id)){
        if(!inherits(id,"list") && !inherits(id,"data.frame")){
            stop("Argument \'id\' should inherit be a list or a data.frame. \n")
        }
        if(length(id) != 1){
            stop("Argument \'id\' should have length 1. \n")
        }
        if(is.null(names(id))){
            stop("The first element of argument \'id\' should be named. \n")
        }
        if(!is.vector(id[[1]])){
            stop("The first element of argument \'id\' should be a vector. \n")
        }
        if(length(id[[1]])!=NROW(proba)){
            stop("The first element of argument \'id\' should be a vector of same length as the number of rows of argument \'proba\'. \n")
        }
        if(any(duplicated(id[[1]]))){
            stop("The first element of argument \'id\' should not contain any duplicated value. \n")
        }
    }

    ## ** create object
    out <- list(formula = stats::reformulate(termlabels = "0", response = variable),
                proba = proba,
                IF = IF,
                id = id,
                type = "IPTW")

    ## ** export
    class(out) <- append("IPWbox",class(out))
    return(out)
}

## * nobs.IPWbox
##' @export
nobs.IPWbox <- function(object, ...){
    NROW(object$proba)
}

## * drop1.IPWbox
##' @export
drop1.IPWbox <- function(object, scope, method = "trim", ...){
    
    ## ** normalize user input
    method <- match.arg(method, c("trim","threshold"))

    ## ** update object
    
    if(all(lengths(object$IF)==0)){
        if(method=="trim"){
            object$proba[object$proba<1/scope] <- Inf
        }else if(method=="threshold"){
            object$proba[object$proba<1/scope] <- 1/scope
        }
    }else{
        
        for(iT in 1:NCOL(object$proba)){
            iIndex.drop <- which(object$proba[,iT]<1/scope)
            if(method=="trim"){
                object$proba[iIndex.drop,iT] <- Inf
                object$IF[[iT]][,iIndex.drop] <- 0
            }else if(method=="threshold"){
                object$proba[iIndex.drop,iT] <- 1/scope
                ## object$IF[[iT]][,iIndex.drop] <- ???
            }
        }

    }

    ## ** export
    return(object)
}


##----------------------------------------------------------------------
### IPWbox.R ends here
