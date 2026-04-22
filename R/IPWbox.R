### IPWbox.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 16 2026 (15:54) 
## Version: 
## Last-Updated: apr 22 2026 (19:32) 
##           By: Brice Ozenne
##     Update #: 192
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * IPWbox (documentation)
##' @title Encapsulate Weights
##' @description Object containing user-defined weights that can be used with the \code{ate} function.
##'
##' @param variable [character vector] name of the treatment variable or (time,censoring) variables .
##' @param proba [matrix] value of the weight for each observation (rows) and treatment group or timepoints (columns).
##' @param IF [array] value of the influence function for each observation (rows) relative to each weight (columns,slice).
##' @param indicator [matrix] indicator of the observed treatment assignement. Should have the same dimension as argument proba.
##' @param id [list of a single vector ] optional id variable.
##'
##' @return An IPWbox object
##' 
##' @details Enable to specified user-defined treatment and censoring weights with the ate function.
##' For instance if one wants to truncate the weights or restrict the population to non-extreme weights.
##' 
#' @seealso
#' \code{\link{drop1.IPWbox}} to modify the weights.
##' 
##' @examples
##'
#' library(data.table)
#' library(survival)
#' \dontrun{
#' estimator <- "all"
#' }
#' \dontshow{
#' estimator <- "IPTW"
#' }
#' #############################################
#' #### Survival settings without censoring ####
#' #### ATE with glm                        ####
#' #############################################
#' 
#' #### generate data ####
#' n <- 250
#' set.seed(2)
#' dt <- sampleData(n, outcome="binary")
#' dt[, id := paste0("id",1:.N)]
#' dt[, X2 := as.numeric(X2)]
#' 
#' #### working models ####
#' m.event <- glm(formula = Y ~ X1+X2+X3, data=dt, family = "binomial")
#' m.treatment <- glm(X3~X1+X2+X5+X8, data=dt,family=binomial(link="logit"))
#'
#' #### compute the ATE ####
#' ## A) usual
#' ateFit1 <- ate(m.event, treatment = m.treatment, 
#'                data = dt, store = c(weights = TRUE),
#'                estimator = estimator)
#' model.tables(ateFit1, estimator = "all")
#' colSums(weights(ateFit1))
#' 
#'
#' ## B) same except that the IPTW sum up to exactly NROW(dt)
#' ateFit2 <- ate(m.event, treatment = m.treatment,
#'                data = dt, store = c(weights = TRUE), data.index = 1:NROW(dt),
#'                estimator = estimator)
#' model.tables(ateFit2, estimator = "all")
#' colSums(weights(ateFit2))
#' 
#' ## C) same as version 1 but 'manually' providing weights values + influence function
#' manualWeights <- weights(ateFit1, type = "IPTWbox")
#' ateFit3 <- ate(m.event, treatment = manualWeights,
#'                data = dt, store = c(weights = TRUE),
#'                estimator = estimator)
#' model.tables(ateFit3, estimator = "all")
#' colSums(weights(ateFit3))
#'
#' \dontrun{
#' ## D) set observations with high weight to 0 weight (only affect IPTW and AIPTW not GFORMULA)
#' manualWeights.trim <- drop1(manualWeights, method = "trim", scope = 5)
#' ateFit.trim <- ate(m.event, treatment = manualWeights.trim,
#'                    data = dt, store = c(weights = TRUE),
#'                    estimator = estimator)
#' model.tables(ateFit.trim, estimator = "all")
#' colSums(weights(ateFit.trim))
#'
#' ## E) set observations with high weight to 5 weight (only affect IPTW and AIPTW not GFORMULA)
#' manualWeights.threshold <- drop1(manualWeights, method = "threshold", scope = 5)
#' ateFit.threshold <- ate(m.event, treatment = manualWeights.threshold,
#'                         data = dt, store = c(weights = TRUE),
#'                         estimator = estimator)
#' model.tables(ateFit.threshold, estimator = "all")
#' colSums(weights(ateFit.threshold))
#'
#' ## F) evaluate on a subset
#' index.rm <- which(rowSums(weights(ateFit1, type = "IPTW"))<=5)
#' dt.rm <- dt[index.rm]
#' ateFit.rm <- ate(m.event, treatment = m.treatment,
#'                  data = dt.rm, data.index = index.rm, store = c(weights = TRUE),
#'                  estimator = estimator)
#' model.tables(ateFit.rm, estimator = "all")
#' colSums(weights(ateFit.rm))
#' }

## * IPWbox (code)
##' @export
IPWbox <- function(variable, proba, IF = NULL, indicator = NULL, id = NULL){

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
    if(!is.null(indicator)){
        if(!is.matrix(indicator)){
            stop("Argument \'indicator\' should be a matrix. \n")
        }
        if(any(dim(indicator)!=dim(proba))){
            stop("Argument \'indicator\' should be a matrix with the same dimension as argument \'proba\'. \n")
        }
        if(any(indicator %in% 0:1 == FALSE)){
            stop("Argument \'indicator\' should only contain 0/1 or TRUE/FALSE values. \n")
        }
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
                indicator = indicator,
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

## * nobs.IPWbox
##' @export
print.IPWbox <- function(x, ...){
    cat("\t",x$type," for ",NROW(x$proba)," observations with respect to ",all.vars(x$formula)[1],".\n \n",sep="")

    ## ** compute IPW 
    weights <- x$indicator / x$proba
    if(any(is.na(x$proba) & x$indicator==0)){
        weights[is.na(x$proba) & x$indicator==0] <- 0
    }
    if(any(is.infinite(x$proba) & x$indicator==0)){
        weights[is.infinite(x$proba) & x$indicator==0] <- 0
    }
    ## ** print
    dfprint <- data.frame(levels = c(colnames(x$proba),"all"),
                          range = c(apply(weights, MARGIN = 2, FUN = function(iCol){paste0("[",formatC(min(iCol[iCol!=0]), format = "f", digits = 2),"; ",
                                                                                           formatC(max(iCol[iCol!=0]), format = "f", digits = 2),"]")}),
                                    paste0("[",formatC(min(weights[weights!=0]),format = "f", digits = 2),"; ",formatC(max(weights[weights!=0]), format = "f", digits = 2),"]")),                          
                          total = c(paste0(formatC(colSums(weights), format = "f", digits = 2), " (",colSums(weights==0 & x$indicator>0)," zeros)"),""))

    rownames(dfprint) <- NULL
    print(dfprint, row.names = FALSE)
    
    return(invisible(dfprint))    
    
}

## * drop1.IPWbox
##' @title Update Large IPW values
##' @description Update IPW larger than a user-defined value, either setting them to 0 or to a fixed value.
##'
##' @param object [IPWbox] output of IPWbox.
##' @param scope [numeric,>0] value beyond which the weight should be modified.
##' @param method [character] either set the weight to 0 (\code{"trim"}) by setting the corresponding probability to infinity
##' or set the weight to the value of argument scope (\code{"threshold"}).
##' @param normalize [logical] if \code{TRUE} will rescale the 'non-modified' weights (i.e. whose value is below scope)
##' to match the initial population size. The sum by columns of the weights is then unchanged between before and after the application of drop1.
##' @param ... not used. For compatibility with the generic method
##'
##' @return An IPWbox object
##'
##' @details WARNING: how to update the influence function is not well established and a pragmatic choice has been made: \itemize{
##' \item \code{"trim"}: influence function for the weights set to 0 is also set to 0.
##' \item \code{"threshold"}: influence function is unchanged
##' }
##' but there is approach has not be proven to control the type 1 error.
##' @export
drop1.IPWbox <- function(object, scope, method = "trim", normalize = TRUE, ...){
    
    ## ** compute IPW 
    weights <- object$indicator / object$proba
    if(any(is.na(object$proba) & object$indicator==0)){
        weights[is.na(object$proba) & object$indicator==0] <- 0
    }
    if(any(is.infinite(object$proba) & object$indicator==0)){
        weights[is.infinite(object$proba) & object$indicator==0] <- 0
    }
    n.obs <- NROW(weights)

    ## ** normalize user input
    ## *** method
    method <- match.arg(method, c("trim","threshold"))

    ## *** normalize
    if(normalize){
        if(is.null(object$indicator)){
            stop("Cannot normalize without knowledge of the treatment membership. \n",
                 "Consider specifying the \'indicator\' argument when calling IPWbox. \n")
        }
    }

    ## ** update object
    ls.indexDrop <- apply(object$proba, MARGIN = 2, FUN = function(iCol){which(iCol<1/scope)})
    for(iT in 1:NCOL(object$proba)){
        iIndex.drop <- ls.indexDrop[[iT]]
        if(length(iIndex.drop)>0 && method=="trim"){
            ## update estimate
            object$proba[iIndex.drop,iT] <- Inf
            ## update influence function
            if(all(lengths(object$IF)==0)){
                object$IF[[iT]][,iIndex.drop] <- 0
            }
        }else if(length(iIndex.drop)>0 && method=="threshold"){
            ## update estimate
            object$proba[iIndex.drop,iT] <- (1/scope)
            ## update influence function
            if(all(lengths(object$IF)==0)){
                ## object$IF[[iT]][,iIndex.drop] <- ???
            }
        }

    }

    ## ** compute new IPW and normalize
    ## lengths(ls.indexDrop)>0: some weights should be modified for normalization to not be 1
    ## lengths(ls.indexDrop)<n.obs: some weights should not be fixed to 0 or scope so something can be normalized
    if(normalize & any(lengths(ls.indexDrop)>0 & lengths(ls.indexDrop)<n.obs)){
        if(method == "trim"){
            ## find normalization values
            newweights <- object$indicator/object$proba
            if(any(is.na(object$proba) & object$indicator==0)){
                newweights[is.na(object$proba) & object$indicator==0] <- 0
            }
            if(any(is.infinite(object$proba) & object$indicator==0)){
                newweights[is.infinite(object$proba) & object$indicator==0] <- 0
            }
            normalize.values <- colSums(weights)/colSums(newweights)
        }else if(method == "threshold"){
            indicatorDrop <- do.call(cbind,lapply(ls.indexDrop, function(iIndex.drop){1:n.obs %in% iIndex.drop}))
            normalize.values <- (colSums(weights) - colSums(object$indicator*indicatorDrop) * scope) / colSums(weights*(indicatorDrop==FALSE))
            if(length(iIndex.drop) < NROW(object$proba) && any(weights[-iIndex.drop,iT] * normalize.values[iT] > scope)){
                message("Due to the normalization step some of the weights are greater than the value of argument \'scope\'.\n",
                        "Consider re-applying the drop1 method to the output. \n")
            }
        }

        for(iT in which(lengths(ls.indexDrop)>0 & lengths(ls.indexDrop)<n.obs)){
            iIndex.drop <- ls.indexDrop[[iT]]
            object$proba[-iIndex.drop,iT] <- object$proba[-iIndex.drop,iT] / normalize.values[iT]

            if(length(object$IF[[iT]])>0){
                object$IF[[iT]][,-iIndex.drop] <- object$IF[[iT]][,-iIndex.drop] / normalize.values[iT]
            }
        }

            
            ## update influence function
            

    }

    ## ** export
    return(object)
}


##----------------------------------------------------------------------
### IPWbox.R ends here
