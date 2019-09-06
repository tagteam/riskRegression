### predictRiskIID.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 28 2019 (14:38) 
## Version: 
## Last-Updated: sep  6 2019 (11:12) 
##           By: Brice Ozenne
##     Update #: 64
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - predictRiskIID
#' @title  Extract corrected i.i.d. decomposition
#' @description  Extract corrected i.i.d. decomposition8 from a gaussian linear model.
#' @name predictRiskIID
#'
#' @param object A fitted model from which to extract i.i.d. decomposition of the fitted values.
#' @param newdata A data frame containing predictor variable combinations.
#' @param times A vector of times in the range of the response variable,
#' e.g. for which the cumulative incidences event probabilities are computed.
#' Disregarded when argument \code{object} is a logistic regression.
#' @param average.iid Should the i.i.d. decomposition be averaged over individual observations in argument \code{newdata}?
#' @param factor  When \code{average.iid} is TRUE, enables to perform weighted averages instead of averages.
#' @param ... for compatibility with the generic method.
#'
#' @details Argument \code{factor} must be a matrix with the same number of rows as argument \code{newdata}.
#' For given column, the i.i.d. decomposition is multiplied by the column (i.e. each individual contribution is weighted) before taking the average.
#' So the first element in the first column corresponds to the weight given to the first individual.
#' Therefore for each column, an averaged i.i.d. decomposition is returned.
#' Note that a one-column matrix with only 1 is equivalent to an average, i.e. not specifying the argument \code{factor}.
#' 
#' @seealso \code{\link{iidCox}} to obtain the iid decomposition for Cox models.
#'
#' @return \code{average.iid} is FALSE: an array containing
#' the influence of each individual from the training dataset at each time on each prediction.
#' \code{average.iid} is TRUE: a list where each element contains a "weighted" averaged i.i.d. decomposition.
#' @export
predictRiskIID <- function(object, newdata, average.iid, factor, ...){
  UseMethod("predictRiskIID",object)
}

## * predictRiskIID.default
#' @rdname predictRiskIID
#' @export
predictRiskIID.default <- function(object,
                                   newdata,
                                   average.iid,
                                   factor,
                                   ...){
    stop(paste0("No method available for extracting the iid decomposition for the predictions from objects in class: ",class(object),".
                 But, you can write it yourself or ask the package manager."),call.=FALSE)
}

## * predictRiskIID.glm
#' @rdname predictRiskIID
#' @export
predictRiskIID.glm <- function(object,
                               newdata,
                               average.iid,
                               factor = NULL,
                               ...){
    
    if (object$family$family=="binomial"){

        if(!is.null(factor)){
            n.obs <- stats::nobs(object)
            if(!is.matrix(factor)){
                stop("Argument \'factor\' must be a matrix. \n")
            }
            if(NROW(factor) != NROW(newdata)){
                stop("Argument \'factor\' must have the same number of rows as argument \'newdata\'. \n")
            }
            attr(average.iid, "factor") <- factor
        }
        
        resPred <- predictGLM(object,
                              newdata = newdata,
                              average.iid = average.iid)

        ## hidden argument: enable to ask for the prediction of Y==1 or Y==0
        level <- list(...)$level
        if(!is.null(level)){
            matching.Ylevel <- table(object$data[[all.vars(formula(object))[1]]],
                                     object$y)
            all.levels <- rownames(matching.Ylevel)
            level <- match.arg(level, all.levels)

            index.level <- which(matching.Ylevel[level,]>0)
            if(length(index.level)==2){
                stop("Unknown value for the outcome variable \n")
            }else if(index.level == 1){
                attr(resPred,"iid") <- - attr(resPred,"iid")
            }
        }
        ## convert to list format for compatibility with the other predictRiskIID
        if(!is.null(factor) && NCOL(factor)>1){
            return(lapply(1:NCOL(factor), function(iF){attr(resPred,"iid")[,iF,drop=FALSE]}))
        }else{
            return(list(attr(resPred,"iid")))
        }
        
    }else{
        stop("Currently only the binomial family is implemented for extracting the iid decomposition of the predictions from a glm object.")
    }
    
}

## * predictRiskIID."cox"
#' @rdname predictRiskIID
#' @export
predictRiskIID.coxph <- function(object, newdata, average.iid, factor = NULL, times, ...){
    if(!is.null(factor)){
        n.times <- length(times)
        n.obs <- NROW(newdata)
        if(!is.matrix(factor)){
            stop("Argument \'factor\' must be a matrix. \n")
        }
        if(NROW(factor) != NROW(newdata)){
            stop("Argument \'factor\' must have the same number of rows as argument \'newdata\'. \n")
        }
        attr(average.iid, "factor") <- sapply(1:NCOL(factor),function(iCol){
            list(matrix(factor[,iCol,drop=FALSE], nrow = n.obs, ncol = n.times, byrow = FALSE))
        })
    }
    
    resPred <- predictCox(object,
                          newdata = newdata,
                          times = times,
                          iid = !average.iid,
                          average.iid = average.iid,
                          type = "survival")
    ## there is a minus because predictor.cox return iid(survival) = -iid(risk)

    if(average.iid){
        if(!is.null(factor)){
            return(lapply(resPred$survival.average.iid, function(x){-x}))
        }else{
            return(list(-resPred$survival.average.iid))
        }
    }else{
        return(-resPred$survival.iid)
    }    
}

#' @rdname predictRiskIID
#' @export
predictRiskIID.cph <- predictRiskIID.coxph

#' @rdname predictRiskIID
#' @export
predictRiskIID.phreg <- predictRiskIID.coxph

## * predictRiskIID.CSC
#' @rdname predictRiskIID
#' @export
predictRiskIID.CauseSpecificCox <- function(object,
                                            newdata,
                                            average.iid,
                                            factor = NULL,
                                            times,
                                            ...){

    if(!is.null(factor)){
        n.times <- length(times)
        n.obs <- NROW(newdata)
        n.cause <- length(object$cause)
        if(!is.matrix(factor)){
            stop("Argument \'factor\' must be a matrix. \n")
        }
        if(NROW(factor) != NROW(newdata)){
            stop("Argument \'factor\' must have the same number of rows as argument \'newdata\'. \n")
        }
        attr(average.iid, "factor") <- factor
    }
    resPred <- predict(object,
                       newdata = newdata,
                       times = times,
                       iid = !average.iid,
                       average.iid = average.iid,
                       ...)
    if(average.iid){
        if(!is.null(factor)){
            return(resPred$absRisk.average.iid)
        }else{
            return(list(resPred$absRisk.average.iid))
        }
    }else{
        return(resPred$absRisk.iid)
    }
}

## * add existing methods to riskRegression.options
riskRegression.options(method.predictRiskIID = as.character(utils::methods("predictRiskIID")))


######################################################################
### predictRiskIID.R ends here
