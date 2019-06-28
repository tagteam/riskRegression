### predictRiskIID.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 28 2019 (14:38) 
## Version: 
## Last-Updated: jun 28 2019 (16:03) 
##           By: Brice Ozenne
##     Update #: 28
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * predictRiskIID 
predictRiskIID <- function(object, newdata, average.iid, ...){
  UseMethod("predictRiskIID",object)
}

## * predictRiskIID.default
predictRiskIID.default <- function(object,
                                   newdata,
                                   average.iid,
                                   ...){
    stop(paste0("No method available for extracting the iid decomposition for the predictions from objects in class: ",class(object),".
                 But, you can write it yourself or ask the package manager."),call.=FALSE)
}

## * predictRiskIID.glm
predictRiskIID.glm <- function(object,
                               newdata,
                               average.iid,
                               factor){
    
    if (object$family$family=="binomial"){

        if(!is.matrix(factor)){
            stop("Argument \'factor\' must be a matrix. \n")
        }
        if(NROW(factor) != stats::nobs(object)){
            stop("Argument \'factor\' must have the same number of rows as the dataset used to fit the object. \n")
        }

        attr(average.iid, "factor") <- factor
        resPred <- predictGLM(object,
                              newdata = newdata,
                              average.iid = average.iid)
        return(attr(resPred,"iid"))
    }else{
        stop("Currently only the binomial family is implemented for extracting the iid decomposition of the predictions from a glm object.")
    }
    
}

## * predictRiskIID."cox"
predictRiskIID.coxph <- function(object, newdata, times, average.iid, factor = NULL, ...){
    if(!is.null(factor)){
        n.times <- length(times)
        n.obs <- NROW(newdata)
        if(!is.matrix(factor)){
            stop("Argument \'factor\' must be a matrix. \n")
        }
        if(NROW(factor) != coxN(object)){
            stop("Argument \'factor\' must have the same number of rows as the dataset used to fit the object. \n")
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
                          type = "survival",
                          ...)
    ## there is a minus because predictor.cox return iid(survival) = -iid(risk)

    if(average.iid){
        if(!is.null(factor)){
            return(lapply(resPred$survival.average.iid, function(x){-x}))
        }else{
            return(-resPred$survival.average.iid)
        }
    }else{
        return(-resPred$survival.iid)
    }
    
}

predictRisk.cph <- predictRisk.coxph
predictRisk.phreg <- predictRisk.coxph

## * predictRiskIID.CSC
predictRiskIID.CauseSpecificCox <- function(object,
                                            newdata,
                                            times,
                                            average.iid,
                                            factor = NULL,
                                            ...){

    if(!is.null(factor)){
        n.times <- length(times)
        n.obs <- NROW(newdata)
        n.cause <- length(object$cause)
        if(!is.matrix(factor)){
            stop("Argument \'factor\' must be a matrix. \n")
        }
        if(NROW(factor) != coxN(object)[1]){
            stop("Argument \'factor\' must have the same number of rows as the dataset used to fit the object. \n")
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
        return(resPred$absRisk.average.iid)
    }else{
        return(resPred$absRisk.iid)
    }
}

######################################################################
### predictRiskIID.R ends here
