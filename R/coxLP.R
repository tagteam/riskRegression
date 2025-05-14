### coxLP.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:32) 
## Version: 
## Last-Updated: May 14 2025 (17:14) 
##           By: Thomas Alexander Gerds
##     Update #: 33
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * coxLP
#' @title Compute the linear predictor of a Cox model
#' @description Compute the linear predictor of a Cox model
#' @name coxLP 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param center should the linear predictor be computed after centering the covariates
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @details In case of empty linear predictor returns a vector of 0 with the same length as the number of rows of the dataset

#' @rdname coxLP
#' @export
coxLP <- function(object, data, center){
  UseMethod("coxLP") 
} 

## ** coxLP.cph
#' @rdname coxLP
#' @method coxLP cph
#' @export
coxLP.cph <- function(object, data, center){
  
    coef <- stats::coef(object)
    n.varLP <- length(coef)

    if(n.varLP==0){
        if(is.null(data)){
            Xb <- rep(0, coxN(object))
        }else{
            Xb <- rep(0, NROW(data))
        }
    }else if(is.null(data)){ ## training dataset
    
        Xb <- object$linear.predictors
    
        if(center[[1]] == FALSE && n.varLP != 0){
            Xb <- Xb + sum(coxCenter(object)*coef)
        }
    
    }else{ ## new dataset
        ## if(all(object$Design$assume %in% c("category","asis","interaction","strata"))){
        ## if(!("rcspline" %in% object$Design$assume)){
        Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
        ## }else{ ## does not work with splines
        ## X <- model.matrix(object, data = data)
        ## Xb <- as.vector(X %*% coef)
        ## }
        
      
      if(center == FALSE){
        Xb <- Xb + sum(coxCenter(object)*coef)
      }
      
  }
  
  return(unname(Xb))
}

## ** coxLP.coxph
#' @rdname coxLP
#' @method coxLP coxph
#' @export
coxLP.coxph <- function(object, data, center){
  
  coef <- stats::coef(object)
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    
    Xb <- object$linear.predictors
    
    if(center[[1]] == FALSE && n.varLP != 0){
      Xb <- Xb + sum(coxCenter(object)*coef)
    }
    
  }else{ ## new dataset
      if(n.varLP>0){
      is.strata <- attr(object$terms, "special")$strata

      
      if(length(is.strata)>0){
        object.strata <- object[["strata"]]
        object[["strata"]] <- NULL # solve a bug in survival:::predict.coxph when fitting the model with strata

        Xb <- try(rowSums(stats::predict(object, newdata = as.data.frame(data), 
                                         type = "terms")), silent = TRUE)
        if(inherits(x=Xb,what="try-error")){ ## Fix an error when the dataset used to fit the object is removed from the global environment
          ## survival:::predict.coxph search for it and read (at least) the status variable
          txt <- paste0("survival::predict.coxph returns the following error:\n",
                        as.character(Xb),
                        "It seems that the dataset used to fit the model is no more compatible with the model,\n",
                        "probably because it has been modified afterwards.\n",
                        "coxLP.coxph will try to reconstruct the original dataset and continue the execution.\n")
          warning(txt)
          ## So avoid an error, the following code re-create the original dataset
          object[["strata"]] <- object.strata
          object$call$data <- reconstructData(object)
          object[["strata"]] <- NULL
          
          Xb <- rowSums(stats::predict(object, newdata = as.data.frame(data), 
                                       type = "terms"))
        }
      }else{ 
        Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
      }
      if(center == FALSE){
        Xb <- Xb + sum(coxCenter(object)*coef)
      }
      
    }else{ 
      Xb <- rep(0, NROW(data)) 
    } 
  }

  return(unname(Xb))
}

## ** coxLP.phreg
#' @rdname coxLP
#' @method coxLP phreg
#' @export
coxLP.phreg <- function(object, data, center){

    coef <- stats::coef(object)
    n.varLP <- length(coef)

    if(n.varLP>0){
        if(is.null(data)){ ## training dataset
            X <- object$X
        }else{
            X <- model.matrix(object, data = data)
        }
  
        Xb <- as.vector(X[,names(coef), drop = FALSE] %*% coef)
    
        if(center){
            Xb <- Xb - sum(coxCenter(object)*coef)
        }
    
    }else{
        if(is.null(data)){
            Xb <- rep(0, coxN(object))
        }else{
            Xb <- rep(0, NROW(data))
        }
    } 
  
  
  return(Xb)
}


## ** coxLP.prodlim
#' @rdname coxLP
#' @method coxLP prodlim
#' @export
coxLP.prodlim <- function(object, data, center){

    if(is.null(data)){
        Xb <- rep(0, coxN(object))
    }else{
        Xb <- rep(0, NROW(data))
    }
    return(Xb)
}


## ** coxLP.GLMnet
#' @rdname coxLP
#' @method coxLP GLMnet
#' @export
coxLP.GLMnet <- function(object, data, center = FALSE){
    if(is.null(data)){
        if(inherits(object$fit,"cv.glmnet")){
            return(as.numeric(predict(object$fit$glmnet.fit, s = object$selected.lambda, newx = object$x)))
        } else{
            return(as.numeric(predict(object$fit, s = object$selected.lambda, newx = object$x)))
        }
    } else{
        ## ff <- stats::formula(stats::delete.response(prodlim::strip.terms(object$terms,specials = c("unpenalized","strata"),arguments = NULL)))
        ff <- stats::formula(stats::delete.response(object$terms))
        newdata <- Publish::specialFrame(formula = ff,
                                data = data,
                                strip.specials = c("strata","unpenalized"),
                                strip.arguments = NULL,
                                specials = c("strata","unpenalized"),
                                unspecials.design = TRUE,
                                specials.design = TRUE,
                                response = FALSE)
        if (NCOL(newdata$unpenalized)>0){
            newX <- cbind(newdata$design,newdata$unpenalized)
        }else{
            newX <- newdata$design
        }
        if(inherits(object$fit,"cv.glmnet")){
            return(as.numeric(predict(object$fit$glmnet.fit, s = object$selected.lambda, newx = newX)))
        } else{
            return(as.numeric(predict(object$fit, s = object$selected.lambda, newx = newX)))
        }
    }
}

######################################################################
### coxLP.R ends here
