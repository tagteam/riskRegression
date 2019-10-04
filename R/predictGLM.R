### predictGLM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2018 (13:29) 
## Version: 
## Last-Updated: okt  4 2019 (15:55) 
##           By: Brice Ozenne
##     Update #: 44
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * predictGLM
##' @title Compute the influence function for the prediction.
##' @description Compute the influence function for the prediction from a linear or logistic model.
##' @name predictGLM
##'
##' @param object glm model.
##' @param newdata [data.frame] dataset containing the covariate to condition on.
##' @param average.iid [logical] Should the influence function be averaged over the empirical distribution.
##' 
##' @examples
##' \dontrun{
##' library(lava)
##' m <- lvm(Y~X1+X2+X3)
##'
##' set.seed(10)
##' d <- lava::sim(m, 1e2)
##'
##' ## check for lm
##' e.lm <- lm(Y~X1+X2+X3, data = d)
##' test <- predictGLM(e.lm, newdata = d)
##' test.av <- predictGLM(e.lm, newdata = d, average.iid = TRUE)
##'
##' GS <- lava::estimate(e.lm, f = function(p,data){
##'   p["(Intercept)"] + d[,"X1"] * p["X1"] + d[,"X2"] * p["X2"] + d[,"X3"] * p["X3"]
##' })
##' range(test[,1]-GS$coef)
##' range(attr(test,"iid")-t(GS$iid))
##' range(colMeans(attr(test,"iid"))-attr(test.av,"iid"))
##'
##' ## check for glm
##' e.glm <- glm(I(Y>0)~X1+X2+X3, data = d, family = binomial(link = "logit"))
##' test <- predictGLM(e.glm, newdata = d)
##' test.av <- predictGLM(e.glm, newdata = d, average.iid = TRUE)
##'
##' GS <- lava::estimate(e.glm, f = function(p,data){
##'   lava::expit(p["(Intercept)"] + d[,"X1"] * p["X1"] + d[,"X2"] * p["X2"] + d[,"X3"] * p["X3"])
##' })
##' range(test[,1]-GS$coef)
##' range(attr(test,"iid")-t(GS$iid))
##' range(colMeans(attr(test,"iid"))-attr(test.av,"iid"))
##' }


## * predictGLM (code)
##' @rdname predictGLM
predictGLM <- function(object, newdata, average.iid = FALSE){

    if(any(class(object) %in% c("lm","glm") == FALSE)){
        stop("Can only hand lm or glm models \n")
    }
    
    n.obs <- NROW(newdata)
    out <- cbind(predict(object, type = "response", newdata = newdata, se = FALSE))

    ## ** prepare average.iid
    if(average.iid){
        if(is.null(attr(average.iid,"factor"))){
            factor <- list(matrix(1, nrow = n.obs, ncol = 1))
        }else{
            factor <- attr(average.iid, "factor")
            if(is.matrix(factor)){
                factor <- list(factor)
            }
            if(!is.list(factor)){
                stop("Attribute \'factor\' for argument \'average.iid\' must be a list \n")
            }
            if(any(sapply(factor, is.matrix)==FALSE)){
                stop("Attribute \'factor\' for argument \'average.iid\' must be a list of matrices \n")
            }
            if(any(sapply(factor, function(iF){NROW(iF)==NROW(newdata)})==FALSE)){
                stop("Attribute \'factor\' for argument \'average.iid\' must be a list of matrices with ",NROW(newdata)," rows \n")
            }
        }
        n.factor <- NCOL(factor)
    }
    
    ## ** compute influence function of the coefficients using lava
    iid.beta <- lava::iid(object)

    ## ** chain rule
    newX <- model.matrix(stats::formula(object), newdata)

    if(identical(class(object)[1],"lm") || object$family$link[[1]]=="identity"){
        if(average.iid){
            E.X <- apply(factor, 2, function(iFactor){ ## iFactor <- factor[,1]
                colMeans(colMultiply_cpp(newX, scale = iFactor))
            })
            attr(out, "iid") <- iid.beta %*% E.X            
        }else{            
            attr(out, "iid") <- t(apply(newX, 1, function(iRow){ ## iRow <- newX[1,]
                iid.beta %*% cbind(iRow)
            }))
        }
    }else if(object$family$link=="logit"){
        ## 1/(1+exp(-Xbeta)) - risk.i
        ## newX %*% coef(object) - Xbeta
        Xbeta <- predict(object, type = "link", newdata = newdata, se=FALSE)
        if(average.iid){
            attr(out, "iid") <- lapply(factor, function(iFactor){
                iE.X <- apply(iFactor, 2, function(iiFactor){ ## iiFactor <- factor[[1]][,1]
                    colMeans(colMultiply_cpp(newX, scale = iiFactor * exp(-Xbeta)/(1+exp(-Xbeta))^2))
                })
                return(iid.beta %*% iE.X)
            })
            if(is.null(attr(average.iid,"factor"))){
                attr(out, "iid") <- attr(out, "iid")[[1]]
            }
        }else{
            attr(out, "iid") <- t(sapply(1:n.obs, function(iObs){ ## iObs <- 1
                iid.beta %*% cbind(newX[iObs,]) * exp(-Xbeta[iObs])/(1+exp(-Xbeta[iObs]))^2
            }))
        }
    }else {
        stop("Cannot handle ",object$family$link," \n",
             "Only handle the following link function: identity, logit \n")
    }

    return(out)
}

######################################################################
### predictGLM.R ends here
