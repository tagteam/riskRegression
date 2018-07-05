### iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jul  5 2018 (13:29) 
## Version: 
## Last-Updated: jul  5 2018 (14:57) 
##           By: Brice Ozenne
##     Update #: 16
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .predictGLM
#' @title Compute the influence function for the prediction.
#' @description Compute the influence function for the prediction from a linear or logistic model.
#' @name .predictGLM
#' 
#' @examples
#' \dontrun{
#' library(lava)
#' m <- lvm(Y~X1+X2+X3)
#'
#' set.seed(10)
#' d <- lava::sim(m, 1e2)
#'
#' ## check for lm
#' e.lm <- lm(Y~X1+X2+X3, data = d)
#' test <- .predictGLM(e.lm, newdata = d)
#'
#' GS <- lava::estimate(e.lm, f = function(p,data){
#'   p["(Intercept)"] + d[,"X1"] * p["X1"] + d[,"X2"] * p["X2"] + d[,"X3"] * p["X3"]
#' })
#' range(test[,1]-GS$coef)
#' range(attr(test,"iid")-t(GS$iid))
#'
#' ## check for glm
#' e.glm <- glm(I(Y>0)~X1+X2+X3, data = d, family = binomial(link = "logit"))
#' test <- .predictGLM(e.glm, newdata = d)
#'
#' GS <- lava::estimate(e.glm, f = function(p,data){
#'   lava::expit(p["(Intercept)"] + d[,"X1"] * p["X1"] + d[,"X2"] * p["X2"] + d[,"X3"] * p["X3"])
#' })
#' range(test[,1]-GS$coef)
#' range(attr(test,"iid")-t(GS$iid))
#' }


## * .predictGLM (code)
#' @rdname .predictGLM
.predictGLM <- function(object, newdata){

    if(any(class(object) %in% c("lm","glm") == FALSE)){
        stop("Can only hand lm or glm models \n")
    }
    
    n.obs <- NROW(newdata)
    out <- cbind(predict(object, type = "response", newdata = newdata, se = FALSE))
         
    ## ** compute influence function of the coefficients using lava
    iid.beta <- lava::iid(object)

    ## ** chain rule
    newX <- model.matrix(stats::formula(object), newdata)
    if(identical(class(object),"lm") || object$family$link=="identity"){
        attr(out, "iid") <- t(apply(newX, 1, function(iRow){ ## iRow <- newX[1,]
            iid.beta %*% cbind(iRow)
        }))
    }else if(object$family$link=="logit"){
        ## 1/(1+exp(-Xbeta)) - risk.i
        ## newX %*% coef(object) - Xbeta
        Xbeta <- predict(object, type = "link", newdata = newdata, se=FALSE)
        attr(out, "iid") <- t(sapply(1:n.obs, function(iObs){ ## iObs <- 1
            iid.beta %*% cbind(newX[iObs,]) * exp(-Xbeta[iObs])/(1+exp(-Xbeta[iObs]))^2
        }))                    
    }else {
        stop("Cannot handle ",object$family$link," \n",
             "Only handle the following link function: identity, logit \n")
    }

    ## ** export
    return(out)
}

######################################################################
### iid.R ends here
