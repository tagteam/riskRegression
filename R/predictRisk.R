### predictRisk.R ---
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (09:02)
## Version:
## last-updated: May 14 2025 (15:33) 
##           By: Thomas Alexander Gerds
##     Update #: 619
#----------------------------------------------------------------------
##
### Commentary:
#
# methods for extracting predictions from various objects
#
### Change Log:
#----------------------------------------------------------------------
##
### Code:
# --------------------------------------------------------------------
#' @title Extrating predicting risks from regression models
#'
#' @description Extract event probabilities from fitted regression models and machine learning objects.
#' The function predictRisk is a generic function, meaning that it invokes
#' specifically designed functions depending on the 'class' of the first
#' argument. See \code{\link{predictRisk}}.
#' @name predictRisk
#'
#' @aliases predictRisk predictRisk.CauseSpecificCox
#' predictRisk.riskRegression predictRisk.FGR
#' predictRisk.prodlim predictRisk.rfsrc predictRisk.aalen
#' predictRisk.riskRegression predictRisk.ARR predictRisk.cox.aalen
#' predictRisk.coxph predictRisk.cph predictRisk.default
#' predictRisk.matrix predictRisk.pecCtree
#' predictRisk.pecCforest predictRisk.prodlim predictRisk.psm
#' predictRisk.selectCox predictRisk.survfit predictRisk.randomForest
#' predictRisk.lrm predictRisk.glm
#' predictRisk.rpart predictRisk.gbm
#' predictRisk.flexsurvreg
#'
#' @param object A fitted model from which to extract predicted event
#' probabilities.
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted event probabilities.
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
#' @param diag when \code{FALSE} the hazard/cumlative hazard/survival for all observations at all times is computed,
#' otherwise it is only computed for the i-th observation at the i-th time.
#' @param iid Should the iid decomposition be output using an attribute?
#' @param average.iid Should the average iid decomposition be output using an attribute?
#' @param product.limit If \code{TRUE} the survival is computed using the product limit estimator.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
#' @param truncate If \code{TRUE} truncates the predicted risks to be in the range [0, 1]. For now only implemented for the Cause Specific Cox model. 
#' @param \dots Additional arguments that are passed on to the current method.
#'
#' @return For binary outcome a vector with predicted risks. For survival outcome with and without
#' competing risks
#' a matrix with as many rows as \code{NROW(newdata)} and as many
#' columns as \code{length(times)}. Each entry is a probability and in
#' rows the values should be increasing.
#'
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#'
#' @details
#' In uncensored binary outcome data there is no need to choose a time point.
#'
#' When operating on models for survival analysis (without competing risks) the function still
#' predicts the risk, as 1 - S(t|X) where S(t|X) is survival chance of a subject characterized
#' by X.
#'
#' When there are competing risks (and the data are right censored) one needs
#' to specify both the time horizon for prediction (can be a vector) and the
#' cause of the event. The function then extracts the absolute risks F_c(t|X)
#' aka the cumulative incidence of an event of type/cause c until time t for a
#' subject characterized by X. Depending on the model it may or not be possible
#' to predict the risk of all causes in a competing risks setting. For example. a
#' cause-specific Cox (CSC) object allows to predict both cases whereas a Fine-Gray regression
#' model (FGR) is specific to one of the causes.
#'
#' @keywords survival
#'
#' @examples
#' ## binary outcome
#' library(rms)
#' set.seed(7)
#' d <- sampleData(80,outcome="binary")
#' nd <- sampleData(80,outcome="binary")
#' fit <- lrm(Y~X1+X8,data=d)
#' predictRisk(fit,newdata=nd)
#' 
#' ## survival outcome
#' # generate survival data
#' library(prodlim)
#' set.seed(100)
#' d <- sampleData(100,outcome="survival")
#' d[,X1:=as.numeric(as.character(X1))]
#' d[,X2:=as.numeric(as.character(X2))]
#' # then fit a Cox model
#' library(rms)
#' cphmodel <- cph(Surv(time,event)~X1+X2,data=d,surv=TRUE,x=TRUE,y=TRUE)
#' # or via survival
#' library(survival)
#' coxphmodel <- coxph(Surv(time,event)~X1+X2,data=d,x=TRUE,y=TRUE)
#'
#' # Extract predicted survival probabilities
#' # at selected time-points:
#' ttt <- quantile(d$time)
#' # for selected predictor values:
#' ndat <- data.frame(X1=c(0.25,0.25,-0.05,0.05),X2=c(0,1,0,1))
#' # as follows
#' predictRisk(cphmodel,newdata=ndat,times=ttt)
#' predictRisk(coxphmodel,newdata=ndat,times=ttt)
#'
#' ## simulate learning and validation data
#' set.seed(10)
#' learndat <- sampleData(80,outcome="survival")
#' valdat <- sampleData(10,outcome="survival")
#' ## use the learning data to fit a Cox model
#' library(survival)
#' fitCox <- coxph(Surv(time,event)~X6+X2,data=learndat,x=TRUE,y=TRUE)
#' ## suppose we want to predict the survival probabilities for all subjects
#' ## in the validation data at the following time points:
#' ## 0, 1, 2, 3, 4
#' psurv <- predictRisk(fitCox,newdata=valdat,times=seq(0,4,1))
#' ## This is a matrix with event probabilities (1-survival)
#' ## one column for each of the 5 time points
#' ## one row for each validation set individual
#' 
#' ## competing risks
#' library(survival)
#' library(riskRegression)
#' library(prodlim)
#' set.seed(8)
#' train <- sampleData(80)
#' test <- sampleData(10)
#' cox.fit  <- CSC(Hist(time,event)~X1+X6,data=train,cause=1)
#' predictRisk(cox.fit,newdata=test,times=seq(1:10),cause=1)
#'
#' ## with strata
#' cox.fit2  <- CSC(list(Hist(time,event)~strata(X1)+X6,
#'                       Hist(time,cause)~X1+X6),data=train)
#' predictRisk(cox.fit2,newdata=test,times=seq(1:10),cause=1)
#'
#' @export
predictRisk <- function(object,newdata,...){
    UseMethod(generic = "predictRisk",object = object)
}

## * predictRisk.default
##' @export
#' @rdname predictRisk
#' @method predictRisk default
predictRisk.default <- function(object,newdata,times,cause,...){
    stop(paste0("No method available for evaluating predicted probabilities from objects in class: ",paste0(class(object),collapse = ", "),". But, you can write it yourself or ask the package manager."),call.=FALSE)
}

## * predictRisk.double
##' @export
##' @rdname predictRisk
##' @method predictRisk double
predictRisk.double <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

## * predictRisk.integer
##' @export
##' @rdname predictRisk
##' @method predictRisk integer
predictRisk.integer <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

## * predictRisk.factor
##' @export
##' @rdname predictRisk
##' @method predictRisk factor
predictRisk.factor <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    as.numeric(object)
}

## * predictRisk.numeric
##' @export
##' @rdname predictRisk
##' @method predictRisk numeric
predictRisk.numeric <- function(object,newdata,times,cause,...){
    stopifnot(NROW(object)==NROW(newdata))
    object
}

## * predictRisk.glm
##' @export
##' @rdname predictRisk
##' @method predictRisk glm
predictRisk.glm <- function(object, newdata, iid = FALSE, average.iid = FALSE,...){

    dots <- list(...)

    if (object$family$family=="binomial"){

        n.obs <- NROW(newdata)
        out <- predict(object, type = "response", newdata = newdata, se = FALSE)

        if(iid || average.iid){
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

            ## ** compute influence function of the coefficients using wglm (more accurate than lava)
            ## iid.beta <- lava::iid(object)
            wobject <- list(times = "1", fit = list("1"=object))
            class(wobject) <- "wglm"
            object.score <- lava::score(wobject, times = "1", simplify = FALSE, indiv = TRUE)
            object.info <- lava::information(wobject, times = "1", simplify = FALSE)
            iid.beta <- object.score[["1"]] %*% solve(object.info[["1"]])

            ff.rhs <- stats::delete.response(stats::terms(stats::formula(object)))
            newX <- model.matrix(ff.rhs, newdata)
            Xbeta <- predict(object, type = "link", newdata = newdata, se = FALSE)

            if(!identical(colnames(iid.beta),colnames(newX))){
                warning("Mismatch between the column names of the design matrix and the iid. \n",
                        "Difference: \"",paste(union(setdiff(colnames(iid.beta),colnames(newX)), setdiff(colnames(newX),colnames(iid.beta))),collapse = "\" \""),"\" \n",
                        "Could be due to a factor variable with an empty level. \n")
                newX <- newX[,colnames(iid.beta),drop=FALSE]
            }

            ## ** chain rule
            if(average.iid){
                attr(out,"average.iid") <- lapply(factor, function(iFactor){
                    iE.X <- apply(iFactor, 2, function(iiFactor){ ## iiFactor <- factor[[1]][,1]
                        colMeans(colMultiply_cpp(newX, scale = iiFactor * exp(-Xbeta)/(1+exp(-Xbeta))^2))
                    })
                    return(iid.beta %*% iE.X)
                })
                if(is.null(attr(average.iid,"factor"))){
                    attr(out,"average.iid") <- attr(out,"average.iid")[[1]]
                }
            }
            if(iid){
                attr(out,"iid") <- iid.beta %*% t(colMultiply_cpp(newX, scale = exp(-Xbeta)/(1+exp(-Xbeta))^2))
            }
        }

        ## ** set correct level
        ## hidden argument: enable to ask for the prediction of Y==1 or Y==0
        level <- dots$level
        type <- dots$type ## hidden argument for ate
        if(identical(type,"survival")){
            stop("Unkown argument \'type\' for predictRisk.glm: use argument \'level\' instead. \n")
        }

        if(!is.null(level)){
            matching.Ylevel <- table(object$data[[all.vars(formula(object))[1]]],
                                     object$y)
            all.levels <- rownames(matching.Ylevel)
            level <- match.arg(level, all.levels)

            index.level <- which(matching.Ylevel[level,]>0)
            if(length(index.level) > 1){
                stop("Unknown value for the outcome variable \n")
            }else if(index.level == 1){
                out <- 1 - out
                if(iid){
                    attr(out,"iid") <- - attr(out,"iid")
                }
                if(average.iid){
                    if(is.list(attr(out,"average.iid"))){
                        attr(out,"average.iid") <- lapply(attr(out,"average.iid"), function(iIID){-iIID})
                        names(attr(out,"average.iid")) <- names(factor)
                    }else{
                        attr(out,"average.iid") <- - attr(out,"average.iid")
                    }
                }
            }
        }
        ## print(sum(abs(attr(out,"average.iid")[[1]])))
        return(out)
    } else {
        stop("Currently only the binomial family is implemented for predicting a status from a glm object.")
    }
}

## * predictRisk.multinom
##' @export
##' @rdname predictRisk
##' @method predictRisk multinom
predictRisk.multinom <- function(object, newdata, iid = FALSE, average.iid = FALSE, cause = NULL, ...){
    n.obs <- NROW(newdata)
    n.class <- length(object$lev)
    n.coef <- length(coef(object))
    n.coefperY <- n.coef/(n.class-1)
    out <- predict(object, newdata = newdata, type = "probs")

    if(n.class == 2){
        out <- cbind(1-out, out)
        colnames(out) <- object$lev
    }

    ## ** set correct level
    ## hidden argument: enable to ask for the prediction of a specific level
    level <- list(...)$level
    if(!is.null(cause)){
        if(is.null(level)){
            level <- cause
        }else if(!identical(cause,level)){
            stop("Argument \'cause\' and argument \'level\' should take the same value. \n")
        }
    }
    if(!is.null(level)){
        if(length(level)>1){
            stop("Argument \'level\' must have length 1 \n")
        }
        if(level %in% object$lev == FALSE){
            stop("Argument \'level\' must be one of \"",paste0(object$lev,collapse="\" \""),"\" \n")
        }
        out <- out[,level]
    }else if(iid || average.iid){
        stop("Argument \'level\' must be specified when exporting the iid decomposition")
    }

    if(iid || average.iid){
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

        ## ** compute influence function of the coefficients
        oldX <- model.matrix(object)
        beta <- coef(object)
        oldY <- object$fitted.values + object$residuals
        if(n.class == 2){
            beta <- matrix(beta, nrow = 1)
            oldY <- cbind(1-oldY,oldY)
        }
        oldXbeta <- oldX %*% t(beta)
        oldeXbeta <- exp(oldXbeta)
        ## range(out - colScale_cpp(cbind(1,exp(Xbeta)), scale = 1+rowSums(exp(Xbeta))))

        ## > Likelihood
        ## \prod_i \prod_k p_k(X_i)^y_{ki}
        ## > Log-likelihood
        ## \sum_i \sum_k y_{ki}\log(p_k(X_i)) = \sum_i \sum_k>1 y_{ki}\log(p_k(X_i)) + (1-\sum_k>1 y_{ki})\log(p_1(X_i))
        ##                                     = \sum_i \sum_k>1 y_{ki}\log(p_k(X_i)/p_1(X_i)) + log(p_1(X_i))
        ##                                     = \sum_i \sum_k>1 y_{ki} X_i \beta_k - log(1 + \sum_k>1 exp(X_i \beta_k))
        ## LL <- prod(out^Y)
        ## ll <- rowSums(Y[,-1,drop=FALSE] * Xbeta) - log(1 + rowSums(eXbeta))
        ## > Score
        ## /\beta_kj = \sum_i y_{ki} X_ij - X_ij exp(X_i \beta_k) / (1 + \sum_k>1 exp(X_i \beta_k))
        score <- do.call(cbind,lapply(2:n.class, function(iY){colMultiply_cpp(oldX, scale = (oldY[,iY,drop=FALSE] - oldeXbeta[,iY-1,drop=FALSE] / (1+rowSums(oldeXbeta))))}))
        ## > information
        informationM1 <- stats::vcov(object)

        iid.beta <- score %*% informationM1
        iid.beta.level <- lapply(2:n.class, function(iClass){iid.beta[,(iClass-2) * n.coefperY + 1:n.coefperY,drop=FALSE]})
        names(iid.beta.level) <- object$lev[-1]

        ## ** chain rule
        newX <- model.matrix(stats::formula(object), newdata)
        Xbeta <- cbind(0,newX %*% t(beta))
        eXbeta_rowSum <- rowSums(exp(Xbeta))
        colnames(Xbeta) <- object$lev

        if(average.iid){
            attr(out,"average.iid") <- lapply(factor, function(iFactor){
                apply(iFactor, 2, function(iiFactor){ ## exp(XB_j)/(1+\sum_k exp(XB_j))

                    if(level != object$lev[1]){
                        ## if level l which is not the reference level:  exp(XB_l)/(1+\sum_j exp(XB_j))
                        ## derivative of the numerator:  IF_l X exp(XB_l)/(1+\sum_j exp(XB_j))
                        ## derivative of the denumerator: - \sum_k IF_k X exp(XB_k+XB_l)/(1+\sum_j exp(XB_j))^2
                        iRes <- Reduce("+",lapply(object$lev[-1], function(iClass){ ## iClass <- "1"
                            - iid.beta.level[[iClass]] %*% colMeans(colMultiply_cpp(newX, scale = iiFactor * exp(Xbeta[,iClass]+Xbeta[,level])/eXbeta_rowSum^2))
                        }))
                        iRes <- iRes + iid.beta.level[[level]] %*% colMeans(colMultiply_cpp(newX, scale = iiFactor * exp(Xbeta[,level])/eXbeta_rowSum))
                    }else{
                        ## if reference level: 1/(1+\sum_j exp(XB_j))
                        ## derivative of the numerator: 0 because the numerator is 1
                        ## derivative of the denumerator: - \sum_k IF_k X exp(XB_k)/(1+\sum_j exp(XB_j))^2
                        iRes <- Reduce("+",lapply(object$lev[-1], function(iClass){ ## iClass <- "1"
                            - iid.beta.level[[iClass]] %*% colMeans(colMultiply_cpp(newX, scale = iiFactor * exp(Xbeta[,iClass])/eXbeta_rowSum^2))
                        }))
                    }

                    return(iRes)
                })
            })

            if(is.null(attr(average.iid,"factor"))){
                attr(out,"average.iid") <- attr(out,"average.iid")[[1]]
            }
        }
        if(iid){
            attr(out,"iid") <- Reduce("+",lapply(object$lev[-1], function(iClass){ ## iClass <- "1"
                - iid.beta.level[[iClass]] %*% t(colMultiply_cpp(newX, scale = exp(Xbeta[,iClass]+Xbeta[,level])/eXbeta_rowSum^2))
            }))
            if(which(level %in% object$lev)>0){
                attr(out,"iid") <- attr(out,"iid") + iid.beta.level[[level]] %*% t(colMultiply_cpp(newX, scale = exp(Xbeta[,level])/eXbeta_rowSum))
            }
        }
    }

    ## ** export
    return(out)
}
## * predictRisk.formula
##' @export
##' @rdname predictRisk
##' @method predictRisk formula
predictRisk.formula <- function(object,newdata,...){
    ff <- update.formula(object,"NULL~.")
    if (length(all.vars(ff))==1){
        p <- stats::model.frame(ff,newdata)[[1]]
        p
    } else{
        fit <- glm(object,data=newdata,family="binomial")
        predictRisk(fit,newdata=newdata,...)
    }
}

## * predictRisk.BinaryTree
##' @export
##' @rdname predictRisk
##' @method predictRisk BinaryTree
predictRisk.BinaryTree <- function(object,newdata,...){
    requireNamespace("party")
    treeresponse <- party::treeresponse
    sapply(treeresponse(object,newdata=newdata),function(x)x[1])
}

## * predictRisk.lrm
##' @export
##' @rdname predictRisk
##' @method predictRisk lrm
predictRisk.lrm <- function(object,newdata,...){
  as.numeric(stats::predict(object,newdata=newdata,type="fitted"))
}

## * predictRisk.rpart
##' @export
##' @rdname predictRisk
##' @method predictRisk rpart
predictRisk.rpart <- function(object,newdata,...){
  requireNamespace("rpart",quietly=FALSE)
  p <- as.numeric(stats::predict(object,newdata=newdata))
  p
}

## * predictRisk.randomForest
##' @export
##' @rdname predictRisk
##' @method predictRisk randomForest
predictRisk.randomForest <- function(object,newdata,...){
  as.numeric(stats::predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
}

## * predictRisk.matrix
##' @export
##' @rdname predictRisk
##' @method predictRisk matrix
predictRisk.matrix <- function(object,newdata,times,cause,...){
    if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
        stop(paste("Prediction matrix has wrong dimensions: ",
                   NROW(object),
                   " rows and ",
                   NCOL(object),
                   " columns.\n But requested are predicted probabilities for ",
                   NROW(newdata),
                   " subjects (rows) in newdata and ",
                   length(times),
                   " time points (columns)",
                   sep=""))
    }
    object
}

## * predictRisk.aalen
##' @export
##' @rdname predictRisk
##' @method predictRisk aalen
predictRisk.aalen <- function(object,newdata,times,...){
  if (is.null(object$B.iid)) {
    stop("Need resample.iid=1 when fitting the aalen model to make",
         " predictions. Please re-fit the aalen model with resample.idd=1.")
  }
  
  out <- timereg::predict.aalen(object=object,
                                newdata=newdata,
                                times=times,
                                se=FALSE,
                                ...)$S0
  return(1 - out)
}

## * predictRisk.cox.aalen
##' @export
##' @rdname predictRisk
##' @method predictRisk cox.aalen
predictRisk.cox.aalen <- function(object,newdata,times,...){
    out <- timereg::predict.cox.aalen(object=object,
                                      newdata=newdata,
                                      times=times,
                                      se=FALSE,
                                      ...)$S0
    return(1 - out)
}


## * predictRisk.comprisk
##' @export
##' @rdname predictRisk
##' @method predictRisk comprisk
predictRisk.comprisk <- function(object, newdata, times, ...) {
  out <- timereg::predict.comprisk(object=object,
                                   newdata=newdata,
                                   times=times,
                                   se=FALSE,
                                   ...)$P1
  return(out)
}


## * predictRisk.coxph
##' @export
##' @rdname predictRisk
##' @method predictRisk coxph
predictRisk.coxph <- function(object,
                              newdata,
                              times,
                              product.limit = FALSE,
                              diag = FALSE,
                              iid = FALSE,
                              average.iid = FALSE,
                              ...){
    dots <- list(...)
    type <- dots$type ## hidden argument for ate
    outPred <- predictCox(object=object,
                          newdata=newdata,
                          times=times,
                          se = FALSE,
                          iid = iid,
                          diag = diag,
                          average.iid = average.iid,
                          store = dots$store, ## hidden argument
                          product.limit = product.limit,
                          type="survival")

    if(identical(type,"survival")){
        out <- outPred$survival
    }else{
        out <- 1-outPred$survival
    }
    if(iid){
        if(identical(type,"survival")){
            attr(out,"iid") <- outPred$survival.iid
        }else{
            attr(out,"iid") <- -outPred$survival.iid
        }
    }
    if(average.iid){
        if(identical(type,"survival")){
            attr(out,"average.iid") <- outPred$survival.average.iid
        }else{
            if(is.list(outPred$survival.average.iid)){
                attr(out,"average.iid") <- lapply(outPred$survival.average.iid, function(iIID){-iIID})
            }else{
                attr(out,"average.iid") <- -outPred$survival.average.iid
            }
        }
    }

    return(out)
}



## * predictRisk.coxph.penal
##' @export
##' @rdname predictRisk
##' @method predictRisk coxph.penal
predictRisk.coxph.penal <- function(object,newdata,times,...){
    # assume that only one cluster/sparse penalty is allowed 
    frailhistory <- object$history[[1]]$history
    if (length(frailhistory) == 0){
        predictRisk.coxph(object,newdata,times,...)
    }else{
        frailVar <- frailhistory[NROW(frailhistory),1]
        linearPred <- predict(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
        basehaz <- basehaz(object)
        bhTimes <- basehaz[,2]
        bhValues <- basehaz[,1]
        survPred <- do.call("rbind",lapply(1:NROW(newdata),function(i){
            (1+frailVar*bhValues*exp(linearPred[i]))^{-1/frailVar}
        }))
        where <- prodlim::sindex(jump.times=bhTimes,eval.times=times)
        p <- cbind(1,survPred)[,where+1,drop = FALSE]
        if ((miss.time <- (length(times) - NCOL(p)))>0)
            p <- cbind(p,matrix(rep(NA,miss.time*NROW(p)),nrow=NROW(p)))
        if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
        1-p
    }
}


## * predictRisk.cph
##' @export
##' @rdname predictRisk
##' @method predictRisk cph
predictRisk.cph <- predictRisk.coxph

## * predictRisk.selectCox
##' @export
##' @rdname predictRisk
##' @method predictRisk selectCox
predictRisk.selectCox <- function(object,newdata,times,...){
    predictRisk(object[[1]],newdata=newdata,times=times,...)
}


## * predictRisk.prodlim
##' @export
##' @rdname predictRisk
##' @method predictRisk prodlim
predictRisk.prodlim <- function(object,
                                newdata,
                                times,
                                cause,
                                diag = FALSE,
                                iid = FALSE,
                                average.iid = FALSE,...){
    ## require(prodlim)
    if (object$model[[1]]=="competing.risks" && missing(cause)){
        stop(paste0("Cause is missing. Should be one of the following values: ",paste(attr(object$model.response,"states"),collapse=", ")))
    }

    if(diag == FALSE && iid == FALSE && average.iid == FALSE){
        p <- predict(object=object,
                     cause=cause,
                     type="cuminc",
                     newdata=newdata,
                     times=times,
                     mode="matrix",
                     level.chaos=1)
        ## if the model has no covariates
        ## then all cases get the same prediction
        ## in this exceptional case we return a vector
        if (NROW(p)==1 && NROW(newdata)>=1)
            p <- as.vector(p)
        ## p[is.na(p)] <- 0
        if (is.null(dim(p))){
            if (length(p)!=length(times))
                stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
        }
        else{
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
        }

    }else{ ## NOTE: can compute prediction at individual specific times (argument diag = TRUE, e.g. time 1 for id=1 and time 5 for id = 2)
        ##       can also return the influence function
        ## but does not handle competing risk
        dots <- list(...)
        type <- dots$type ## hidden argument    
        product.limit <- dots$product.limit ## hidden argument    
        iPred <- predictCox(object, diag = diag, newdata = newdata, times = times,
                            type = "survival", product.limit = product.limit, iid = iid, average.iid = average.iid,
                            store = dots$store) ## hidden argument

        if(is.null(type) || "risk" == type){
            p <- 1-iPred$survival
            if(iid){
                attr(p,"iid") <- -iPred$survival.iid
            }
            if(average.iid){
                if(is.list(iPred$survival.average.iid)){
                    attr(p,"average.iid") <- lapply(iPred$survival.average.iid, function(iM){-iM})
                }else{
                    attr(p,"average.iid") <- -iPred$survival.average.iid
                }
            }
        }else if("survival" == type){
            p <- iPred$survival
            if(iid){
                attr(p,"iid") <- iPred$survival.iid
            }
            if(average.iid){
                attr(p,"average.iid") <- iPred$survival.average.iid
            }
        }else{
            stop("Incorrect argument \'type\': should be \"survival\" or \"risk\". \n")
        }
    }

    ## ** export
    return(p)

}

## * predictRisk.survfit
##' @rdname predictRisk
##' @method predictRisk survfit
##' @export
predictRisk.survfit <- function(object,newdata,times,...){
    p <- predict.survfit(object,newdata=newdata,times=times,type="cuminc",bytimes=TRUE,fill="last")
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predict.survfit <- function(object,newdata,times,bytimes=TRUE,type="cuminc",fill="last",...){
    if (length(class(object))!=1 || class(object)[[1]]!="survfit" || object$type[[1]] !="right")
        stop("Predictions only available \nfor class 'survfit', possibly stratified Kaplan-Meier fits.\n For class 'cph' Cox models see survest.cph.")
    if (missing(newdata))
        npat <- 1
    else
        if (is.data.frame(newdata))
            npat <- nrow(newdata)
    else stop("If argument `newdata' is supplied it must be a dataframe." )
    ntimes <- length(times)
    maxtime <- max(object$time)
    sfit <- summary(object,times=pmin(times,maxtime))
    if (is.na(fill))
        Fill <- function(x,len){x[1:len]}
    else if (fill=="last")
        Fill <- function(x,len){
            y <- x[1:len]
            y[is.na(y)] <- x[length(x)]
            y}
    else stop("Argument fill must be the string 'last' or NA.")
    if (is.null(object$strata)){
        pp <- Fill(sfit$surv,ntimes)
        p <- matrix(rep(pp,npat),
                    ncol=ifelse(bytimes,ntimes,npat),
                    nrow=ifelse(bytimes,npat,ntimes),
                    byrow=bytimes)
    }
    else{
        covars <- attr(terms(eval.parent(object$call$formula)),"term.labels")
        if (!all(match(covars,names(newdata),nomatch=FALSE)))
            stop("Not all strata defining variables occur in newdata.")
        ## FIXME there are different ways to build strata levels
        ## how can we test which one was used???
        stratdat <- lapply(covars,function(x)newdata[[x]])
        names(stratdat) <- covars
        NewStratVerb <- survival::strata(stratdat)
        NewStrat <- interaction(stratdat,sep=" ")
        levs <- levels(sfit$strata)
        #    print(levs)
        #    print(levels(NewStrat))
        #    print(levels(NewStratVerb))
        if (!all(choose <- match(NewStratVerb,levs,nomatch=F))
            &&
            !all(choose <- match(NewStrat,levs,nomatch=F)))
            stop("Not all strata levels in newdata occur in fit.")
        survlist <- split(sfit$surv,sfit$strata)
        pp <- lapply(survlist[choose],Fill,ntimes)
        p <- matrix(unlist(pp,use.names=FALSE),
                    ncol=ifelse(bytimes,ntimes,npat),
                    nrow=ifelse(bytimes,npat,ntimes),
                    byrow=bytimes)
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    if (type=="cuminc") 1-p else p
}


## * predictRisk.psm
##' @export
##' @rdname predictRisk
##' @method predictRisk psm
predictRisk.psm <- function(object,newdata,times,...){
    if (length(times)==1){
        p <- rms::survest(object,times=c(0,times),newdata=newdata,what="survival",conf.int=FALSE)[,2]
    }else{
        p <- rms::survest(object,times=times,newdata=newdata,what="survival",conf.int=FALSE)
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    1-p
}

## * predictRisk.ranger
##' @export
##' @rdname predictRisk
##' @method predictRisk ranger
predictRisk.ranger <- function(object, newdata, times, cause, ...){
    if (object$treetype == "Regression") stop("Don't know how to predict risks based on regression trees.")
    xvars <- object$forest$independent.variable.names
    newdata <- subset(newdata,select=xvars)
    if (missing(times)||is.null(times)){
        if (object$treetype=="Classification") {
            stop("You must set the argument probability=TRUE in the call to ranger in order to use predictRisk.ranger.")
            pred <- stats::predict(object,data=newdata,predict.all=TRUE,...)$predictions
            p <- rowMeans(pred == 2)
        } else {
            p <- stats::predict(object,data=newdata,importance="none",...)$predictions[,2,drop=TRUE]
        }
        if (length(p) != NROW(newdata))
            stop(paste("\nPrediction vector has wrong length:\nRequested length of newdata: ",NROW(newdata)," \nProvided prediction vector: ",length(p),"\n\n",sep=""))
        p
    }else{
        if (object$treetype=="Survival") {
            ptemp <- 1-stats::predict(object,data=newdata,...)$survival
            pos <- prodlim::sindex(jump.times=ranger::timepoints(object),eval.times=times)
            if (NROW(newdata) == 1)
                p <- matrix(c(0,ptemp),nrow = 1)[,pos+1,drop=FALSE]
            else
                p <- cbind(0,ptemp)[,pos+1,drop=FALSE]
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
            p
        }else{
            message("Sorry. Hope Marvin does this soon.")
        }
    }
}

## * predictRisk.rfsrc
##' @export
##' @rdname predictRisk
##' @method predictRisk rfsrc
predictRisk.rfsrc <- function(object, newdata, times, cause, ...){
    if (missing(times)||is.null(times)){
        p <- stats::predict(object,newdata=newdata,importance="none",...)$predicted
        if (NCOL(p)>1)
        p <- as.numeric(p[,2,drop=TRUE])
        p
    }else{
        if (object$family[[1]]=="surv") {
            ptemp <- 1-stats::predict(object,newdata=newdata,importance="none",...)$survival
            pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
            p <- cbind(0,ptemp)[,pos+1,drop=FALSE]
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
            p
        }else{
            if (is.character(cause)) cause <- as.numeric(cause)
            if (!is.numeric(cause)) stop("cause is not numeric")
            cif <- stats::predict(object,newdata=newdata,importance="none",...)$cif[,,cause,drop=TRUE]
            # if necessary restore matrix format after dropping third dimension of array
            if (NROW(newdata)==1) {
                cif <- matrix(cif,nrow=1)
            }
            pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
            p <- cbind(0,cif)[,pos+1,drop=FALSE]
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
            p
        }
    }
}

## * predictRisk.FGR
##' @export
##' @rdname predictRisk
##' @method predictRisk FGR
predictRisk.FGR <- function(object,newdata,times,cause,...){
    ## require(cmprsk)
    ## predict.crr <- cmprsk:::predict.crr
    p <- predict(object=object,newdata=newdata,times=times)
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

## * predictRisk.riskRegression
##' @export
##' @rdname predictRisk
##' @method predictRisk riskRegression
predictRisk.riskRegression <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata,times=times)
  pos <- prodlim::sindex(jump.times=temp$time,eval.times=times)
  ## if (type=="survival")
      ## p <- cbind(1,1-temp$cuminc)[,pos+1,drop=FALSE]
  ## else
  p <- cbind(0,temp$risk)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}

## * predictRisk.ARR
##' @export
##' @rdname predictRisk
##' @method predictRisk ARR
predictRisk.ARR <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata,times=times)
  pos <- prodlim::sindex(jump.times=temp$time,eval.times=times)
  p <- cbind(0,temp$P1)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}

## * predictRisk.CauseSpecificCox
##' @export
##' @rdname predictRisk
##' @method predictRisk CauseSpecificCox
predictRisk.CauseSpecificCox <- function (object, newdata, times, cause,
                                          product.limit = TRUE, diag = FALSE, iid = FALSE, average.iid = FALSE, truncate = FALSE, ...) {
    dots <- list(...)
    type <- dots$type
    if(is.null(type)){
        type <- "absRisk"
    }

    ## prediciton
    outPred <- predict(object=object,
                       newdata=newdata,
                       times=times,
                       cause=cause,
                       keep.strata=FALSE,
                       se = FALSE,
                       iid = iid,
                       diag = diag,
                       average.iid = average.iid,
                       product.limit = product.limit,
                       store = dots$store, ## hidden argument
                       type = type)
    
    out <- outPred[[type]]
    if (truncate){
      out <- apply(out, c(1, 2), function(x) max(min(x,1),0))
    }
    if(iid){
        attr(out,"iid") <- outPred[[paste0(type,".iid")]]
    }
    if(average.iid){
        attr(out,"average.iid") <- outPred[[paste0(type,".average.iid")]]
    }
    return(out)
}


##' S3-wrapper for S4 function penalized
##'
##' S3-wrapper for S4 function penalized
##' @export
##' @param formula Communicated outcome and explanatory variables. See examples.
##' @param data Data set in which formula is to be interpreted
##' @param type String specifying the type of penalization. Should match one of the following values:
##' \code{"ridge"}, \code{"lasso"}, \code{"elastic.net"}.
##' @param lambda1 Lasso penalty
##' @param lambda2 ridge penalty
##' @param fold passed to \code{penalized::profL1}
##' @param ... Arguments passed to penalized
##' @examples
##' library(prodlim)
##' \dontrun{
##' ## too slow
##' if (require("penalized",quietly=TRUE)){
##' library(penalized)
##' set.seed(8)
##' d <- sampleData(200,outcome="binary")
##' newd <- sampleData(80,outcome="binary")
##' fitridge <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="ridge",
##'                 standardize=TRUE, model="logistic",trace=FALSE)
##' fitlasso <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="lasso",
##'                 standardize=TRUE, model="logistic",trace=FALSE)
##' # fitnet <- penalizedS3(Y~X1+X2+pen(7:8), data=d, type="elastic.net",
##' # standardize=TRUE, model="logistic",trace=FALSE)
##' predictRisk(fitridge,newdata=newd)
##' predictRisk(fitlasso,newdata=newd)
##' # predictRisk(fitnet,newdata=newd)
##' Score(list(fitridge),data=newd,formula=Y~1)
##' Score(list(fitridge),data=newd,formula=Y~1,split.method="bootcv",B=2)
##' data(nki70) ## S4 fit
##' fitS4 <- penalized(Surv(time, event), penalized = nki70[,8:77],
##'                  unpenalized = ~ER+Age+Diam+N+Grade, data = nki70,
##'                  lambda1 = 1)
##' fitS3 <- penalizedS3(Surv(time,event)~ER+Age+Diam+pen(8:77)+N+Grade,
##'                      data=nki70, lambda1=1)
##' ## or
##' penS3 <- penalizedS3(Surv(time,event)~ER+pen(TSPYL5,Contig63649_RC)+pen(10:77)+N+Grade,
##'                      data=nki70, lambda1=1)
##' ## also this works
##' penS3 <- penalizedS3(Surv(time,event)~ER+Age+pen(8:33)+Diam+pen(34:77)+N+Grade,
##'                     data=nki70, lambda1=1)
##' }}
##' @export
penalizedS3 <- function(formula,
                        data,
                        type="elastic.net",
                        lambda1,
                        lambda2,
                        fold,
                        ...){
    responseFormula <- stats::update(formula,~1)
    responsevars <- all.vars(responseFormula)
    if (length(responsevars)==1){
        FRAME  <- Publish::specialFrame(formula,
                                        data,
                                        specials=c("pen","unpen"),
                                        strip.specials=c("pen","unpen"),
                                        strip.unspecials="pen",
                                        strip.arguments=NULL,
                                        ## strip.alias=list("pen"="penalized","unpen"="unpenalized"),
                                        drop.intercept=TRUE,
                                        specials.design=TRUE,
                                        response=TRUE,
                                        na.action=options()$na.action)
        response <- FRAME$response[[1]]
        pen <- FRAME$pen
        unpen <- FRAME$unpen
        if (is.null(unpen))
            args <- list(response=response,penalized=pen,data=data,...)
        else
            args <- list(response=response,penalized=pen,unpenalized=unpen,data=data,...)
        #modelframe <- stats::model.frame(formula=formula,data=data,na.action=na.fail)
    }else{
        # {{{ distangle the formula
        EHF <- prodlim::EventHistory.frame(formula=formula,
                                           data=data,
                                           specials=c("pen","unpen"),
                                           stripSpecials=c("pen","unpen"),
                                           stripUnspecials="pen",
                                           # stripAlias=list("penalized"="pen","unpenalized"="unpen"),
                                           specialsDesign=TRUE)
        response <- EHF$event.history
        pen <- EHF$pen
        unpen <- EHF$unpen
        if (is.null(unpen))
            args <- list(response=response,penalized=pen,data=data,...)
        else
            args <- list(response=response,penalized=pen,unpenalized=unpen,data=data,...)
    }
    # {{{ find optimal L1
    if ((tolower(type) %in% c("lasso","elastic.net")) && missing(lambda1)){
        if (missing(fold)) fold <- 1:NROW(data)
        las1 <- do.call(penalized::profL1,c(args,list(fold=fold)))
        lambda1 <- do.call(penalized::optL1,c(args,list(fold=las1$fold)))$lambda
    }
    # }}}
    # {{{ find optimal L2
    if ((tolower(type) %in% c("ridge","elastic.net")) && missing(lambda2)){
        if (missing(fold)) fold <- 1:NROW(data)
        lambda2 <- do.call(penalized::optL2,
                           c(args,list(fold=fold)))$lambda
    }
    # }}}
    # {{{ call S4 method
    ## unpenalized terms are communicated via
    ## the left hand side of response
    fitS4 <- switch(type,"ridge"={
        do.call(penalized::penalized,c(args,list(lambda1=0, lambda2=lambda2)))
    },
    "lasso"={
        do.call(penalized::penalized,c(args,list(lambda1=lambda1, lambda2=0)))
    },
    "elastic.net"={
        do.call(penalized::penalized,c(args,list(lambda1=lambda1, lambda2=lambda2)))
    })
    if (is.infinite(fitS4@loglik)){
        print(fitS4)
        stop()
    }
    # }}}
    # {{{ deliver S3 object
    fit <- list(fitS4=fitS4,call=match.call())
    fit$terms <- terms(formula)
    class(fit) <- "penfitS3"
    fit
    # }}}
}

## * predictRisk.penfitS3
##' @export
##' @rdname predictRisk
##' @method predictRisk penfitS3
predictRisk.penfitS3 <- function(object,
                                 newdata,
                                 times,
                                 ...){
    penfit <- object$fitS4
    if (missing(newdata)) stop("Argument 'newdata' is missing")
    if (NROW(newdata) == 0) stop("No (non-missing) observations")
    rhs <- as.formula(delete.response(object$terms))
    responseFormula <- stats::update(object$terms,~1)
    responsevars <- all.vars(responseFormula)
    if (length(responsevars)==1){
        dummy.formula=stats::update.formula(rhs,paste0(responsevars,"~."))
        FRAME  <- Publish::specialFrame(formula=responseFormula,
                                        data=newdata,
                                        specials=c("pen","unpen"),
                                        strip.specials=c("pen","unpen"),
                                        strip.unspecials="pen",
                                        strip.arguments=NULL,
                                        ## strip.alias=list("pen"="penalized","unpen"="unpenalized"),
                                        drop.intercept=TRUE,
                                        specials.design=TRUE,
                                        response=TRUE,
                                        na.action=options()$na.action)
        response <- FRAME$response[[1]]
        pen <- FRAME$pen
        unpen <- FRAME$unpen
        if (is.null(unpen))
            args <- list(penfit,penalized=pen,data=newdata,...)
        else
            args <- list(penfit,penalized=pen,unpenalized=unpen,data=newdata,...)
        #modelframe <- stats::model.frame(formula=formula,data=data,na.action=na.fail)
    }else{
        newdata$dummy.time=rep(1,NROW(newdata))
        newdata$dummy.event=rep(1,NROW(newdata))
        dummy.formula=stats::update.formula(rhs,"Hist(dummy.time,dummy.event)~.")
        EHF <- prodlim::EventHistory.frame(formula=dummy.formula,
                                           data=newdata,
                                           specials=c("pen","unpen"),
                                           stripSpecials=c("pen","unpen"),
                                           stripUnspecials="pen",
                                           specialsDesign=TRUE)
        pen <- EHF$pen
        unpen <- EHF$unpen
        args <- list(penfit,penalized=pen)
    }
    if (length(unpen)>0) args <- c(args,list(unpenalized=unpen))
    p <- do.call(penalized::predict,args)
    if (penfit@model=="cox"){
        if (missing(times)) stop("Need time points for predicting risks from a penalized Cox regression model.")
        ttt <- p@time
        p <- cbind(0,1-p@curves)[,prodlim::sindex(jump.times=ttt,eval.times=times)+1,drop=FALSE]
        if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
        }
    }
    p
}

##' @title SmcFcs
##' @description TODO
##'
##' @param formula TODO
##' @param data TODO
##' @param m TODO
##' @param method TODO
##' @param fitter TODO
##' @param fit.formula TODO
##' @param ... TODO
##'
##' # @export
SmcFcs  <- function(formula,data,m=5,method,fitter="glm",fit.formula,...){
    requireNamespace("smcfcs",quietly=FALSE)
    this <- as.character(formula)
    Yname <- this[2]
    Xnames <- this[3]
    sform <- paste0(Yname,"~",Xnames)
    if (length(unique(data[[Yname]]))!=2)
        stop("Outcome must be binary")
    Xvars <- all.vars(formula)
    if (missing(fit.formula)) fit.formula <- formula
    Avars <- all.vars(fit.formula)
    Xvars <- union(Xvars,Avars)
    Xdata <- subset(data,select=Xvars)
    if (missing(method)){
        method <- sapply(Xdata,function(x){
            if (any(is.na(x))){
                if(length(unique(x))==2){
                    "logreg"
                }else{
                    if(is.factor(x)){
                        "mlogit"
                    } else{
                        "norm"
                    }
                }
            } else{
                ""
            }
        })
    }
    idata.list <- smcfcs::smcfcs(smformula=sform,
                                 originaldata=Xdata,
                                 m=m,
                                 smtype="logistic",...,
                                 method=method)$impDatasets
    res <- lapply(idata.list,function(d){
        do.call(fitter,list(fit.formula,data=d,family="binomial"))
    })
    class(res) <- "SmcFcs"
    res
}

##' @export
predictRisk.SmcFcs <- function(object,newdata,...){
    p <- Reduce("+",lapply(object,function(x){
        predictRisk(x,newdata=newdata)
    }))/length(object)
    p
}

## * predictRisk.SuperPredictor
##' @export
##' @rdname predictRisk
##' @method predictRisk SuperPredictor
predictRisk.SuperPredictor  <- function(object,newdata,...){
    p <- SuperLearner::predict.SuperLearner(object=object,newdata=newdata)$pred
    ## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        ## stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

## * predictRisk.gbm
##' @export
##' @rdname predictRisk
##' @method predictRisk gbm
predictRisk.gbm <- function(object, newdata, times, ...) {
    requireNamespace("gbm")
    n.trees <- object$n.trees
    traindata <-  gbm::reconstructGBMdata(object)
    p <- matrix(0, NROW(newdata), length(times))
    xb.train <- predict(object ,newdata = traindata, n.trees = n.trees)
    H2 <- gbm::basehaz.gbm(t = traindata[, as.character(object$call$formula[[2]][[2]])],
                           delta = traindata[, as.character(object$call$formula[[2]][[3]])],
                           f.x = xb.train, t.eval = times)
    xb.test <- predict(object, newdata = newdata , n.trees = n.trees )
    for (i in 1:length(times)) p[,i] <- exp(-H2[i] * exp(xb.test))
    p[,times==0] <- 1
    1 - p
}
## * predictRisk.flexsurvreg
##' @export
##' @rdname predictRisk
##' @method predictRisk flexsurvreg
predictRisk.flexsurvreg <- function(object, newdata, times, ...) {
    requireNamespace("flexsurv")
    newdata <- data.frame(newdata)
    p <- matrix(0, NROW(newdata), length(times))
    term <- attr(terms(as.formula(object$call$formula)), "term.labels")
    sm <- summary(object, newdata = newdata[, term, drop = FALSE], t = times, start = 0, B = 0) #no confidence interval simulations
    for (i in 1:NROW(newdata)){
        p[i,] <- sm[[i]][,2]
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop("Prediction failed")
    1 - p
}

##' @export
##' @rdname predictRisk
##' @method predictRisk Hal9001
predictRisk.Hal9001 <- function(object,
                                newdata,
                                times,
                                cause,
                                ...){
    stopifnot(object$family=="cox")
    newdata$dummy.time=rep(1,NROW(newdata))
    newdata$dummy.event=rep(1,NROW(newdata))
    rhs <- as.formula(delete.response(object$terms))
    dummy.formula=stats::update.formula(rhs,"Hist(dummy.time,dummy.event)~.")
    EHF <- prodlim::EventHistory.frame(formula=dummy.formula,
                                       data=newdata,
                                       specials = NULL,
                                       unspecialsDesign=TRUE)
    newdata$dummy.time = NULL
    newdata$dummy.event = NULL
    # blank Cox object obtained with riskRegression:::coxModelFrame
    info <- object$surv_info
    hal_pred <- as.numeric(predict(object$fit,new_data=EHF$design))
    L0 <- riskRegression::baseHaz_cpp(starttimes = info$start,
                                      stoptimes = info$stop,
                                      status = info$status,
                                      eXb = info$eXb,
                                      strata = 1,
                                      nPatients = NROW(info$stop),
                                      nStrata = 1,
                                      emaxtimes = max(info$stop),
                                      predtimes = times, 
                                      cause = 1,
                                      Efron = TRUE,
                                      reverse = FALSE)
    hal_Surv <- exp(-hal_pred%o%L0$cumhazard)
    p <- 1-hal_Surv
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    }
    p
}


## * predictRisk.GLMnet
##' @rdname predictRisk
##' @method predictRisk GLMnet
##' @export
predictRisk.GLMnet <- function(object,
                               newdata,
                               times,
                               product.limit = FALSE,
                               diag = FALSE,
                               ...){
    has_survival <- inherits(object$fit,"coxnet")
    dots <- list(...)
    type <- dots$type ## hidden argument for ate
    if (has_survival){
        pred <- 1-predictCox(object=object,
                             newdata=newdata,
                             times=times,
                             iid = FALSE,
                             confint = FALSE,
                             diag = diag,
                             average.iid = FALSE,
                             product.limit = product.limit,
                             type="survival")$survival
        if(identical(type,"survival")){
            pred <- 1-pred
        }
        pred
    }
}

## * predictRisk.GLMnet2
##' @rdname predictRisk
##' @method predictRisk GLMnet2
##' @export
predictRisk.GLMnet2 <- function(object,newdata,times=NA,...) {
    args <- list(...)
    # check if user has supplied a lambda value
    slambda <- args$lambda
    # check if object has saved a selected lambda value
    if (length(slambda) == 0)
        slambda <- attr(object,"selected.lambda")
    if (length(slambda) == 0 || length(slambda)>1 || !(is.numeric(slambda))){
        stop("You must choose a single numeric lambda value for predictRisk ... ")
    }
    pos.lambda <- match(slambda,object$fit$lambda,nomatch = 0)
    if (pos.lambda == 0){
        stop("The fitted model was not fitted with the specified penalty parameter (lambda)")
    }
    lambda=cv=NULL
    # library(glmnet)
    # requireNamespace(c("prodlim","glmnet"))
    # predict.cv.glmnet <- utils::getFromNamespace("predict.cv.glmnet","glmnet")
    # predict.glmnet <- utils::getFromNamespace("predict.glmnet","glmnet")
    rhs <- as.formula(delete.response(object$terms))
    if (length(info <- object$surv_info) == 0){
        xnew <- model.matrix(rhs,data=newdata)
        if (is.null(slambda) && object$cv){
            p <- predict(object$fit,newx=xnew,type = "response", s="lambda.min")
        }
        else if (pos.lambda == 0 && !object$cv){
            if (length(object$lambda) == 1){
                p <- predict(object$fit,newx=xnew,type = "response", s=object$lambda)
            }
            else {
                stop("Object fitted with multiple lambdas. You must pick one lambda for predictRisk!")
            }
        }
        else {
            p <- predict(object$fit,newx=xnew,type = "response", s=slambda)
        }
    } else {
        # convert covariates to dummy variables
        newdata$dummy.time=rep(1,NROW(newdata))
        newdata$dummy.event=rep(1,NROW(newdata))
        dummy.formula=stats::update.formula(rhs,"prodlim::Hist(dummy.time,dummy.event)~.")
        EHF <- prodlim::EventHistory.frame(formula=dummy.formula,data=newdata,specials = NULL,unspecialsDesign=TRUE)
        newdata$dummy.time = NULL
        newdata$dummy.event = NULL
        # blank Cox object obtained with riskRegression:::coxModelFrame
        if (pos.lambda == 0 && object$cv){
            GLMnet_pred <- c(exp(predict(object$fit,newx=EHF$design,type = "link", s="lambda.min")))
            lambda <- object$fit$lambda.min ## is needed for train_eXb
        }
        else if (pos.lambda == 0 && !object$cv){
            if (length(object$lambda) == 1){
                GLMnet_pred <- c(exp(predict(object$fit,newx=EHF$design,type = "link", s=object$lambda)))
                lambda <- object$lambda
            }
            else {
                stop("Object fitted with multiple lambdas. You must pick a single value for lambda.")
            }
        } else {
            if (all((pos.lambda)>0)){
                GLMnet_pred <- c(exp(predict(object$fit,newx=EHF$design,type = "link", s=slambda)))
                lambda <- slambda
            }
            else {
                stop(paste0("The fitted model was not fitted with the following penalty parameters (lambdas): ",
                            paste0(slambda[pos.lambda == 0],collapse = ", ")))
            }
        }
        train_eXb <- c(exp(predict(object$fit,newx=object$sorted_x_train,type = "link", s=lambda)))
        L0 <- riskRegression::baseHaz_cpp(starttimes = info$start,
                                          stoptimes = info$stop,
                                          status = info$status,
                                          eXb = train_eXb,
                                          strata = 1,
                                          nPatients = NROW(info$stop),
                                          nStrata = 1,
                                          emaxtimes = max(info$stop),
                                          predtimes = sort(unique(info$stop)),
                                          cause = 1,
                                          Efron = TRUE,
                                          reverse = FALSE)$cumhazard
        GLMnetSurv <- exp(-GLMnet_pred%o%L0)
        where <- sindex(jump.times=unique(info$stop),eval.times=times)
        p <- cbind(0,1-GLMnetSurv)[,1+where]
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    }
    p
}


## * predictRisk.singleEventCB
##' @export
##' @rdname predictRisk
##' @method predictRisk singleEventCB
predictRisk.singleEventCB <- function(object, newdata, times, cause, ...) {
    requireNamespace("casebase")
    if (object$call[[1]]!="glm"){
        #get all covariates excluding intercept and time
        coVars=colnames(object$originalData$x)
        #coVars is used in lines 44 and 50
        newdata=data.matrix(drop(subset(newdata, select=coVars)))
    }
    # if (missing(cause)) stop("Argument cause should be the event type for which we predict the absolute risk.")
    # the output of absoluteRisk is an array with dimension dependening on the length of the requested times:
    # case 1: the number of time points is 1
    #         dim(array) =  (length(time), NROW(newdata), number of causes in the data)
    if (length(times) == 1) {
        a <- casebase::absoluteRisk(object, newdata = newdata, time = times)
        p <- matrix(a, ncol = 1)
    } else {
        # case 2 a) zero is included in the number of time points
        if (0 %in% times) {
            # dim(array) =  (length(time)+1, NROW(newdata)+1, number of causes in the data)
            a <- casebase::absoluteRisk(object, newdata = newdata, time = times)
            p <- t(a)
        } else {
            # case 2 b) zero is not included in the number of time points (but the absoluteRisk function adds it)
            a <- casebase::absoluteRisk(object, newdata = newdata, time = times)
            ### we need to invert the plot because, by default, we get cumulative incidence
            #a[, -c(1)] <- 1 - a[, -c(1)]
            ### we remove time 0 for everyone, and remove the time column
            a <- a[-c(1), -c(1)] ### a[-c(1), ] to keep times column, but remove time 0 probabilities
            # now we transpose the matrix because in riskRegression we work with number of
            # observations in rows and time points in columns
            p <- t(a)
        }
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    }
    p
}


## GrpSurv <- function(formula,data,...){
    ## requireNamespace(c("grpreg","prodlim"))
    ## EHF = prodlim::EventHistory.frame(formula,data,unspecialsDesign = TRUE,specials = NULL)
    ## fit = grpreg::grpsurv(X = EHF$design,y = EHF$event.history,...)
    ## fit = list(fit = fit,terms = terms(formula),call=match.call())
    ## class(fit) = c("GrpSurv",class(fit))
    ## fit
## }

## predictRisk.GrpSurv <- function(object, newdata, times, cause, ...){
    ## newdata$dummy.time=rep(1,NROW(newdata))
    ## newdata$dummy.event=rep(1,NROW(newdata))
    ## rhs <- as.formula(delete.response(object$terms))
    ## dummy.formula=stats::update.formula(rhs,"Hist(dummy.time,dummy.event)~.")
    ## EHF <- prodlim::EventHistory.frame(formula=dummy.formula,
                                       ## data=newdata,
                                       ## specials = NULL,
                                       ## unspecialsDesign=TRUE)
    ## newdata$dummy.time = NULL
    ## newdata$dummy.event = NULL
    ## p <- predict(object$fit, EHF$design, type="survival")
    ## p <- 1-sapply(p,function(f)f(times))
    ## if (length(times) == 1){
        ## p = cbind(p)
    ## }else(p = t(p))
    ## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
        ## stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    ## }
    ## p
## }

## XgbSurv <- function(formula,data,...){
    ## requireNamespace(c("survXgboost","xgboost","prodlim"))
    ## EHF = prodlim::EventHistory.frame(formula,data,unspecialsDesign = TRUE,specials = NULL)
    ## fit = survXgboost::xgb.train.surv(data = EHF$design,label = ifelse(EHF$event.history[,2] == 1, EHF$event.history[,1], -EHF$event.history[,1]),...)
    ## fit = list(fit = fit,terms = terms(formula),call=match.call())
    ## class(fit) = c("XgbSurv",class(fit))
    ## fit
## }

## predictRisk.XgbSurv <- function(object, newdata, times, cause, ...){
    ## newdata$dummy.time=rep(1,NROW(newdata))
    ## newdata$dummy.event=rep(1,NROW(newdata))
    ## rhs <- as.formula(delete.response(object$terms))
    ## dummy.formula=stats::update.formula(rhs,"Hist(dummy.time,dummy.event)~.")
    ## EHF <- prodlim::EventHistory.frame(formula=dummy.formula,
                                       ## data=newdata,
                                       ## specials = NULL,
                                       ## unspecialsDesign=TRUE)
    ## newdata$dummy.time = NULL
    ## newdata$dummy.event = NULL
    ## p <-predict(object = object$fit, newdata = EHF$design, type = "surv", times = times)
    ## p <- 1-p
    ## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
        ## stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    ## }
    ## p
## }

## * predictRisk.wglm
#' @rdname predictRisk
#' @method predictRisk wglm
#' @export
predictRisk.wglm <- function(object, newdata, times = NULL, 
                             product.limit = NULL, diag = FALSE, iid = FALSE, average.iid = FALSE, ...){

    dots <- list(...)
    if(is.null(dots$se)){ ## hidden se argument
        se <- "robust"
    }else{
        se <- match.arg(se, c("robust","robust-wknown"))
    }

    ## ** extract information and normalize arguments
    if(is.null(times)){
        times <- object$times
    }else{
        if(any(times %in% object$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }

    if(is.null(product.limit)){
        product.limit <- object$product.limit
    }
    n.times <- length(times)
    n.sample <- sum(coxN(object))
    n.newdata <- NROW(newdata)

    if(average.iid){
        if(is.null(attr(average.iid,"factor"))){
            factor <- list(matrix(1, nrow = n.newdata, ncol = n.times))
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
            for(iFactor in 1:length(factor)){ ## iFactor <- 1
                ## when only one column and diag = FALSE, use the same weights at all times
                if((diag == FALSE) && (NCOL(factor[[iFactor]])==1) && (n.times > 1)){
                    factor[[iFactor]] <- matrix(factor[[iFactor]][,1],
                                                nrow = NROW(factor[[iFactor]]),
                                                ncol = n.times, byrow = FALSE)
                }
                ## check dimensions
                if(any(dim(factor[[iFactor]])!=c(n.newdata, diag + (1-diag)*n.times))){
                    stop("Attribute \"factor\" of argument \'average.iid\' must be a list of matrices of size ",n.newdata,",",diag + (1-diag)*n.times," \n")
                }
            }
        }
        n.factor <- length(factor)
    }

    ## hidden argument: enable to ask for the prediction of Y==1 or Y==0
    level <- list(...)$level
    type <- dots$type ## hidden argument for ate
    if(identical(type,"survival")){
        stop("Unkown argument \'type\' for predictRisk.wglm: use argument \'level\' instead. \n")
    }

    ## ** prepare output
    out <- matrix(NA, nrow = n.newdata, ncol = n.times)
    if(iid){
        attr(out,"iid") <- array(NA, dim = c(n.sample, n.times, n.newdata))
    }
    if(average.iid){
        attr(out,"average.iid") <- lapply(1:n.factor, function(x){matrix(NA, nrow = n.sample, ncol = n.times)})
        if(!is.null(names(factor))){
            names(attr(out,"average.iid")) <- names(factor)
        }
    }

    if(iid || average.iid){
        if(se=="robust"){
            object.iid <- lava::iid(object, simplify = FALSE, times = times)
        }else if(se == "robust-wknown"){
            object.iid <- lapply(object$fit[match(times, object$times)],lava::iid)
        }
    }

    
    for(iTime in 1:n.times){
        iTime2 <- which(object$times == times[iTime])
        iObject <- object$fit[[iTime2]]
        iFormula <- stats::formula(iObject)
        
        ## ** identify correct level
        if(!is.null(level)){
            matching.Ylevel <- table(iObject$data[[all.vars(formula(iObject))[1]]],
                                     iObject$y)
            all.levels <- rownames(matching.Ylevel)
            level <- match.arg(level, all.levels)

            index.level <- which(matching.Ylevel[level,]>0)
            if(length(index.level) > 1){
                stop("Unknown value for the outcome variable \n")
            }
        }else{
            index.level <- 2
        }

        ## ** point estimate
        if(index.level == 1){
            out[,iTime] <- 1-predict(iObject, type = "response", newdata = newdata, se = FALSE)
        }else{
            out[,iTime] <- predict(iObject, type = "response", newdata = newdata, se = FALSE)
        }
        
        ## ** uncertainty (chain rule)
        if(iid || average.iid){
            iid.beta <- object.iid[[iTime]]
            newX <- model.matrix(delete.response(terms(iFormula)), newdata)
            Xbeta <- predict(iObject, type = "link", newdata = newdata, se = FALSE)

            if(average.iid){
                for(iFactor in 1:n.factor){
                    iE.X <- colMeans(colMultiply_cpp(newX, scale = factor[[iFactor]][,iTime] * exp(-Xbeta)/(1+exp(-Xbeta))^2))
                    if(index.level == 1){
                        attr(out,"average.iid")[[iFactor]][,iTime] <- -iid.beta %*% iE.X
                    }else{
                        attr(out,"average.iid")[[iFactor]][,iTime] <- iid.beta %*% iE.X
                    }
                }
            }
            if(iid){
                if(index.level == 1){
                    attr(out,"iid")[,iTime,] <- -iid.beta %*% t(colMultiply_cpp(newX, scale = exp(-Xbeta)/(1+exp(-Xbeta))^2))
                }else{
                    attr(out,"iid")[,iTime,] <- iid.beta %*% t(colMultiply_cpp(newX, scale = exp(-Xbeta)/(1+exp(-Xbeta))^2))
                }
            }
        }
    }

    ## ** export
    if(average.iid && is.null(attr(average.iid,"factor"))){
        attr(out,"average.iid") <- attr(out,"average.iid")[[1]]
    }
    return(out)
}

#----------------------------------------------------------------------
### predictRisk.R ends here

