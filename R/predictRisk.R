### predictRisk.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (09:02) 
## Version: 
## last-updated: aug 31 2020 (17:14) 
##           By: Brice Ozenne
##     Update #: 332
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
#' @param landmark The starting time for the computation of the cumulative risk.
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
#'\dontrun{
#' library(SuperLearner)
#' set.seed(1)
#' sl = SuperLearner(Y = d$Y, X = d[,-1], family = binomial(),
#'       SL.library = c("SL.mean", "SL.glmnet", "SL.randomForest"))
#'}
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
#' # stratified cox model
#' sfit <- coxph(Surv(time,event)~strata(X1)+X2,data=d,x=TRUE,y=TRUE)
#' predictRisk(sfit,newdata=d[1:3,],times=c(1,3,5,10))
#' 
#' ## simulate learning and validation data
#' learndat <- sampleData(100,outcome="survival")
#' valdat <- sampleData(100,outcome="survival")
#' ## use the learning data to fit a Cox model
#' library(survival)
#' fitCox <- coxph(Surv(time,event)~X1+X2,data=learndat,x=TRUE,y=TRUE)
#' ## suppose we want to predict the survival probabilities for all subjects
#' ## in the validation data at the following time points:
#' ## 0, 12, 24, 36, 48, 60
#' psurv <- predictRisk(fitCox,newdata=valdat,times=seq(0,60,12))
#' ## This is a matrix with event probabilities (1-survival)
#' ## one column for each of the 5 time points
#' ## one row for each validation set individual
#' 
#' # Do the same for a randomSurvivalForest model
#' # library(randomForestSRC)
#' # rsfmodel <- rfsrc(Surv(time,event)~X1+X2,data=learndat)
#' # prsfsurv=predictRisk(rsfmodel,newdata=valdat,times=seq(0,60,12))
#' # plot(psurv,prsfsurv)
#' 
#' ## Cox with ridge option
#' f1 <- coxph(Surv(time,event)~X1+X2,data=learndat,x=TRUE,y=TRUE)
#' f2 <- coxph(Surv(time,event)~ridge(X1)+ridge(X2),data=learndat,x=TRUE,y=TRUE)
#' \dontrun{
#' plot(predictRisk(f1,newdata=valdat,times=10),
#'      riskRegression:::predictRisk.coxph(f2,newdata=valdat,times=10),
#'      xlim=c(0,1),
#'      ylim=c(0,1),
#'      xlab="Unpenalized predicted survival chance at 10",
#'      ylab="Ridge predicted survival chance at 10")
#'}
#' 
#' ## competing risks
#' 
#' library(survival)
#' library(riskRegression)
#' library(prodlim)
#' train <- prodlim::SimCompRisk(100)
#' test <- prodlim::SimCompRisk(10)
#' cox.fit  <- CSC(Hist(time,cause)~X1+X2,data=train)
#' predictRisk(cox.fit,newdata=test,times=seq(1:10),cause=1)
#'
#' ## with strata
#' cox.fit2  <- CSC(list(Hist(time,cause)~strata(X1)+X2,Hist(time,cause)~X1+X2),data=train)
#' predictRisk(cox.fit2,newdata=test,times=seq(1:10),cause=1)
#'
#' @export 
predictRisk <- function(object,newdata,...){
    UseMethod("predictRisk",object)
}

## * predictRisk.default
##' @export
#' @rdname predictRisk
#' @method predictRisk default
predictRisk.default <- function(object,newdata,times,cause,...){
    stop(paste0("No method available for evaluating predicted probabilities from objects in class: ",class(object),". But, you can write it yourself or ask the package manager."),call.=FALSE)
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

            ## ** compute influence function of the coefficients using lava
            iid.beta <- lava::iid(object)

            newX <- model.matrix(stats::formula(object), newdata)
            Xbeta <- predict(object, type = "link", newdata = newdata, se = FALSE)
            
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
        level <- list(...)$level
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
        return(out)
    } else {
        stop("Currently only the binomial family is implemented for predicting a status from a glm object.")
    }
}
## * predictRisk.multinom
predictRisk.multinom <- function(object, newdata, iid = FALSE, average.iid = FALSE,...){
    n.obs <- NROW(newdata)
    n.class <- length(object$lev)
    n.coef <- length(coef(object))
    n.coefperY <- NROW(coef(object))
    out <- predict(object, newdata = newdata, type = "probs")

    if(n.class == 2){
        out <- cbind(1-out, out)
        colnames(out) <- object$lev
    }

    ## ** set correct level
    ## hidden argument: enable to ask for the prediction of a specific level
    level <- list(...)$level
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
                apply(iFactor, 2, function(iiFactor){
                    iRes <- Reduce("+",lapply(object$lev[-1], function(iClass){ ## iClass <- "1"
                        - iid.beta.level[[iClass]] %*% colMeans(colMultiply_cpp(newX, scale = iiFactor * exp(Xbeta[,iClass]+Xbeta[,level])/eXbeta_rowSum^2))
                    }))
                    if(which(level %in% object$lev)>1){
                        iRes <- iRes + iid.beta.level[[level]] %*% colMeans(colMultiply_cpp(newX, scale = iiFactor * exp(Xbeta[,level])/eXbeta_rowSum))
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
    ## require(timereg)
    stop("FIXME")
    time.coef <- data.frame(object$cum)
    ntime <- nrow(time.coef)
    objecttime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    covanames <- names(time.coef)[-(1:2)]
    notfound <- match(covanames,names(newdata),nomatch=0)==0
    if (any(notfound))
        stop("\nThe following predictor variables:\n\n",
             paste(covanames[notfound],collapse=","),
             "\n\nwere not found in newdata, which only provides the following variables:\n\n",
             paste(names(newdata),collapse=","),
             "\n\n")
    time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
    nobs <- nrow(newdata)
    ## hazard <- .C("survest_cox_aalen",
                 ## timehazard=double(ntime*nobs),
                 ## as.double(unlist(time.coef[,-1])),
                 ## as.double(unlist(time.vars)),
                 ## as.integer(ntimevars+1),
                 ## as.integer(nobs),
                 ## as.integer(ntime),PACKAGE="pec")$timehazard
    hazard <- matrix(hazard,ncol=ntime,nrow=nobs,dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
    surv <- pmin(exp(-hazard),1)
    if (missing(times)) times <- sort(unique(objecttime))
    p <- surv[,prodlim::sindex(jump.times=objecttime,eval.times=times)]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    1-p
}

## * predictRisk.cox.aalen
##' @export
##' @rdname predictRisk
##' @method predictRisk cox.aalen
predictRisk.cox.aalen <- function(object,newdata,times,...){
    #  require(timereg)
    ##  The time-constant effects first
    stop("FIXME")
    const <- c(object$gamma)
    names(const) <- substr(dimnames(object$gamma)[[1]],6,nchar(dimnames(object$gamma)[[1]])-1)
    constant.part <- t(newdata[,names(const)])*const
    constant.part <- exp(colSums(constant.part))
    ##  Then extract the time-varying effects
    time.coef <- data.frame(object$cum)
    ntime <- nrow(time.coef)
    objecttime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
    nobs <- nrow(newdata)
    ## time.part <- .C("survest_cox_aalen",timehazard=double(ntime*nobs),as.double(unlist(time.coef[,-1])),as.double(unlist(time.vars)),as.integer(ntimevars+1),as.integer(nobs),as.integer(ntime),PACKAGE="pec")$timehazard
    time.part <- matrix(time.part,ncol=ntime,nrow=nobs)
    ## dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
    surv <- pmin(exp(-time.part*constant.part),1)
    if (missing(times)) times <- sort(unique(objecttime))
    p <- surv[,prodlim::sindex(objecttime,times)]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    1-p
}

    
## * predictRisk.coxph
##' @export
##' @rdname predictRisk
##' @method predictRisk coxph
predictRisk.coxph <- function(object, newdata, times,
                              product.limit = FALSE, diag = FALSE, iid = FALSE, average.iid = FALSE, ...){
    dots <- list(...)
    type <- dots$type ## hidden argument for ate
    store.iid <- dots$store ## hidden argument for ate
    
    if(product.limit){
        outPred <- predictCoxPL(object=object,
                                newdata=newdata,
                                times=times,
                                iid = iid,
                                average.iid = average.iid,
                                keep.times=FALSE,
                                diag = diag,
                                store.iid = store.iid,
                                type="survival")
    }else{
        outPred <- predictCox(object=object,
                              newdata=newdata,
                              times=times,
                              iid = iid,
                              diag = diag,
                              average.iid = average.iid,
                              store.iid = store.iid,
                              type="survival")
    }
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

## * predictRisk.coxphTD
##' @export
##' @rdname predictRisk
##' @method predictRisk coxphTD
predictRisk.coxphTD <- function(object,newdata,times,landmark,...){
    stopifnot(attr(object$y,"type")=="counting")
    bh <- survival::basehaz(object,centered=TRUE)
    Lambda0 <- bh[,1]
    etimes <- bh[,2]
    lp <- predict(object,newdata=newdata,type="lp")
    p <- do.call("cbind",lapply(times,function(ttt){
        index <- prodlim::sindex(eval.times=c(landmark,landmark+ttt),jump.times=etimes)
        Lambda0.diff <- c(0,Lambda0)[1+index[2]] - c(0,Lambda0)[1+index[1]]
        1-exp(-Lambda0.diff * exp(lp))
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    return(p)
}


## * predictRisk.CSCTD
##' @export
##' @rdname predictRisk
##' @method predictRisk CSCTD
predictRisk.CSCTD <- function(object,newdata,times,cause,landmark,...){
    stopifnot(attr(object$models[[1]]$y,"type")=="counting")
    if (missing(cause)) cause <- object$theCause
    else{
        ## cause <- prodlim::checkCauses(cause,object$response)
        cause <- unique(cause)
        if (!is.character(cause)) cause <- as.character(cause)
        fitted.causes <- prodlim::getStates(object$response)
        if (!(all(cause %in% fitted.causes))){
            stop(paste0("Cannot find requested cause(s) in object\n\n",
                        "Requested cause(s): ",
                        paste0(cause,collapse=", "),
                        "\n Available causes: ",
                        paste(fitted.causes,collapse=", "),"\n"))
        }
    }
    causes <- object$causes
    index.cause <- which(causes == cause)
    bh <- lapply(1:length(object$models),function(m){
        survival::basehaz(object$models[[m]],centered=TRUE)
    })
    lp <- lapply(1:length(object$models),function(m){
        predict(object$models[[m]],
                newdata=newdata,
                type="lp")
    })
    p <- do.call("cbind",lapply(times,function(ttt){
        hazard <- lapply(1:length(object$models),function(m){
            index <- sindex(eval.times=c(landmark,landmark+ttt),jump.times=bh[[m]][,2])
            bh.m <- bh[[m]][,1]
            exp(lp[[m]])*(c(0,bh.m)[1+index[2]] - c(0,bh.m)[1+index[1]])
        })
        surv <- exp(- Reduce("+",hazard))
        surv*hazard[[index.cause]]
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    return(p)
}



## * predictRisk.coxph.penal
##' @export
##' @rdname predictRisk
##' @method predictRisk coxph.penal
predictRisk.coxph.penal <- function(object,newdata,times,...){
  frailhistory <- object$history$'frailty(cluster)'$history
  frailVar <- frailhistory[NROW(frailhistory),1]
  linearPred <- predict(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
  basehaz <- basehaz(object)
  bhTimes <- basehaz[,2]
  bhValues <- basehaz[,1]
  survPred <- do.call("rbind",lapply(1:NROW(newdata),function(i){
    (1+frailVar*bhValues*exp(linearPred[i]))^{-1/frailVar}
  }))
  where <- prodlim::sindex(jump.times=bhTimes,eval.times=times)
  p <- cbind(1,survPred)[,where+1]
  if ((miss.time <- (length(times) - NCOL(p)))>0)
    p <- cbind(p,matrix(rep(NA,miss.time*NROW(p)),nrow=NROW(p)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  1-p
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
predictRisk.prodlim <- function(object,newdata,times,cause,...){
    ## require(prodlim)
    if (object$model[[1]]=="competing.risks" && missing(cause))
        stop(paste0("Cause is missing. Should be one of the following values: ",paste(attr(object$model.response,"states"),collapse=", ")))
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
    p
}

## * predictRisk.survfit
##' @export
##' @rdname predictRisk
##' @method predictRisk survfit
predictRisk.survfit <- function(object,newdata,times,...){
    p <- predict.survfit(object,newdata=newdata,times=times,type="cuminc",bytimes=TRUE,fill="last")
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

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
    sfit <- summary(object,times=times)
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
        stratdat <- newdata[,covars,drop=FALSE]
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
    xvars <- object$forest$independent.variable.names
    newdata <- subset(newdata,select=xvars)
    if (missing(times)||is.null(times)){
        if (object$treetype=="Classification") {
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
            p <- cbind(1,ptemp)[,pos+1,drop=FALSE]
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
        p <- as.numeric(stats::predict(object,newdata=newdata,importance="none",...)$predicted[,2,drop=TRUE])
        ## class(p) <- "predictRisk"
        p
    }else{
        if (object$family=="surv") {
            ptemp <- 1-stats::predict(object,newdata=newdata,importance="none",...)$survival
            pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
            p <- cbind(1,ptemp)[,pos+1,drop=FALSE]
            if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
            p
        }else{
            if (is.character(cause)) cause <- as.numeric(cause)
            if (!is.numeric(cause)) stop("cause is not numeric")
            cif <- stats::predict(object,newdata=newdata,importance="none",...)$cif[,,cause,drop=TRUE]
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
                                          product.limit = TRUE, diag = FALSE, iid = FALSE, average.iid = FALSE, ...) {
    dots <- list(...)
    type <- dots$type
    if(is.null(type)){
        type <- "absRisk"
    }
    store.iid <- dots$store.iid
    
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
                       store.iid = store.iid,
                       type = type)

    out <- outPred[[type]]
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
##' }
##'\dontrun{
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
##' }
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

        newdata$dummy.time=1
        newdata$dummy.event=1
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
    requireNamespace("smcfcs")
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
    newdata <- data.frame(newdata)
    p <- matrix(0, NROW(newdata), length(times))
    term <- attr(terms(as.formula(object$call$formula)), "term.labels")
    sm <- summary(object, newdata = newdata[, term], t = times, start = 0, B = 0) #no confidence interval simulations
    for (i in 1:NROW(newdata)){
        p[i,] <- sm[[i]][,2]
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
        stop("Prediction failed")
    1 - p
}

## * predictRisk.singleEventCB
##' @export
##' @rdname predictRisk
##' @method predictRisk singleEventCB
predictRisk.singleEventCB <- function(object, newdata, times, cause, ...) {
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


#----------------------------------------------------------------------
### predictRisk.R ends here

