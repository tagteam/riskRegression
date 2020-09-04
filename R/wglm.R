### wglm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  1 2020 (14:58) 
## Version: 
## Last-Updated: sep  4 2020 (18:22) 
##           By: Brice Ozenne
##     Update #: 257
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * wglm
#' @title Logistic Regression Using IPCW 
#' @description Logistic regression over multiple timepoints
#' where right-censoring is handle using inverse probability of censoring weighting.
#'
#' @param regressor.event [formula] a formula with empty left hand side and the covariates for the logistic regression on the right hand side.
#' @param formula.censor [formula] a formula used to fit the censoring model.
#' @param times [numeric vector] time points at which to model the probability of experiencing an event.
#' @param data [data.frame] dataset containing the time at which the event occured, the type of event, and regressors used to fit the censoring and logistic models.
#' @param cause [character or numeric] the cause of interest. Defaults to the first cause. 
#' @param fitter [character] routine to fit the Cox regression models.
#' @param product.limit [logical] if \code{TRUE} the survival is computed using the product limit estimator.
#'
#' @details First, a Cox model is fitted (argument formula.censor)
#' and the censoring probabilities are computed relative to each timepoint (argument times) to obtain the censoring weights.
#' Then, for each timepoint, a logistic regression is fitted with the appropriate censoring weights
#' and where the outcome is the indicator of having experience the event of interest (argument cause) at or before the timepoint.
#' 
#' @return an object of class \code{"wglm"}.
#'
#' @examples
#' set.seed(10)
#' n <- 250
#' tau <- 1:5
#' d <- sampleData(n, outcome = "competing.risks")
#' d$Y <- (d$event == 1)*(d$time <= tau[3])
#' d0 <- d[event!=0] ## remove censoring
#'
#' ## no censoring
#' e0.wglm <- wglm(regressor.event = ~ X1, formula.censor = Surv(time,event==0) ~ X1,
#'                times = tau, data = d0)
#' e0.glm <- glm(Y ~ X1, family = binomial, data = d0)
#'
#' ## censoring
#' e.wglm <- wglm(regressor.event = ~ X1, formula.censor = Surv(time,event==0) ~ X1,
#'                times = tau, data = d)
#' @export
wglm <- function(regressor.event, formula.censor, times, data, cause = NA,
                 fitter = "coxph", product.limit = FALSE, ...){
    
    tol <- 1e-12
    varSurv <- riskRegression:::SurvResponseVar(formula.censor)
    newname <- paste0(varSurv$status,".",times)
    IPCW.newname <- paste0("IPCW.",times)
    n.times <- length(times)
    
    ## ** check arguments
    fitter <- match.arg(fitter,c("coxph","cph","phreg"))
    if(any(newname %in% names(data))){
        stop("Argument \'data\' should not have a column named \"",paste0(newname[newname %in% names(data)],collapse="\" \""),"\"\n",
             "This name is used internally \n.")
    }
    if(any(IPCW.newname %in% names(data))){
        stop("Argument \'data\' should not have a column named \"",paste0(IPCW.newname[IPCW.newname %in% names(data)],collapse="\" \""),"\"\n",
             "This name is used internally \n.")
    }
    if(any(is.na(data))){
        warning("Argument \'data\' contains missing values. \n")
    }

    ## ** fit censoring model
    object.censor <- do.call(fitter, list(formula = formula.censor,data = data, x=TRUE,y=TRUE))

    ## ** extract information
    recoverCensor <- object.censor$y[,"status"]
    recoverTime <- object.censor$y[,"time"]
    recoverStatus <- data[[varSurv$status]]
    if(is.na(cause)){
        cause <- sort(recoverStatus[recoverCensor!=1])[1]
    }
    allCauses <- sort(unique(recoverStatus[recoverCensor!=1]))
    n.censor <- sapply(times, function(iTime){sum((recoverTime<iTime)*recoverCensor)})
    
    ## ** create glm object
    out <- list(fit = vector(mode = "list", length = n.times))
    for(iTime in 1:n.times){
        iFormula.glm <- update(regressor.event, paste0(newname[iTime],"~."))
        data[[newname[iTime]]] <- (data[[varSurv$status]]==cause)*(recoverTime<=times[iTime])

        if(n.censor[iTime]>0){
            iPred <- predictRisk(object.censor, diag = TRUE, newdata = data, times = pmin(recoverTime, times[iTime]) - tol,
                                 type = "survival", product.limit = product.limit)
            data[[IPCW.newname[[iTime]]]] <- (1-(recoverCensor==1)*(recoverTime<times[iTime]))/iPred[,1]
        
            suppressWarnings(out$fit[[iTime]] <- do.call(stats::glm, list(formula = iFormula.glm, family = binomial(link = "logit"), data = data, weights = data[[IPCW.newname[[iTime]]]])))
            ## Warning message:
            ##             In eval(family$initialize) : non-integer #successes in a binomial glm!
            out$fit[[iTime]]$time.prior.weights <- pmin(recoverTime, times[iTime]) - tol
        }else{
            IPCW.newname[[iTime]] <- NA
            out$fit[[iTime]] <- do.call(stats::glm, list(formula = iFormula.glm, family = binomial(link = "logit"), data = data))
        }
    }

    ## ** export
    out$call <- match.call()
    out$formula <- update(as.formula(paste0(varSurv$status,"~.")),regressor.event)
    out$n <- NROW(data)
    out$n.censor <- n.censor
    out$cox <- object.censor
    out$data <- data
    out$var.outcome <- varSurv$status
    out$var.time <- varSurv$time
    out$name.IPCW <- IPCW.newname
    out$name.outcome <- newname
    out$times <- times
    out$theCause <- cause
    out$causes <- allCauses
    out$product.limit <- product.limit
    class(out) <- append("wglm",class(out))
    return(out)
}

## * nobs.wglm
nobs.wglm <- function(object){
    return(object$n)
}
## * formula.wglm
formula.wglm <- function(object){
    return(object$formula)
}

## * coef.wglm
coef.wglm <- function(object, times = NULL, simplifies = TRUE, ...){
    if(is.null(times)){
        times <- object$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)
    
    M.coef <- do.call(rbind,setNames(lapply(object$fit, coef),object$time))
    return(M.coef[object$times %in% times,,drop=simplifies])
}

## * summary.wglm
summary.wglm <- function(object, print = TRUE, se = "robust", times = NULL, ...){

    se <- match.arg(se, c("robust","model-wknown","robust-wknown"))
    if(is.null(times)){
        times <- object$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)
    if(se == "robust"){
        object.iid <- lava::iid(object, simplifies = FALSE, times = times)
    }else if(se == "robust-wknown"){
        object.iid <- lapply(object$fit[match(times, object$times)],lava::iid)
    }
    
    out <- setNames(vector(mode = "list", length = n.times), times)
    if(print){
        if(length(object$causes)>1){
            text.CR <- paste0("for cause ",object$theCause)
        }else{
            text.CR <- ""
        }
        if(any(object$n.censor>0)){
            text.IPCW <- "IPCW "
        }else{
            text.IPCW <- ""
        }
        cat("     ",text.IPCW,"logistic regression ",text.CR,": \n",sep="")
    }
    for(iTime in 1:n.times){
        iTime2 <- which(object$times == times[iTime])
        suppressWarnings(out[[iTime]] <- summary(object$fit[[iTime2]])$coef)
        if(se %in% c("robust","robust-wknown")){
            out[[iTime]][,"Std. Error"] <- sqrt(diag(crossprod(object.iid[[iTime]])))
            out[[iTime]][,"z value"] <- out[[iTime]][,"Estimate"]/out[[iTime]][,"Std. Error"]
            out[[iTime]][,"Pr(>|z|)"] <- 2*(1-pnorm(abs(out[[iTime]][,"z value"])))
        }
        if(print){
            cat("----------------------------------------------------------------------------------\n")
            cat("  > time: ",times[iTime],"\n",sep="")
            if(!is.na(eval(object$name.IPCW[[iTime2]]))){                
                print(as.call(list(quote(glm),eval(object$fit[[iTime2]]$call$formula),family = quote(binomial(link="logit")), weights = eval(object$name.IPCW[[iTime2]]))))
            }else{
                print(as.call(list(quote(glm),eval(object$fit[[iTime2]]$call$formula),family = quote(binomial(link="logit")))))
            }
            cat("\n")
            print(out[[iTime]])
            cat("----------------------------------------------------------------------------------\n")
        }
    }

    return(invisible(out))
}

## * print.wglm
print.wglm <- function(object, times = NULL, ...){
    if(is.null(times)){
        times <- object$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)
    
    if(length(object$causes)>1){
        text.CR <- paste0("for cause ",object$theCause)
    }else{
        text.CR <- ""
    }
    if(any(object$n.censor>0)){
        text.IPCW <- "IPCW "
    }else{
        text.IPCW <- ""
    }
    cat("     ",text.IPCW,"logistic regression ",text.CR,": \n",sep="")
    for(iTime in which(times %in% object$times)){
        iTime2 <- which(object$times == times[iTime])
        if(!is.na(eval(object$name.IPCW[[iTime2]]))){                
            object$fit[[iTime2]]$call <- as.call(list(quote(glm),eval(object$fit[[iTime2]]$call$formula),family = quote(binomial(link="logit")), weights = eval(object$name.IPCW[[iTime2]])))
        }else{
            object$fit[[iTime2]]$call <- as.call(list(quote(glm),eval(object$fit[[iTime2]]$call$formula),family = quote(binomial(link="logit"))))
        }
        cat("----------------------------------------------------------------------------------\n")
        cat("  > time: ",times[iTime],"\n",sep="")
        print(object$fit[[iTime2]])
        cat("----------------------------------------------------------------------------------\n")
    }
}

## * score.wglm
## same as in lava
#' @export
score.wglm <- function(x, indiv = FALSE, times = NULL, simplifies = TRUE, ...){
    if(is.null(times)){
        times <- x$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)

    out <- setNames(vector(mode = "list", length = n.times), times)
    for(iTime in 1:n.times){
        iTime2 <- which(x$times == times[iTime])
        iObject <- x$fit[[iTime2]]

        X <- stats::model.matrix(iObject)
        pi <- stats::predict(iObject, type = "response")
        Y <- iObject$y
        W <- iObject$prior.weights
        out[[iTime]] <- colMultiply_cpp(X, scale = W*(Y - pi))
        colnames(out[[iTime]]) <- colnames(X)
        if(indiv==FALSE){out[[iTime]] <- colSums(out[[iTime]])}
    }

    if(n.times == 1 && simplifies){
        return(out[[1]])
    }else{
        return(out)
    }
}

## * information.wglm
## same as in lava
#' @export
information.wglm <- function(x, times = NULL, simplifies = TRUE, ...){
    if(is.null(times)){
        times <- x$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)

    out <- setNames(vector(mode = "list", length = n.times), times)
    for(iTime in 1:n.times){
        iTime2 <- which(x$times == times[iTime])
        iObject <- x$fit[[iTime2]]

        X <- stats::model.matrix(iObject)
        pi <- stats::predict(iObject, type = "response")
        W <- iObject$prior.weights
        out[[iTime]] <- t(colMultiply_cpp(X, scale = W*pi*(1-pi))) %*% X
        rownames(out[[iTime]]) <- colnames(out[[iTime]])
    }

    if(n.times == 1 && simplifies){
        return(out[[1]])
    }else{
        return(out)
    }

}

## * iid.wglm
#' @export
iid.wglm <- function(x, times = NULL, simplifies = TRUE, ...){
    if(is.null(times)){
        times <- x$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)

    ls.score <- score(x, times = times, simplifies = FALSE, indiv = TRUE)
    ls.info <- information(x, times = times, simplifies = FALSE)

    out <- setNames(vector(mode = "list", length = n.times), times)
    for(iTime in 1:n.times){
        iTime2 <- which(x$times == times[iTime])
        iObject <- x$fit[[iTime2]]

        ## ** compute the uncertainty related to the weights
        ## S score, I information matrix, H hessian, X design matrix, Y outcome, \pi fitted probabilities
        ## \beta regresssion parameters (logit) \eta nuisance parameters (censoring)
        ## The estimating equation of \beta is
        ## S = 0 = sum_i W_i X_i (Y_i - \pi_i(\beta))
        ##       = n\int W(O,\Prob_n) X(O) (Y(O) - \pi(O,\beta(\Prob_n))) d\Prob_n(O)
        ## dS(\Prob_h)/dh =   n\int dW(O, \Prob_h)/dh X(O) (Y(O) - \pi(O,\beta(\Prob_h))) d\Prob_h(O)
        ##                  - n\int W(O, \eta(\Prob_h)) X(O) d\pi(O,\beta(\Prob_h))/dh d\Prob_h(O)
        ##                  + n\int W(O, \eta(\Prob_h)) X(O) (Y(O) - \pi(O,\beta(\Prob_h))) d(d\Prob_h(O)/dh)
        ## dS(\Prob_h)/dh =   n(\int dW(O,\Prob_h)/dh X(O) (Y(O) - \pi(O)) d\Prob_h(O))
        ##                  - n(\int W(0)) X(O) d\pi(0,\beta)/d\beta d\Prob_h(O)) * (d\beta(\Prob_h)/d\Prob_h)
        ##                  + n\int W(0) X(O) (Y(O) - \pi(O)) d(d\Prob_n(O)/dh)
        ## \Prob_h = (1-h)\Prob + h \delta_{O_i}
        ## dS(\Prob_h)/dh|h = 0 = n\int IF_{W_O}(O_i) X(O) (Y(O) - \pi(O)) d\Prob_h(O) \IF_\eta(O_i)
        ##                        + dS/d\beta IF_\beta(O_i)
        ##                        + (-S + S_i) where S=0 (since it is at ML)
        ## IF_\beta(O_i) = (S_i + n\int IF_{W_O}(O_i) X(O) (Y(O) - \pi(O)) d\Prob_h(O)) / (-dS/d\beta)
        ##               = (S_i + n*AIF_{W,X*(Y-pi)}(O_i) / I
        ## NOTE: IF_{W_O}(O_i) = dW(O,\eta)/d\eta \IF_\eta(O_i)
        if(x$n.censor[iTime]>0){
            n.obs <- stats::nobs(iObject)
            X <- stats::model.matrix(iObject)
            pi <- stats::predict(iObject, type = "response")
            Y <- iObject$y

            factor <- TRUE        
            attr(factor,"factor") <- lapply(apply(colMultiply_cpp(X, (Y - pi)), 2, list), function(iVec){cbind(iVec[[1]])})
            iPred <- predictRisk(x$cox, diag = TRUE, newdata = x$data, times = iObject$time.prior.weights,
                                 type = "survival", product.limit = FALSE, average.iid = factor, store.iid = "full")
        }
        
        ## ** assemble uncertainty
        ## (S+dS/dW)/I
        if(x$n.censor[iTime]>0){
            out[[iTime]] <- (ls.score[[iTime]] + do.call(cbind,attr(iPred, "average.iid"))*NROW(x$data)) %*% solve(ls.info[[iTime]])
        }else{
            out[[iTime]] <- ls.score[[iTime]] %*% solve(ls.info[[iTime]])
        }
        ## ## WRONG VERSION
        ## factor1 <- colMultiply_cpp(X, (Y - pi)) %*% iVcov
        ## factor2 <- - ls.score[[iTime]] %*% crossprod(iVcov) %*% t(colMultiply_cpp(X, scale = -pi*(1-pi))) %*% X
        
        ## factor <- TRUE        
        ## attr(factor,"factor") <- lapply(apply(factor1 + factor2, 2, list), function(iVec){cbind(iVec[[1]])})
        ## iPred <- predictRisk(x$cox, diag = TRUE, newdata = x$data, times = iObject$time.prior.weights,
        ##                      type = "survival", product.limit = FALSE, average.iid = factor)
        ## out[[iTime]] <- out[[iTime]] - do.call(cbind,attr(iPred, "average.iid"))*NROW(x$data)

    }
    
    ## ** export
    if(n.times == 1 && simplifies){
        return(out[[1]])
    }else{
        return(out)
    }
}


## * predictRisk.wglm
#' @export
#' @rdname predictRisk
#' @method predictRisk wglm
predictRisk.wglm <- function(object, newdata, times = NULL, se = "robust",
                             product.limit = FALSE, diag = FALSE, iid = FALSE, average.iid = FALSE, ...){

    se <- match.arg(se, c("robust","robust-wknown"))

    ## ** extract information and normalize arguments
    if(is.null(times)){
        times <- object$times
    }else{
        if(any(times %in% object$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)
    n.sample <- stats::nobs(object)
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
        n.factor <- NCOL(factor)
    }

    ## hidden argument: enable to ask for the prediction of Y==1 or Y==0
    level <- list(...)$level

    ## ** prepare output
    out <- matrix(NA, nrow = n.newdata, ncol = n.times)
    if(iid){
        attr(out,"iid") <- array(NA, dim = c(n.sample, n.times, n.newdata))
    }
    if(average.iid){
        attr(out,"average.iid") <- lapply(n.factor, function(x){matrix(NA, nrow = n.sample, ncol = n.times)})
        if(!is.null(names(factor))){
            names(attr(out,"average.iid")) <- names(factor)
        }
    }

    if(iid || average.iid){
        if(se=="robust"){
            object.iid <- lava::iid(object, simplifies = FALSE, times = times)
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

######################################################################
### wglm.R ends here
