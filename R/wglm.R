### wglm.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  1 2020 (14:58) 
## Version: 
## Last-Updated: sep  3 2020 (11:41) 
##           By: Brice Ozenne
##     Update #: 159
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * wglm
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

    ## ** fit censoring model
    object.censor <- do.call(fitter, list(formula = formula.censor,data = data, x=TRUE,y=TRUE))

    ## ** extract information
    recoverCensor <- object.censor$y[,"status"]
    recoverTime <- object.censor$y[,"time"]
    recoverStatus <- data[[varSurv$status]]
    if(is.na(cause)){
        cause <- sort(recoverStatus[recoverCensor!=1])[1]
    }
    
    ## ** create glm object
    out <- list(fit = vector(mode = "list", length = n.times))
    object.censor <- iidCox(object.censor)
    for(iTime in 1:n.times){
        iPred <- predictRisk(object.censor, diag = TRUE, newdata = data, times = pmin(recoverTime, times[iTime]) - tol,
                             type = "survival", product.limit = product.limit, iid = TRUE, store.iid = "full")
        data[[IPCW.newname[[iTime]]]] <- (1-(recoverCensor==1)*(recoverTime<times[iTime]))/iPred[,1]

        data[[newname[iTime]]] <- (data[[varSurv$status]]==cause)*(recoverTime<=times[iTime])
        iFormula.glm <- update(regressor.event, paste0(newname,"~."))
        ## Warning message:
        ##             In eval(family$initialize) : non-integer #successes in a binomial glm!
        suppressWarnings(out$fit[[iTime]] <- do.call(stats::glm, list(formula = iFormula.glm, family = binomial(link = "logit"), data = data, weights = data[[IPCW.newname[[iTime]]]])))
        out$fit[[iTime]]$iid.prior.weights <- attr(iPred,"iid")
        out$fit[[iTime]]$time.prior.weights <- pmin(recoverTime, times[iTime]) - tol
    }

    ## ** export
    out$cox <- object.censor
    out$data <- data
    out$name.IPCW <- IPCW.newname
    out$name.outcome <- newname
    out$times <- times
    out$cause <- cause
    out$product.limit <- product.limit
    class(out) <- append("wglm",class(out))
    return(out)
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
summary.wglm <- function(object, print = TRUE, times = NULL, ...){
    if(is.null(times)){
        times <- object$times
    }else{
        if(any(times %in% x$times == FALSE)){
            stop("Incorrect specification of argument \'times\' \n",
                 "Should be one of \"",paste0(times,collapse="\" \""),"\" \n")
            
        }
    }
    n.times <- length(times)

    out <- setNames(vector(mode = "list", length = n.times), times)
    if(print){cat("     IPCW logistic regression: \n")}
    for(iTime in 1:n.times){
        iTime2 <- which(object$times == times[iTime])
        suppressWarnings(out[[iTime]] <- summary(object$fit[[iTime2]]))
        out[[iTime]]$call <- as.call(list(quote(glm),eval(out[[iTime]]$call$formula),family = quote(binomial(link="logit")), weights = eval(object$name.IPCW[[iTime2]])))
        if(print){
            cat("----------------------------------------------------------------------------------\n")
            cat("  > time: ",times[iTime],"\n",sep="")
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

    cat("     IPCW logistic regression: \n")
    for(iTime in which(times %in% object$times)){
        object$fit[[iTime]]$call <- as.call(list(quote(glm),eval(object$fit[[iTime]]$call$formula),family = quote(binomial(link="logit")), weights = eval(object$name.IPCW[[iTime]])))
        cat("----------------------------------------------------------------------------------\n")
        cat("  > time: ",times[iTime],"\n",sep="")
        print(object$fit[[iTime]])
        cat("----------------------------------------------------------------------------------\n")
    }
}

## * score.wglm
## same as in lava
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
        out[[iTime]] <- t(colMultiply_cpp(X, scale = -W*pi*(1-pi))) %*% X
        rownames(out[[iTime]]) <- colnames(out[[iTime]])
    }

    if(n.times == 1 && simplifies){
        return(out[[1]])
    }else{
        return(out)
    }

}

## * iid.wglm
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

        ## ** extract information
        n.obs <- stats::nobs(iObject)
        X <- stats::model.matrix(iObject)
        pi <- stats::predict(iObject, type = "response")
        Y <- iObject$y

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

        factor <- TRUE        
        attr(factor,"factor") <- lapply(apply(colMultiply_cpp(X, (Y - pi)), 2, list), function(iVec){cbind(iVec[[1]])})
        iPred <- predictRisk(x$cox, diag = TRUE, newdata = x$data, times = iObject$time.prior.weights,
                             type = "survival", product.limit = FALSE, average.iid = factor)

        ## ** assemble uncertainty
        ## (S+dS/dW)/I
        out[[iTime]] <- (ls.score[[iTime]] + do.call(cbind,attr(iPred, "average.iid"))*NROW(x$data)) %*% solve(ls.info[[iTime]])

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

######################################################################
### wglm.R ends here
