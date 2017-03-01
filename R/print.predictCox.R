### print.predictCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 15 2017 (17:36) 
## Version: 
## last-updated: Feb 28 2017 (16:27) 
##           By: Thomas Alexander Gerds
##     Update #: 104
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Print predictions from a Cox model
#' @description Print predictions from a Cox model
#' 
#' @param x object obtained with the function \code{predictCox}.
#' @inheritParams predictCox
#' @param ci Logical. If \code{TRUE} display the confidence intervals for the predictions.
#' @param reduce.data Logical. If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param conf.level confidence level of the interval.
#' @param digit integer indicating the number of decimal places.
#' @param ... Passed to print.
#' 
#' @examples
#' library(survival)
#' library(rms)
#' 
#' set.seed(10)
#' d <- sampleData(1e2, outcome = "survival")
#' m.cox <- coxph(Surv(time,event)~ X1 + X2 + X3,
#'                data = d, x = TRUE, y = TRUE)
#' predictCox(m.cox)
#'
#' pred <- predictCox(m.cox, newdata = d[1:5,],
#'                    times = 1:5, type = "survival")
#' pred
#' 
#' pred.data <- predictCox(m.cox, newdata = d[1:4,],
#'              times = 1:5, type = "survival", keep.newdata = TRUE)
#' pred.data
#' print(pred.data, reduce.data = TRUE)
#'
#' pred.ci <- predictCox(m.cox, newdata = d[1:5,],
#'                       times = 1:5,se = TRUE)
#' print(pred.ci, ci = TRUE)
#' 
#' m.cox <- coxph(Surv(time,event)~ strata(X1) + strata(X2) + X3 + X6,
#'                data = d, x = TRUE, y = TRUE)
#' pred.cox <- predictCox(m.cox, newdata = d[c(1:5,10,50),],
#'                        time = 1:5)
#' pred.cox
#' 
#' m.cox <- cph(Surv(time,event)~ strat(X1) + strat(X2) + X3 + X6,
#'              data = d, x = TRUE, y = TRUE)
#' pred.cox <- predictCox(m.cox, newdata = d[c(1:5,10,50),],
#'                        time = 1:5)
#' pred.cox
#' 
#' pred.dataci <- predictCox(m.cox, newdata = d[1:5,],
#'                        times = 1:5, keep.newdata = TRUE, se = TRUE)
#' pred.dataci
#' print(pred.dataci, ci = TRUE)
#'
#' @method print predictCox
#' @export
print.predictCox <- function(x,
                             type = NULL,
                             ci = FALSE,
                             reduce.data = FALSE,
                             conf.level = 0.95,
                             digit = 3, ...){

    ## initialise type
    possible.type <- c("hazard","cumhazard","survival")
    if(is.null(type)){        
        type <- possible.type[possible.type %in% names(x)]
    }else if(any(possible.type %in% names(x) == FALSE)){
        stop(paste(possible.type[possible.type %in% names(x) == FALSE], collapse = " ")," are not in the object \n")
    }

    ## check the presence of the standard errors
    if(ci && any(paste0(type,".se") %in% names(x)==FALSE)){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set argment \'se\' to TRUE when calling predictCox \n")
    }

    ## remove covariates that have the same value for all observations
    newdata <- copy(x$newdata)
    if(!is.null(newdata) && reduce.data){
        test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
        if(any(test)){
            newdata[, (names(test)[test]):=NULL]
        }        
    }

    ## 
    res <- resDisplay <- list()
    
    for(iType in type){
        ls.resPrint <- predict2print(outcome = x[[iType]],
                                     outcome.se = if(ci){x[[paste0(iType,".se")]]}else{NULL},
                                     newdata = newdata,
                                     strata = x$strata,
                                     times = x$times,
                                     digit = digit,
                                     name.outcome = iType,
                                     conf.level = conf.level,
                                     lower = 0, upper = if(iType == "survival"){1}else{Inf})
        res[[iType]] <- ls.resPrint$res
        resDisplay[[iType]] <- ls.resPrint$resDisplay            

    }

    ## rename the element of the list for the display
    names(resDisplay) <- paste0(type, " (columns: ",
                                if(!is.null(x$strata)){"strata|"},
                                if(!is.null(newdata)){"covariate|"},
                                "time, rows: individuals",
                                if(ci){paste0(", [.;.] ",100*conf.level,"% confidence interval)")}else{")"}
                                )

    ## display
    print(noquote(resDisplay), ...)

    ## export
    return(invisible(res))    
}


predict2print <- function(outcome, outcome.se, newdata, strata, times,
                          digit, name.outcome, conf.level, lower, upper){

    n.obs <- NROW(outcome)
    n.time <- NCOL(outcome)
    if(!is.null(times)){
        colnames(outcome) <- times
    }else{
        time.names <- 1:n.time
    }
    
    if(is.null(outcome.se)){
        res <- outcome
        resDisplay <- round(outcome, digit)
    }else{
        resDisplay <- upper <- lower <- matrix(NA, nrow = n.obs, ncol = n.time, dimnames = dimnames(outcome))
            
        lower[] <- pmax(lower,outcome + qnorm((1-conf.level)/2) * outcome.se)
        upper[] <- pmin(upper,outcome + qnorm(1-(1-conf.level)/2) * outcome.se)

        resDisplay[] <- paste0(formatC(outcome, format = "f", digits = digit),
                             " [",formatC(lower, format = "f", digits = digit),
                             ";",formatC(upper, format = "f", digits = digit),"]")

        res <- list(outcome, lower, upper)
        names(res) <- c(name.outcome, "lower", "upper")
    }
    if(!is.null(newdata)){
        index.numeric <- which(unlist(lapply(newdata, is.numeric)))
        
        resDisplay <- cbind(newdata,resDisplay)            
    }
    if(!is.null(strata)){
        resDisplay <- cbind(strata = as.character(strata),resDisplay)            
    }

    return(list(res = res,
                resDisplay = data.table::as.data.table(resDisplay)))
}



#----------------------------------------------------------------------
### print.predictCox.R ends here
