### autoplot.predictCSC.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 27 2017 (10:47) 
## Version: 
## last-updated: Apr 27 2026 (14:03) 
##           By: Brice Ozenne
##     Update #: 162
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * autoplot.predictCSC (documentation)
#' @title Plot Predictions From a Cause-specific Cox Proportional Hazard Regression
#' @description Plot predictions from a Cause-specific Cox proportional hazard regression.
#' @name autoplot.predictCSC
#' 
#' @param object,x Object obtained with the function \code{predictCox}.
#' @param ci [logical] If \code{TRUE} display the confidence intervals for the predictions.
#' @param band [logical] If \code{TRUE} display the confidence bands for the predictions.
#' @param group.by [character] The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}. 
#' @param reduce.data [logical] If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param plot [logical] Should the graphic be plotted.
#' @param digits [integer] Number of decimal places.
#' @param smooth [logical] Should a smooth version of the risk function be plotted instead of a simple function?
#' @param alpha [numeric, 0-1] Transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param atRisk [logical] Should the number at risk be diplayed at the bottom of the plot? 
#' @param xlim,ylim [numeric vector of length 2] limits for the x and y axes.
#' @param ... Additional parameters to cutomize the display.
#' 
#' @return Invisible. A list containing:
#' \itemize{
#' \item plot: the ggplot object.
#' \item data: the data used to create the plot.
#' }
#'
#' @seealso
#' \code{\link{predict.CauseSpecificCox}} to compute risks based on a CSC model.

## * autoplot.predictCSC (examples)
#' @examples
#' library(survival)
#' library(rms)
#' library(ggplot2)
#' library(prodlim)
#' 
#' #### simulate data ####
#' set.seed(10)
#' d <- sampleData(1e2, outcome = "competing.risks")
#' seqTau <- c(0,unique(sort(d[d$event==1,time])), max(d$time))
#' 
#' #### CSC model ####
#' m.CSC <- CSC(Hist(time,event)~ X1 + X2 + X6, data = d)
#' 
#' pred.CSC <- predict(m.CSC, newdata = d[1:2,], time = seqTau, cause = 1, band = TRUE)
#' autoplot(pred.CSC, alpha = 0.2)
#' 
#' #### stratified CSC model ####
#' m.SCSC <- CSC(Hist(time,event)~ strata(X1) + strata(X2) + X6,
#'               data = d)
#' pred.SCSC <- predict(m.SCSC, time = seqTau, newdata = d[1:4,],
#'                      cause = 1, keep.newdata = TRUE, keep.strata = TRUE)
#' autoplot(pred.SCSC, group.by = "strata")

## * autoplot.predictCSC (code)
#' @rdname autoplot.predictCSC
#' @method autoplot predictCSC
#' @export
autoplot.predictCSC <- function(object,
                                ci = object$se,
                                band = object$band,
                                plot = TRUE,
                                smooth = FALSE,
                                digits = 2,
                                alpha = NA,
                                group.by = "row",
                                reduce.data = FALSE,
                                atRisk = FALSE,
                                xlim = NULL,
                                ylim = NULL,
                                ...){
  
    ## initialize and check
    group.by <- match.arg(group.by, c("row","covariates","strata", names(object$newdata)))
  
    if(group.by[[1]] == "covariates" && ("newdata" %in% names(object) == FALSE)){
        stop("argument \'group.by\' cannot be \"covariates\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling the predictCox function \n")
    }
    if(group.by[[1]] == "strata" && ("strata" %in% names(object) == FALSE)){
        stop("argument \'group.by\' cannot be \"strata\" when strata is missing in the object \n",
             "set argment \'keep.strata\' to TRUE when calling the predictCox function \n")
    }
  
    if(ci[[1]] && (object$se[[1]]==FALSE || is.null(object$conf.level))){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set arguments \'se\' and \'confint\' to TRUE when calling the predict.CauseSpecificCox function \n")
    }
    if(band[[1]] && (object$band[[1]]==FALSE  || is.null(object$conf.level))){
        stop("argument \'band\' cannot be TRUE when the quantiles for the confidence bands have not been computed \n",
             "set arguments \'band\' and \'confint\' to TRUE when calling the predict.CauseSpecificCox function \n")
    }
    if(!object$baseline & any(rank(object$times) != 1:length(object$times))){
        stop("Invalid object. The prediction times must be strictly increasing \n")
    }
    type <- ifelse("absRisk" %in% names(object),"absRisk","survival")
    ## dots <- list(...)
    ## if(length(dots)>0){
    ##     txt <- names(dots)
    ##     txt.s <- if(length(txt)>1){"s"}else{""}
    ##     stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    ## }
    
    ## ** reshape data
    if(object$baseline){

        object$absRisk <- rbind(object[[type]])
        if(ci){
            object[[paste0(type,".lower")]] <- rbind(object[[paste0(type,".lower")]])
            object[[paste0(type,".upper")]] <- rbind(object[[paste0(type,".upper")]])
        }
        if(band){
            object[[paste0(type,".lowerBand")]] <- rbind(object[[paste0(type,".lowerBand")]])
            object[[paste0(type,".upperBand")]] <- rbind(object[[paste0(type,".upperBand")]])
        }
        
        ## *** add time 0
        if((0 %in% object$times == FALSE) && (is.null(xlim) || xlim[1]<=0)){
                
            n.strata <- ifelse(is.null(object$strata), 1 ,length(levels(object$strata)))
            object[[type]] <- cbind(matrix(type=="survival", nrow = 1, ncol = n.strata),
                                    object[[type]]) ## 0 if absRisk 1 if survival
            if(ci){
                object[[paste0(type,".lower")]] <- cbind(matrix(type=="survival", nrow = 1, ncol = n.strata),
                                                         object[[paste0(type,".lower")]])
                object[[paste0(type,".upper")]] <- cbind(matrix(type=="survival", nrow = 1, ncol = n.strata),
                                                         object[[paste0(type,".upper")]])
            }
            if(band){
                object[[paste0(type,".lowerBand")]] <- cbind(matrix(type=="survival", nrow = 1, ncol = n.strata),
                                                             object[[paste0(type,".lowerBand")]])
                object[[paste0(type,".upperBand")]] <- cbind(matrix(type=="survival", nrow = 1, ncol = n.strata),
                                                             object[[paste0(type,".upperBand")]])
            }
            object$times <- c(rep(0, n.strata), object$times)
            if(!is.null(object$strata)){
                object$strata <- c(factor(levels(object$strata), levels = levels(object$strata)), object$strata)
            }                

            if(!is.null(object$newdata)){
                newdata0 <- object$newdata[1:n.strata,,drop=FALSE]
                newdata0$time <- 0
                newdata0$status <- NA
                newdata0$event <- factor(NA, levels = levels(newdata$status))
                if(!is.null(object$strata)){
                    newdata0$strata <- factor(levels(object$strata), levels = levels(object$strata))
                    newdata0$strata.num <- 1:n.strata
                }
                object$newdata <- rbind(newdata0, object$newdata)
            }
            
            }
            ## *** extrapolate beyond last observations when last observation is an event
            ## NOTE: survival may not be 0 when using exponential approximation
            time.beyond <- max(c(xlim, max(object$time)))
            if(is.null(object$strata)){
                time.maxStrata <- max(object$time)
            }else{
                time.maxStrata <- tapply(object$time, object$strata, max)
            }            
            if(any(time.maxStrata < time.beyond & time.beyond < object$lastEventTime)){

                ## identify which strata needs update and what is the last observation in each strata
                if(is.null(object$strata)){
                    strata.beyond <- 1
                    n.strataBeyond <- 1
                    index.strataBeyond <-  length(object$time)
                }else{
                    strata.beyond <- which(time.maxStrata < time.beyond & time.beyond < object$lastEventTime)
                    n.strataBeyond <- length(strata.beyond)
                    index.strataBeyond <-  which(rev(!duplicated(rev(object$strata))))[strata.beyond]
                }
                ## update: add a ficticious observation
                object[[type]] <- cbind(object[[type]], matrix(NA, nrow = 1, ncol = n.strataBeyond)) ## NA instead of object[[type]][,index.strataBeyond] so no point is displayed (but the plot function will extend the line)
                if(ci){
                    object[[paste0(type,".lower")]] <- cbind(object[[paste0(type,".lower")]],
                                                             matrix(object[[paste0(type,".lower")]][,index.strataBeyond], nrow = 1, ncol = n.strataBeyond))
                    object[[paste0(type,".upper")]] <- cbind(object[[paste0(type,".upper")]],
                                                             matrix(object[[paste0(type,".lower")]][,index.strataBeyond], nrow = 1, ncol = n.strataBeyond))
                }
                if(band){
                    object[[paste0(type,".lowerBand")]] <- cbind(object[[paste0(type,".lowerBand")]],
                                                                 matrix(object[[paste0(type,".lowerBand")]][,index.strataBeyond], nrow = 1, ncol = n.strataBeyond))
                    object[[paste0(type,".upperBand")]] <- cbind(object[[paste0(type,".upperBand")]],
                                                                 matrix(object[[paste0(type,".upperBand")]][,index.strataBeyond], nrow = 1, ncol = n.strataBeyond))
                }
                object$times <- c(object$times, rep(time.beyond, n.strataBeyond))
                if(!is.null(object$strata)){
                    object$strata <- c(object$strata, factor(levels(object$strata)[strata.beyond], levels = levels(object$strata)))
                }                

                if(!is.null(object$newdata)){
                    newdataF <- object$newdata[1:n.strataBeyond]
                    newdataF$time <- 0
                    newdataF$status <- NA
                    newdataF$event <- factor(NA, levels = levels(newdata$status))
                    if(!is.null(object$strata)){
                        newdataF$strata[] <- factor(levels(object$strata)[strata.beyond], levels = levels(object$strata))
                        newdataF$strata.num <- strata.beyond
                    }
                    object$newdata <- rbind(object$newdata, newdataF)
                }
            }

            if(ci && !is.null(object$newdata) && !is.na(alpha)){
                object[[paste0(type,".lower")]][,is.na(object$newdata$status)] <- NA
                object[[paste0(type,".upper")]][,is.na(object$newdata$status)] <- NA
                object[[paste0(type,".lower")]][,object$newdata$status==0] <- NA
                object[[paste0(type,".upper")]][,object$newdata$status==0] <- NA
            }
            
            newdata <- NULL

    }else{
    
        newdata <- copy(object$newdata)
        if(!is.null(newdata) && reduce.data[[1]]){
            test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
            if(any(test)){
                newdata[, (names(test)[test]):=NULL]
            }
            status <- NULL
        }
    }

    dataL <- predict2melt(outcome = object[[type]], ci = ci, band = band,
                          outcome.lower = if(ci){object[[paste0(type,".lower")]]}else{NULL},
                          outcome.upper = if(ci){object[[paste0(type,".upper")]]}else{NULL},
                          outcome.lowerBand = if(band){object[[paste0(type,".lowerBand")]]}else{NULL},
                          outcome.upperBand = if(band){object[[paste0(type,".upperBand")]]}else{NULL},
                          newdata = newdata,
                          status = object$newdata$event,
                          strata = object$strata,
                          times = object$times,
                          name.outcome = type,
                          group.by = group.by,
                          digits = digits,
                          diag = object$diag,
                          baseline = object$baseline
                          )

    ## ** display
    if(object$baseline & length(object$lp.var)>0){
        ylab <- ifelse(type=="absRisk","Baseline absolute risk","Baseline survival")
    }else{
        ylab <- ifelse(type=="absRisk","Absolute risk","Survival")
    }

    gg.res <- predict2plot(dataL = dataL,
                           name.outcome = type, # must not contain space to avoid error in ggplot2
                           ci = ci,
                           band = band,
                           group.by = group.by,
                           conf.level = object$conf.level,
                           alpha = alpha,
                           smooth = smooth,
                           xlab = "time",
                           ylab = ylab,
                           atRisk = atRisk,
                           xlim = xlim,
                           ylim = ylim,                           
                           ...
                           )
      

    if(object$baseline){
        if(length(object$lp.var)>0){
            gg.res$plot <- gg.res$plot + ggplot2::ggtitle(paste0("Baseline: ",paste(paste(names(object$lp.var), object$lp.var, sep = "="), collapse = ", ")))
        }
        if(!is.null(object$strata)){
            gg.res$plot <- gg.res$plot + ggplot2::labs(color = "Strata")
        }
    }

  if(plot){
    print(gg.res$plot)
  }
  
  return(invisible(gg.res))
}

#----------------------------------------------------------------------
### autoplot.predictCSC.R ends here
