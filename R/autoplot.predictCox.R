### autoplot.predictCox.R --- 
##----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 17 2017 (10:06) 
## Version: 
## last-updated: apr 27 2026 (18:24) 
##           By: Brice Ozenne
##     Update #: 1955
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * autoplot.predictCox (documentation)
#' @title Plot Predictions From a Cox Model
#' @description Plot predictions from a Cox model.
#' The method \code{plot} is a warper around \code{autoplot} that automatically produces a graphical display
#' but hides the corresponding dataset and ggplot2 object using \code{invisible}.
#' @name autoplot.predictCox
#'  
#' @param object,x Object obtained with the function \code{predictCox}.
#' @param type [character] The type of predicted value to display.
#' Choices are:
#' \code{"hazard"} the hazard function,
#' \code{"cumhazard"} the cumulative hazard function, 
#' or \code{"survival"} the survival function.
#' @param ci [logical] If \code{TRUE} display the confidence intervals for the predictions.
#' @param band [logical] If \code{TRUE} display the confidence bands for the predictions.
#' @param digits [integer] Number of decimal places when displaying the values of the covariates in the caption.
#' @param alpha [numeric, 0-1] Transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param ylab [character] Label for the y axis.
#' @param smooth [logical] Should a smooth version of the risk function be plotted instead of a simple function?
#' @param first.derivative [logical] If \code{TRUE}, display the first derivative over time of the risks/risk differences/risk ratios.
#' (confidence intervals are obtained via simulation).
#' @param group.by [character] The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}.
#' @param reduce.data [logical] If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param atRisk [logical] Should the number at risk be diplayed at the bottom of the plot? 
#' @param xlim,ylim [numeric vector of length 2] limits for the x and y axes.
#' @param ... Additional parameters to cutomize the display: \itemize{
#' \item size.estimate: thickness of the line used to represent the estimated survival/hazard.
#' \item size.point: size of the points used to represent obsevations (if any).
#' \item size.ci: thickness of the line used to represent the confidence intervals for the survival/hazard.
#' \item size.band: thickness of the line used to represent the confidence bands for the survival/hazard.
#' \item size.atRisk: font size used to display the number at risk.
#' \item space.atRisk: vertical space used between the  display of the number at risk across strata. 
#' \item shape.point: shape of the points used to represent obsevations (if any).
#' \item col.ribbon: should shaded areas be colored according to the grouping factor. Otherwise grey. 
#' }
#'
#' @return A (possibly invisible) list containing:
#' \itemize{
#' \item plot: the ggplot object.
#' \item data: the data used to create the plot.
#' }
#' 
#' @seealso
#' \code{\link{predictCox}} to compute cumulative hazard and survival based on a Cox model.


## * autoplot.predictCox (examples)
#' @examples
#' library(survival)
#' library(ggplot2)
#'
#' #### simulate data ####
#' set.seed(10)
#' d <- sampleData(1e2, outcome = "survival")
#' seqTau <- c(0,sort(unique(d$time[d$event==1])), max(d$time))
#' 
#' #### Cox model ####
#' m.cox <- coxph(Surv(time,event)~ X1 + X2 + X3,
#'                 data = d, x = TRUE, y = TRUE)
#'
#' ## display baseline hazard
#' e.basehaz <- predictCox(m.cox)
#' plot(e.basehaz, type = "cumhazard")
#' \dontrun{
#' plot(e.basehaz, type = "cumhazard", size.point = 0) ## without points
#' plot(e.basehaz, type = "cumhazard", smooth = TRUE)
#' plot(e.basehaz, type = "cumhazard", smooth = TRUE, first.derivative = TRUE)
#' }
#'
#' ## display predicted survival
#' \dontrun{
#' pred.cox <- predictCox(m.cox, newdata = d[1:2,],
#'   times = seqTau, type = "survival", keep.newdata = TRUE)
#' plot(pred.cox)
#' plot(pred.cox, smooth = TRUE)
#' plot(pred.cox, group.by = "covariates")
#' plot(pred.cox, group.by = "covariates", reduce.data = TRUE)
#' plot(pred.cox, group.by = "X1", reduce.data = TRUE)
#' }
#' 
#' ## predictions with confidence interval/bands
#' \dontrun{
#' pred.cox <- predictCox(m.cox, newdata = d[1:2,,drop=FALSE],
#'   times = seqTau, type = "survival", band = TRUE, se = TRUE, keep.newdata = TRUE)
#' res <- plot(pred.cox, ci = TRUE, band = TRUE)
#' res$plot + facet_wrap(~row)
#' res2 <- plot(pred.cox, ci = TRUE, band = TRUE, alpha = NA)
#' res2$plot + facet_wrap(~row, labeller = label_both) + guides(color = "none")
#' }
#' 
#' #### Stratified Cox model ####
#' \dontrun{
#' m.cox.strata <- coxph(Surv(time,event)~ strata(X1) + strata(X2) + X3 + X4,
#'                       data = d, x = TRUE, y = TRUE)
#'
#' ## baseline hazard
#' pred.baseline <- predictCox(m.cox.strata, keep.newdata = TRUE, type = "survival")
#' res <- plot(pred.baseline)
#' res$plot + facet_wrap(~strata, labeller = label_both)
#'
#' ## predictions
#' pred.cox.strata <- predictCox(m.cox.strata, newdata = d[1:3,,drop=FALSE],
#'                               time = seqTau, keep.newdata = TRUE, se = TRUE)
#'
#' res2 <- plot(pred.cox.strata, type = "survival", group.by = "strata")
#' res2$plot + facet_wrap(~strata, labeller = label_both) + theme(legend.position="bottom")
#' 
#' ## smooth version
#' plot(pred.cox.strata, type = "survival", group.by = "strata", smooth = TRUE, ci = FALSE)
#' }
#' 
#' #### Cox model with splines ####
#' \dontrun{
#' require(splines)
#' m.cox.spline <- coxph(Surv(time,event)~ X1 + X2 + ns(X6,4),
#'                 data = d, x = TRUE, y = TRUE)
#' grid <- data.frame(X1 = factor(0,0:1), X2 = factor(0,0:1),
#'                    X6 = seq(min(d$X6),max(d$X6), length.out = 100))
#' pred.spline <- predictCox(m.cox.spline, newdata = grid, keep.newdata = TRUE,
#'                           se = TRUE, band = TRUE, centered = TRUE, type = "lp")
#' plot(pred.spline, group.by = "X6")
#' plot(pred.spline, group.by = "X6", alpha = 0.5, col.ribbon = TRUE)
#' 
#' grid2 <- data.frame(X1 = factor(1,0:1), X2 = factor(0,0:1),
#'                     X6 = seq(min(d$X6),max(d$X6), length.out = 100))
#' pred.spline <- predictCox(m.cox.spline, newdata = rbind(grid,grid2), keep.newdata = TRUE,
#'                           se = TRUE, band = TRUE, centered = TRUE, type = "lp")
#' autoplot(pred.spline, group.by = c("X6","X1"), alpha = 0.5)$plot + facet_wrap(~X1)
#' }

## * autoplot.predictCox (code)
#' @rdname autoplot.predictCox
#' @method autoplot predictCox
#' @export
autoplot.predictCox <- function(object,
                                type = NULL,
                                ci = object$se,
                                band = object$band,
                                smooth = NULL,
                                digits = 2,
                                alpha = NA,
                                group.by = NULL,
                                reduce.data = FALSE,
                                ylab = NULL,
                                first.derivative = FALSE,
                                atRisk = FALSE,
                                xlim = NULL,
                                ylim = NULL,
                                ...){
  
    ## initialize and check    
    possibleType <- c("cumhazard","survival","lp")
    possibleType <- possibleType[possibleType %in% names(object)]
    nVar.lp <- length(object$var.lp)
        
    if(is.null(type)){
        if(length(possibleType) == 1){
            type <- possibleType
        }else{
            type <- utils::tail(possibleType,1)
        }
    }else{
        type <- match.arg(type, possibleType)  
    }
    if(is.null(smooth)){
        if(type == "lp"){
            smooth <- 0.5
        }else{
            smooth <- FALSE
        }
    }
    if(is.null(ylab)){
        if(first.derivative){
            ylab <- switch(type,
                           "cumhazard" = if(object$baseline && nVar.lp>0){"Instantaneous baseline hazard"}else{"Instantaneous hazard"},
                           "survival" = if(object$baseline && nVar.lp>0){"Derivative of the baseline survival"}else{"Derivative of the survival"})
        }else{
            ylab <- switch(type,
                           "lp" = "linear predictor",
                           "cumhazard" = if(object$baseline && nVar.lp>0){"Cumulative baseline hazard"}else{"Cumulative hazard"},
                           "survival" = if(object$baseline && nVar.lp>0){"Baseline survival"}else{"Survival"})
        }
    }

    if(type=="lp"){
        group.by <- match.arg(group.by, object$var.lp, several.ok = TRUE)
    }else{
        if(is.null(group.by)){
            group.by <- ifelse(object$baseline & length(object$var.strata)>0, "strata", "row")
        }else if(length(group.by)>1){
            stop("Argument \'group.by\' must have length 1.\n")
        }else{
            group.by <- match.arg(group.by, c("row","covariates","strata",object$var.lp,object$var.strata))
        }
    }
    

    if(group.by[[1]] == "covariates" && ("newdata" %in% names(object)) == FALSE){
        stop("Argument \'group.by\' cannot be \"covariates\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling the predictCox function \n")
    }

    if(group.by[[1]] %in% c("strata") && ("strata" %in% names(object) == FALSE)){
        stop("Argument \'group.by\' cannot be \"strata\" when strata is missing in the object \n",
             "set argment \'keep.strata\' to TRUE when calling the predictCox function \n")
    }
    if(group.by[[1]] %in% c(object$var.lp,object$var.strata) && ("newdata" %in% names(object) == FALSE)){
        stop("Argument \'group.by\' cannot be \"",group.by[[1]],"\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling the predictCox function \n")
    }
    if(ci[[1]]==TRUE && (object$se[[1]]==FALSE || is.null(object$conf.level))){
        stop("Argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set arguments \'se\' and \'confint\' to TRUE when calling the predictCox function \n")
    }
    if(band[[1]] && (object$band[[1]]==FALSE  || is.null(object$conf.level))){
        stop("Argument \'band\' cannot be TRUE when the quantiles for the confidence bands have not been computed \n",
             "set arguments \'band\' and \'confint\' to TRUE when calling the predictCox function \n")
    }
    if(object$nTimes!= 0 && any(rank(object$times) != 1:length(object$times))){
        stop("Invalid object. The prediction times must be strictly increasing \n")
    }
    dots <- list(...)
    if("shape.point" %in% names(dots) && !object$baseline){
        message("Argument \'shape.points\' is ignored when predictCox has been called with \'newdata\' argument. \n",
                "Only relevant when displaying 'baseline' cumulative hazard or 'baseline' survival. \n")
    }else  if("shape.point" %in% names(dots) && is.null(object$newdata)){
        message("Argumet \'shape.points\' is ignored. \n",
                "Consider setting the argument \"keep.newdata\" to TRUE when calling predictCox. \n")
    }
    if(any(atRisk>0) & is.null(object$newdata)){
        if(!object$baseline){
            message("Argument \'atRisk\' is only relevant when displaying baseline cumulative hazard or survival. \n",
                    "(i.e. when argument \'newdata\' was not used when calling predictCox) \n")
            atRisk <- FALSE
        }else{
            message("Argument \'atRisk\' is ignored unless argument \'keep.newdata\' was set to TRUE when calling predictCox. \n")
            atRisk <- FALSE
        }
    }
    
    ## ** reshape data
    if(type == "lp"){ ## linear predictor
        if(first.derivative){
            stop("Argument \'first.derivative\' should be FALSE when argument \'type\' equals \"lp\". \n")
        }
        if(length(object$var.lp)==0){
            stop("No covariate in the linear predictor so nothing to display.\n")
        }
        if(group.by[1] %in% object$var.lp == FALSE){
            stop("The first element of argument \'group.by\' should refer to one of the covariates: \"",paste(object$var.lp, collapse= "\" \""),"\".\n")
        }
        if("time" %in% group.by){
            stop("The argument \'group.by\' should not contain \"time\" as this name is used internally.\n")
        }
        if("row" %in% group.by){
            stop("The argument \'group.by\' should not contain \"row\" as this name is used internally.\n")
        }
        if("lowerCI" %in% group.by){
            stop("The argument \'group.by\' should not contain \"lowerCI\" as this name is used internally.\n")
        }
        if("upperCI" %in% group.by){
            stop("The argument \'group.by\' should not contain \"upperCI\" as this name is used internally.\n")
        }
        dataL <- data.table::as.data.table(object)
        dataL$lp.smooth <- dataL$lp
        dataL[["time"]] <- dataL[[group.by[1]]]
        dataL[["row"]] <- 1
        if(ci){
            data.table::setnames(dataL, old = c("lp.lower","lp.upper"), new = c("lowerCI","upperCI"))
            dataL$lowerCI.smooth <- dataL$lowerCI
            dataL$upperCI.smooth <- dataL$upperCI
        }
        if(band){
            data.table::setnames(dataL, old = c("lp.lowerBand","lp.upperBand"), new = c("lowerBand","upperBand"))
            dataL$lowerBand.smooth <- dataL$lowerBand
            dataL$upperBand.smooth <- dataL$upperBand
        }
        object$infoVar$time <- group.by[1]
        group.by <- if(length(group.by[-1])==0){"row"}else{group.by[-1]}
    }else{

        if(object$baseline){ ## baseline hazard/survival
            
            object[[type]] <- rbind(object[[type]])
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
                                        object[[type]]) ## 0 if cumhazard 1 if survival
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
                    newdata0 <- as.data.frame(object$newdata)[1:n.strata,,drop=FALSE]
                    newdata0$start <- 0
                    newdata0$stop <- 0
                    newdata0$status <- NA
                    if(!is.null(object$strata)){
                        newdata0$strata[] <- factor(levels(object$strata), levels = levels(object$strata))
                        newdata0$strata.num <- 1:n.strata
                    }
                    newdata0$eXb <- NA
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
                    newdataF$start <- 0
                    newdataF$stop <- time.beyond
                    newdataF$status <- NA
                    if(!is.null(object$strata)){
                        newdataF$strata <- factor(levels(object$strata)[strata.beyond], levels = levels(object$strata))
                        newdataF$strata.num <- strata.beyond
                    }
                    newdataF$eXb <- NA
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
            newdata <- data.table::copy(object$newdata) ## can be NULL
            if(!is.null(newdata) && reduce.data[[1]]==TRUE){
                test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
                if(any(test)){
                    newdata[, (names(test)[test]):=NULL]
                }        
            }
            status <- NULL
            if(first.derivative && ci){
                if(is.null(object[[paste0(type,".iid")]])){
                    stop("Set argument \'iid\' to TRUE when calling predictCox to be able to display confidence intervals for the first derivative of the ",type,".\n")
                }
                iIID <- object[[paste0(type,".iid")]]
                ## compute variance-covariance matrix
                attr(first.derivative,"vcov") <- lapply(1:(dim(iIID)[3]), function(iObs){
                    if(dim(iIID)[2]==1){
                        if(type=="lp"){
                            return(sum(iIID[,iObs,]^2))
                        }else{
                            return(sum(iIID[,,iObs]^2))
                        }
                    }else{ ## crossprod among timepoints where the 'survival' is not only NA 
                        return(crossprod(iIID[,,iObs][,colSums(!is.na(iIID[,,iObs]))>0,drop=FALSE]))
                    }
                })
            }
        }

        dataL <- predict2melt(outcome = object[[type]], ci = ci, band = band,
                              outcome.lower = if(ci){object[[paste0(type,".lower")]]}else{NULL},
                              outcome.upper = if(ci){object[[paste0(type,".upper")]]}else{NULL},
                              outcome.lowerBand = if(band){object[[paste0(type,".lowerBand")]]}else{NULL},
                              outcome.upperBand = if(band){object[[paste0(type,".upperBand")]]}else{NULL},
                              newdata = newdata,
                              status = object$newdata$status,
                              strata = object$strata,
                              times = object$times,
                              name.outcome = type,
                              group.by = group.by,
                              digits = digits,
                              diag = object$diag,
                              baseline = object$baseline
                              )
    }

    ## add covariate names to the caption 
    if(group.by == "strata" && !is.null(object$var.strata)){
        attr(group.by,"label") <- paste0(group.by," (",paste(object$var.strata, collapse = ", "),")")
        dataL[ , (attr(group.by,"label")) := .SD$strata]
    }else if(group.by ==  "covariates"){
        attr(group.by, "label") <- attr(dataL, "covariates")
    }else{
        attr(group.by, "label") <- group.by
    }

    ## ** display
    gg.res <- predict2plot(dataL = dataL,
                           name.outcome = type,
                           ci = ci,
                           band = band,
                           group.by = group.by,
                           conf.level = object$conf.level,
                           alpha = alpha,
                           smooth = smooth,
                           xlab = if(is.null(object$infoVar)){"time"}else{object$infoVar$time},
                           ylab = ylab,
                           first.derivative = first.derivative,
                           atRisk = atRisk,
                           xlim = xlim,
                           ylim = ylim,                           
                           ...
                           )

    if(object$baseline){
        if(nVar.lp>0){
            vec.center <- sapply(attr(object$var.lp,"center"), function(iE){ifelse(is.numeric(iE),round(iE,digits),as.character(iE))})
            gg.res$plot <- gg.res$plot + ggplot2::ggtitle(paste0("Baseline: ",paste(paste(names(vec.center), vec.center, sep = "="), collapse = ", ")))
        }
        if(is.null(object$strata) && (ci || band)){
            gg.res$plot <- gg.res$plot + ggplot2::guides(color = "none")
        }
    }
  
    return(gg.res)
}

## * predict2melt
predict2melt <- function(outcome, name.outcome,
                         ci, outcome.lower, outcome.upper,
                         band, outcome.lowerBand, outcome.upperBand,
                         newdata, status, strata, times, 
                         group.by, digits, diag, baseline){

    patterns <- NULL ## [:CRANtest:] data.table

    ## ** unique name for each observation
    if(!is.null(strata) && baseline){ 
        n.strata <- length(levels(strata))
        index.obs <- unlist(tapply(1:length(strata), strata, FUN = identity, simplify = FALSE))
        index.strata <- unlist(tapply(strata, strata, FUN = function(iVec){1:length(iVec)}, simplify = FALSE))
        name.obs <- paste0("strata",formatC(as.numeric(strata), width = log10(n.strata) + 1, flag = 0),
                           "time",formatC(index.strata[order(index.obs)], width = log10(length(strata)) + 1, flag = 0))
    }else{
        name.obs <- paste0("time",formatC(1:length(times), width = log10(length(times)) + 1, flag = 0))
    }

    ## ** move estimate, CI, Band to a long format
    pattern <- paste0(name.outcome,"_")
    if(diag){
        outcome <- t(outcome)
    }
    colnames(outcome) <- paste0(name.outcome,"_",name.obs)
    n.obs <- NROW(outcome)

    if(ci){
        if(diag){
            outcome.lower <- t(outcome.lower)
            outcome.upper <- t(outcome.upper)
        }
        pattern <- c(pattern,"lowerCI_","upperCI_")
        colnames(outcome.lower) <- paste0("lowerCI_",name.obs)
        colnames(outcome.upper) <- paste0("upperCI_",name.obs)
    }
    if(band){
        if(diag){
            outcome.lowerBand <- t(outcome.lowerBand)
            outcome.upperBand <- t(outcome.upperBand)
        }
        pattern <- c(pattern,"lowerBand_","upperBand_")        
        colnames(outcome.lowerBand) <- paste0("lowerBand_",name.obs)
        colnames(outcome.upperBand) <- paste0("upperBand_",name.obs)

    }

    dataW <- data.table::data.table(
                                cbind(row = 1:n.obs, outcome,
                                      outcome.lower, outcome.upper,
                                      outcome.lowerBand,outcome.upperBand)
                         )
    dataL <- data.table::melt(dataW, id.vars = "row",
                              measure = patterns(pattern),
                              variable.name = "time.factor", value.name = gsub("_","",pattern))
    dataL$time <- times[dataL$time.factor]
    if(diag){
        dataL$row <- 1:length(name.obs)
    }else if(baseline){
        if(is.null(strata)){
            dataL$row <- 1
        }else{
            dataL$row <- strata
        }
        
    }

    ## ** add variables
    ## *** strata
    if(!is.null(strata)){ 

        if(diag || baseline){
            dataL[, strata := strata]
        }else{
            dataL[, strata := strata[row]]
        }
        
        ## if(!is.null(attr(strata,"covariates"))){
        ##     dataW[,c(names(attr(strata,"covariates"))) := attr(strata,"covariates") ]
        ## }
    }

    ## *** covariates from newdata
    if(!is.null(newdata)){ ## predicted survival/risk/average treatment effect
        dataL <- data.table::copy(merge(dataL, cbind(row = 1:NROW(newdata), newdata), by = "row"))

        ## round numeric variables
        if(any(sapply(newdata, is.numeric))){
            col.num <- names(newdata)[sapply(newdata, is.numeric)]
            dataL[, (col.num) := lapply(.SD, round, digits =  digits), .SDcols = col.num]
        }
    }
    if(!is.null(newdata) & group.by == "covariates"){
        ## keep name of the covariate for when merging all covariates together
        attr(dataL,"covariates") <- paste(names(newdata), collapse = ", ")
        dataL[, c("covariates") := interaction(.SD,sep = ", "), .SDcols = names(newdata)]
        dataL[, (attr(dataL,"covariates")) := .SD$covariates]
        ## create twice the same as data.table does not like .SDcols = "name with a ,"
    }

    ## *** status + covariates from original data
    if(!is.null(status)){ 
        dataL$status <- status        
    }

    ## ** export
    return(dataL)    
}

## * predict2plot
predict2plot <- function(dataL, name.outcome,
                         ci, band, group.by, smooth,                        
                         conf.level, alpha, xlab, ylab, xlim = NULL, ylim = NULL,
                         smoother = NULL, formula.smoother = NULL, first.derivative = FALSE, atRisk = FALSE,
                         size.estimate = 2, size.point = 3, size.ci = 1.1, col.ribbon = TRUE, size.band = 1.1, size.atRisk = 5, space.atRisk = 0.05, shape.point = c(3,18,17), n.sim = 250){

    .GRP <- .data <- NULL ## [:: for CRAN CHECK::]
    if(first.derivative && (smooth==FALSE)){
        stop("Set argument \'smooth\' to TRUE when \'first.derivative\' is TRUE. \n")
    }
    dataL <- data.table::copy(dataL)

    ## ** set at t- the value of t-1
    vec.outcome <- name.outcome
    if(ci){
        vec.outcome <- c(vec.outcome,"lowerCI","upperCI")
    }
    if(band){
        vec.outcome <- c(vec.outcome,"lowerBand","upperBand")
    }
    group.by2 <- unique(c(group.by,"row"))
    dataL[,c("timeRight") := c(.SD$time[2:.N]-1e-12,.SD$time[.N]+1e-12), by = group.by2] 

    ## ** smooth
    if(smooth>=1){
        requireNamespace("mgcv",quietly=FALSE)
        if(is.null(smoother)){
            tol <- 1e-12
            test.increasing <- all(na.omit(dataL[,diff(.SD[[name.outcome]]),by="row"][[2]])>=-tol)
            test.decreasing <- all(na.omit(dataL[,diff(.SD[[name.outcome]]),by="row"][[2]])<=tol)
            if(!requireNamespace("scam",quietly=TRUE) || (test.increasing==FALSE && test.decreasing == FALSE)){
                formula.smoother <- ~s(time)
                smoother <- mgcv::gam
            }else{ 
                smoother <- function(formula, data){
                    out <- try(do.call(scam::scam, args = list(formula = formula, data = data)))
                    if(inherits(x=out,what="try-error")){
                        out <- do.call(mgcv::gam, args = list(formula = update(formula, ".~s(time)"), data = data))
                    }
                    return(out)
                }
                if(test.increasing){
                    formula.smoother <- ~s(time, bs = "mpi")
                }else if(test.decreasing){
                    formula.smoother <- ~s(time, bs = "mpd")
                }
            }
        }
        if(length(all.vars(formula.smoother))!=1){
            stop("Argument \'formula.smoother\' must contain exactly one variable \n")
        }
        ff <- update(as.formula(paste0(name.outcome,"~.")),formula.smoother)
        if(first.derivative){
            requireNamespace("numDeriv")

            warper <- function(data){ ## data <- dataL[row==1]
                iModel <- do.call(smoother, args = list(formula = ff, data = data))
                return(numDeriv::grad(function(x){predict(iModel, newdata = data.frame(time = x), type = "response")}, x = data$time))
            }

            dataL[, c(paste0(name.outcome,".smooth")) := warper(.SD), by = "row"]
            if(ci){
                warperCI <- function(data, Sigma, n.sim){ ## data <- dataL[row==2] ; Sigma <- attr(first.derivative,"vcov")[[2]]
                    ls.deriv <- lapply(1:n.sim, function(x){
                        data2 <- data.table::copy(data)
                        data2[, c(name.outcome) := mvtnorm::rmvnorm(1, mean = data[[name.outcome]], sigma = Sigma)[1,]]                        
                        iModel <- try(do.call(smoother, args = list(formula = ff,data = data2)))
                        if(inherits(x=iModel,what="try-error")){
                            return(rep(NA, length(data$time)))
                        }else{
                            return(numDeriv::grad(function(x){predict(iModel, newdata = data.frame(time = x), type = "response")}, x = data$time))
                        }
                    })
                    M.CI <- apply(do.call(rbind,ls.deriv), 2, quantile, probs = c((1-conf.level)/2,1-(1-conf.level)/2), na.rm = TRUE)
                    return(as.data.table(t(M.CI)))
                }
                dataL[, c("lowerCI.smooth","upperCI.smooth") := warperCI(.SD, Sigma = attr(first.derivative,"vcov")[[.GRP]], n.sim = n.sim), by = "row"]
            }
            if(band){
                stop("Confidence bands are not available when argument \'first.derivative\' is TRUE \n")
            }
        }else{
            mysmoother <- function(formula, data){
                data.NNA <- data[!is.na(data[[all.vars(formula)[1]]])]
                fit <- do.call(smoother, args = list(formula = ff,data = data.NNA))
                predict(fit, newdata = data)
            }
            dataL[, c(paste0(name.outcome,".smooth")) := do.call(mysmoother, args = list(formula = ff,data = .SD)), by = "row"]
            if(ci){
                ff <- update(as.formula("lowerCI~."),formula.smoother)
                dataL[, c("lowerCI.smooth") := do.call(mysmoother, args = list(formula = ff,data = .SD)), by = "row"]
                ff <- update(as.formula("upperCI~."),formula.smoother)
                dataL[, c("upperCI.smooth") := do.call(mysmoother, args = list(formula = ff,data = .SD)), by = "row"]
            }
            if(band){
                ff <- update(as.formula("lowerBand~."),formula.smoother)
                dataL[, c("lowerBand.smooth") := do.call(mysmoother, args = list(formula = ff,data = .SD)), by = "row"]
                ff <- update(as.formula("upperBand~."),formula.smoother)
                dataL[, c("upperBand.smooth") := do.call(mysmoother, args = list(formula = ff,data = .SD)), by = "row"]
            }
        }
        
    }
    ## ** display
    labelCI <- paste0(conf.level*100,"% pointwise \n confidence interval")
    labelBand <- paste0(conf.level*100,"% simulaneous \n confidence interval \n")

    group.byL <- attr(group.by,"label")
    dataL[,c(group.byL) := as.factor(.SD[[group.byL]])]

    gg.base <- ggplot2::ggplot(data = dataL)

    ## *** confidence band
    if(band){ 
        if(smooth>0){
            if(!is.na(alpha) && col.ribbon){
                gg.base <- gg.base + ggplot2::geom_ribbon(ggplot2::aes(x = .data$time, ymin = .data$lowerBand.smooth, ymax = .data$upperBand.smooth, group = .data[[group.byL]],
                                                                       fill = .data[[group.byL]]), alpha = alpha)
            }else if(!is.na(alpha) && !col.ribbon){
                gg.base <- gg.base + ggplot2::geom_ribbon(ggplot2::aes(x = .data$time, ymin = .data$lowerBand.smooth, ymax = .data$upperBand.smooth, group = .data[[group.byL]],
                                                                       fill = labelBand), alpha = alpha)
            }else{
                gg.base <- gg.base + ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$lowerBand.smooth, group = .data[[group.byL]], color = .data[[group.byL]], linetype = labelBand),
                                                        linewidth = size.band)
                gg.base <- gg.base + ggplot2::geom_line(ggplot2::aes(x = .data$time, y = .data$upperBand.smooth, group = .data[[group.byL]], color = .data[[group.byL]], linetype = labelBand),
                                                        linewidth = size.band)
                gg.base <- gg.base + ggplot2::labs(linetype = "")
            }
        }else{
            if(!is.na(alpha) && col.ribbon){
                gg.base <- gg.base + ggplot2::geom_rect(ggplot2::aes(xmin = .data$time, xmax = .data$timeRight, ymin = .data$lowerBand, ymax = .data$upperBand,
                                                                     fill = .data[[group.byL]], group = .data[[group.byL]]), linetype = 0, alpha = alpha)
            }else if(!is.na(alpha) && !col.ribbon){
                gg.base <- gg.base + ggplot2::geom_rect(ggplot2::aes(xmin = .data$time, xmax = .data$timeRight, ymin = .data$lowerBand, ymax = .data$upperBand,
                                                                     fill = labelBand, group = .data[[group.byL]]), linetype = 0, alpha = alpha)
                gg.base <- gg.base + scale_fill_manual("", values="grey12")        
                gg.base <- gg.base + ggplot2::labs(fill = "")
            }else{
                gg.base <- gg.base + ggplot2::geom_segment(ggplot2::aes(x = .data$time, y = .data$lowerBand, xend = .data$timeRight, yend = .data$lowerBand, color = labelBand, linetype = labelBand, group = .data[[group.byL]]),
                                                           linewidth = size.band)
                gg.base <- gg.base + ggplot2::geom_segment(ggplot2::aes(x = .data$time, y = .data$upperBand, xend = .data$timeRight, yend = .data$upperBand, color = labelBand, linetype = labelBand, group = .data[[group.byL]]),
                                                           linewidth = size.band)
                gg.base <- gg.base + ggplot2::labs(color = "", linetype = "")
            }
        }
    }

    ## *** confidence interval
    if(ci){
        if(smooth>0){
            if(!is.na(alpha) && band){
                gg.base <- gg.base + ggplot2::geom_errorbar(data = dataL[!is.na(dataL$lowerCI.smooth) & !is.na(dataL$upperCI.smooth)],
                                                            ggplot2::aes(x = .data$time, ymin = .data$lowerCI.smooth, ymax = .data$upperCI.smooth, linetype = labelCI, group = .data[[group.byL]]),
                                                            width = size.ci)
                gg.base <- gg.base + ggplot2::labs(linetype = "")
            }else if(!is.na(alpha) && !band & col.ribbon){
                gg.base <- gg.base + ggplot2::geom_ribbon(data = dataL[!is.na(dataL$lowerCI.smooth) & !is.na(dataL$upperCI.smooth)],
                                                          ggplot2::aes(x = .data$time, ymin = .data$lowerCI.smooth, ymax = .data$upperCI.smooth, group = .data[[group.byL]],
                                                                       fill = .data[[group.byL]]),
                                                          alpha = alpha)
            }else if(!is.na(alpha) && !band & !col.ribbon){
                gg.base <- gg.base + ggplot2::geom_ribbon(data = dataL[!is.na(dataL$lowerCI.smooth) & !is.na(dataL$upperCI.smooth)],
                                                          ggplot2::aes(x = .data$time, ymin = .data$lowerCI.smooth, ymax = .data$upperCI.smooth, group = .data[[group.byL]],
                                                                       fill = labelBand),
                                                          alpha = alpha)
            }else{
                gg.base <- gg.base + ggplot2::geom_line(data = dataL[!is.na(dataL$lowerCI.smooth) & !is.na(dataL$upperCI.smooth)],
                                                        ggplot2::aes(x = .data$time, y = .data$lowerCI.smooth, color = .data[[group.byL]], linetype = labelCI, group = .data[[group.byL]]),
                                                        linewidth = size.ci)
                gg.base <- gg.base + ggplot2::geom_line(data = dataL[!is.na(dataL$lowerCI.smooth) & !is.na(dataL$upperCI.smooth)],
                                                        ggplot2::aes(x = .data$time, y = .data$upperCI.smooth, color = .data[[group.byL]], linetype = labelCI, group = .data[[group.byL]]),
                                                        linewidth = size.ci)
                gg.base <- gg.base + ggplot2::labs(linetype = "")
            }
        }else{
            if(!is.na(alpha) && band){
                gg.base <- gg.base + ggplot2::geom_errorbar(data = dataL[!is.na(dataL$lowerCI) & !is.na(dataL$upperCI)],
                                                            ggplot2::aes(x = .data$time, ymin = .data$lowerCI, ymax = .data$upperCI, linetype = labelCI, group = .data[[group.byL]]),
                                                            width = size.ci)
                gg.base <- gg.base + ggplot2::labs(linetype = "")
            }else if(!is.na(alpha) && !band){
                
                ## count the numbers of censored observations after each event
                dataL$following <- stats::ave(is.na(dataL$lowerCI), ## count the numbers of NA
                                              cumsum(!is.na(dataL$lowerCI)), ## between two consecutive TRUE
                                              FUN = function(x){rep(sum(x),length(x))})
                ## update the event time to the one at the latest censored observation before the next jump
                dataL$timeRightFollowing <- dataL$timeRight[(1:NROW(dataL))+dataL$following]
                ## remove lines corresponding to censored observations
                dataL.NNA <- dataL[!is.na(dataL$lowerCI) & !is.na(dataL$upperCI),.SD,.SDcols = c("time","lowerCI","upperCI",group.byL,"timeRightFollowing")]
                ## second dataset corresponding to the the latest censored observation before the next jump (just before jump)
                dataL.NNA2 <- cbind(time = dataL.NNA$timeRightFollowing, dataL.NNA[,.SD,.SDcols = c("lowerCI","upperCI",group.byL,"timeRightFollowing")])                
                ## update display                
                if(col.ribbon){
                    gg.base <- gg.base + ggplot2::geom_ribbon(data = rbind(dataL.NNA,dataL.NNA2),
                                                              ggplot2::aes(x = .data$time, ymin = .data$lowerCI, ymax = .data$upperCI, group = .data[[group.byL]],
                                                                           fill = .data[[group.byL]]),
                                                              alpha = alpha)
                }else{
                    gg.base <- gg.base + ggplot2::geom_ribbon(data = rbind(dataL.NNA,dataL.NNA2),
                                                              ggplot2::aes(x = .data$time, ymin = .data$lowerCI, ymax = .data$upperCI, group = .data[[group.byL]],
                                                                           fill = labelBand),
                                                              alpha = alpha)
                }
            }else{
                gg.base <- gg.base + ggplot2::geom_segment(data = dataL[!is.na(dataL$lowerCI) & !is.na(dataL$upperCI)],
                                                           mapping = ggplot2::aes(x = .data$time, y = .data$lowerCI, xend = .data$timeRight, yend = .data$lowerCI, color = labelCI, linetype = labelCI, group = .data[[group.byL]]),
                                                           linewidth = size.ci)
                gg.base <- gg.base + ggplot2::geom_segment(data = dataL[!is.na(dataL$lowerCI) & !is.na(dataL$upperCI)],
                                                           mapping = ggplot2::aes(x = .data$time, y = .data$upperCI, xend = .data$timeRight, yend = .data$upperCI, color = labelCI, linetype = labelCI, group = .data[[group.byL]]),
                                                           linewidth = size.ci)
                gg.base <- gg.base + ggplot2::labs(color = "", linetype = "")
            }
        }
    }

    ## *** estimate
    if(smooth>0){
        gg.base <- gg.base + ggplot2::geom_line(mapping = ggplot2::aes(x = .data$time, y = .data[[paste0(name.outcome,".smooth")]], group = .data[[group.byL]], color = .data[[group.byL]]),
                                                linewidth = size.estimate)
    }else{
        gg.base <- gg.base + ggplot2::geom_segment(data = dataL[!is.na(dataL[[name.outcome]])],
                                                   mapping = ggplot2::aes(x = .data$timeRight, y = .data[[name.outcome]], xend = .data$time, yend = .data[[name.outcome]], color = .data[[group.byL]], group = .data[[group.byL]]),
                                                   linewidth = size.estimate)

        if("status" %in% names(dataL) &  size.point>0){
            if(is.numeric(dataL$status)){ ## predictCox  (0,1)
                dataL$status <- as.character(dataL$status)
            }
            gg.base <- gg.base + ggplot2::geom_point(data = dataL[!is.na(dataL$status)],
                                                     mapping = ggplot2::aes(x = .data$time, y = .data[[name.outcome]], color = .data[[group.byL]], shape = .data$status, group = .data[[group.byL]]), size = size.point)
            if(is.factor(dataL$status)){
                gg.base <- gg.base + ggplot2::scale_shape_manual(name = "Type of event", breaks = levels(dataL$status), values = shape.point)
            }else{
                gg.base <- gg.base + ggplot2::scale_shape_manual(name = "Type of event", breaks = c(0,1), values = shape.point[1:2], labels = c("Censoring","Event"))
            }
        }else if(size.point>0){
            gg.base <- gg.base + ggplot2::geom_point(data = dataL[!is.na(dataL[[name.outcome]])],
                                                     mapping = ggplot2::aes(x = .data$time, y = .data[[name.outcome]], color = .data[[group.byL]], group = .data[[group.byL]]), size = size.point)
        }
    }

    ## *** fix legend
    if(group.by=="row"){
        gg.base <- gg.base + ggplot2::labs(color="observation") + ggplot2::theme(legend.key.height=unit(0.1,"npc"),
                                                                                 legend.key.width=unit(0.08,"npc"))
        if(!is.na(alpha) && col.ribbon && (band||ci)){
            if(length(unique(dataL$row))==1){
                gg.base <- gg.base + ggplot2::scale_fill_discrete("", labels = ifelse(band, labelBand, labelCI))
            }else{
                gg.base <- gg.base + ggplot2::labs(fill="observation")
            }
        }
        if(length(unique(dataL$row))==1){
            gg.base <- gg.base + ggplot2::scale_color_discrete(guide="none")            
        }
    }

    if((band || ci) && is.na(alpha)){
        indexTempo <- which(c(ci,band)==1)
        if(smooth == FALSE){
            levels.group.byL <- levels(dataL[[group.byL]])
            n.levels.group.byL <- length(levels.group.byL)
            gg.base <- gg.base + ggplot2::scale_color_manual("", breaks = c(c("ci","band")[indexTempo],levels.group.byL),
                                                             labels = c(c(labelCI,labelBand)[indexTempo],paste0(group.byL," ",levels.group.byL)),
                                                             values = c(c("grey","black")[indexTempo],
                                                                        grDevices::hcl(h = seq(15, 375, length = n.levels.group.byL + 1), l = 65, c = 100)[1:n.levels.group.byL]))

        }else{
            gg.base <- gg.base + ggplot2::scale_linetype_manual("", values = c("dotdash","longdash"))
        }
    }
    
    gg.base <- gg.base + ggplot2::labs(x = xlab, y = ylab)

    
    ## *** coordinates
    if(!is.null(xlim)){
        gg.base <- gg.base + ggplot2::coord_cartesian(xlim = xlim)
    }else if(name.outcome != "lp"){
        gg.base <- gg.base + ggplot2::coord_cartesian(xlim = c(0,max(dataL$timeRight)))
    }
    if(!is.null(ylim)){
        gg.base <- gg.base + ggplot2::coord_cartesian(ylim = ylim)
    }
    
    ## *** add atRisk table
    if(any(atRisk>0)){
        if(length(space.atRisk)==1){
            space.atRisk <- rep(space.atRisk,2)
        }
        ggBUILD.base <- ggplot2::ggplot_build(gg.base)
        
        ybreaks <- ggBUILD.base$layout$panel_params[[1]]$y$breaks
        if(!is.null(ylim)){
            yatRisk <- max(min(ylim), -space.atRisk[1])
        }else{
            yatRisk <- round(min(c(unlist(lapply(ggBUILD.base$data, "[[", "y")),
                                   unlist(lapply(ggBUILD.base$data, "[[", "ymin"))), na.rm = TRUE) - 0.05,2)
        }
        if(length(atRisk)==1){
            atRisk <- ggBUILD.base$layout$panel_params[[1]]$x$breaks
        }

        df.atRisk <- do.call(rbind,lapply(atRisk, function(iT){ ## iT <- 1
            if("strata" %in% names(dataL)){
                iOut <- data.frame(time = iT, y = yatRisk - space.atRisk[2] * 0:(length(levels(dataL$strata))-1), strata = levels(dataL$strata), atRisk = tapply(dataL$time>=iT,dataL$strata,sum))
            }else{
                iOut <- data.frame(time = iT, y = yatRisk, atRisk = sum(dataL$time>=iT))
            }
            return(iOut)
        }))

        gg.base <- gg.base + scale_y_continuous(breaks = unique(c(yatRisk, ybreaks)), ## in case yatRisk matches a break
                                                labels = c("at risk",ybreaks)[!duplicated(c(yatRisk, ybreaks))])
        if("strata" %in% names(dataL)){
            gg.base <- gg.base + ggplot2::geom_text(data = df.atRisk, ggplot2::aes(x = .data$time, y = .data$y, label = .data$atRisk, color = .data$strata), size = size.atRisk)
        }else{
            gg.base <- gg.base + ggplot2::geom_text(data = df.atRisk, ggplot2::aes(x = .data$time, y = .data$y, label = .data$atRisk), size = size.atRisk)
        }
        if(is.null(ylim)){
            ggBUILD.base <- ggplot2::ggplot_build(gg.base)
            gg.base <- gg.base + ggplot2::coord_cartesian(ylim = round(range(c(ggBUILD.base$layout$panel_params[[1]]$y$limits,
                                                                               unlist(lapply(ggBUILD.base$data, "[[", "y")),
                                                                               unlist(lapply(ggBUILD.base$data, "[[", "ymin")),
                                                                               unlist(lapply(ggBUILD.base$data, "[[", "yend"))), na.rm = TRUE),2))
        }
        
    }

    ## ** export
    ls.export <- list(plot = gg.base,
                      data = dataL)
    
    return(ls.export)
}



#----------------------------------------------------------------------
### autoplot.predictCox.R ends here
