### autoplot.predictCox.R --- 
##----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 17 2017 (10:06) 
## Version: 
## last-updated: Oct 16 2024 (09:24) 
##           By: Brice Ozenne
##     Update #: 1299
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
#' @name autoplot.predictCox
#'  
#' @param object Object obtained with the function \code{predictCox}.
#' @param type [character] The type of predicted value to display.
#' Choices are:
#' \code{"hazard"} the hazard function,
#' \code{"cumhazard"} the cumulative hazard function, 
#' or \code{"survival"} the survival function.
#' @param ci [logical] If \code{TRUE} display the confidence intervals for the predictions.
#' @param band [logical] If \code{TRUE} display the confidence bands for the predictions.
#' @param plot [logical] Should the graphic be plotted.
#' @param digits [integer] Number of decimal places when displaying the values of the covariates in the caption.
#' @param alpha [numeric, 0-1] Transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param ylab [character] Label for the y axis.
#' @param smooth [logical] Should a smooth version of the risk function be plotted instead of a simple function?
#' @param first.derivative [logical] If \code{TRUE}, display the first derivative over time of the risks/risk differences/risk ratios.
#' (confidence intervals are obtained via simulation).
#' @param group.by [character] The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}.
#' @param reduce.data [logical] If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param ... Additional parameters to cutomize the display.
#'
#' @return Invisible. A list containing:
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
#' autoplot(e.basehaz, type = "cumhazard")
#' \dontrun{
#' autoplot(e.basehaz, type = "cumhazard", size.point = 0) ## without points
#' autoplot(e.basehaz, type = "cumhazard", smooth = TRUE)
#' autoplot(e.basehaz, type = "cumhazard", smooth = TRUE, first.derivative = TRUE)
#' }
#' 
#' ## display baseline hazard with type of event
#' \dontrun{
#' e.basehaz <- predictCox(m.cox, keep.newdata = TRUE)
#' autoplot(e.basehaz, type = "cumhazard")
#' autoplot(e.basehaz, type = "cumhazard", shape.point = c(3,NA))
#' }
#'
#' ## display predicted survival
#' \dontrun{
#' pred.cox <- predictCox(m.cox, newdata = d[1:2,],
#'   times = seqTau, type = "survival", keep.newdata = TRUE)
#' autoplot(pred.cox)
#' autoplot(pred.cox, smooth = TRUE)
#' autoplot(pred.cox, group.by = "covariates")
#' autoplot(pred.cox, group.by = "covariates", reduce.data = TRUE)
#' autoplot(pred.cox, group.by = "X1", reduce.data = TRUE)
#' }
#' 
#' ## predictions with confidence interval/bands
#' \dontrun{
#' pred.cox <- predictCox(m.cox, newdata = d[1:2,,drop=FALSE],
#'   times = seqTau, type = "survival", band = TRUE, se = TRUE, keep.newdata = TRUE)
#' res <- autoplot(pred.cox, ci = TRUE, band = TRUE, plot = FALSE)
#' res$plot + facet_wrap(~row)
#' res2 <- autoplot(pred.cox, ci = TRUE, band = TRUE, alpha = 0.1, plot = FALSE)
#' res2$plot + facet_wrap(~row)
#' }
#' 
#' #### Stratified Cox model ####
#' \dontrun{
#' m.cox.strata <- coxph(Surv(time,event)~ strata(X1) + strata(X2) + X3 + X4,
#'                       data = d, x = TRUE, y = TRUE)
#'
#' ## baseline hazard
#' pred.baseline <- predictCox(m.cox.strata, keep.newdata = TRUE, type = "survival")
#' res <- autoplot(pred.baseline)
#' res$plot + facet_wrap(~strata, labeller = label_both)
#'
#' ## predictions
#' pred.cox.strata <- predictCox(m.cox.strata, newdata = d[1:3,,drop=FALSE],
#'                               time = seqTau, keep.newdata = TRUE, se = TRUE)
#'
#' res2 <- autoplot(pred.cox.strata, type = "survival", group.by = "strata", plot = FALSE)
#' res2$plot + facet_wrap(~strata, labeller = label_both) + theme(legend.position="bottom")
#' 
#' ## smooth version
#' autoplot(pred.cox.strata, type = "survival", group.by = "strata", smooth = TRUE, ci = FALSE)
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
#' autoplot(pred.spline, group.by = "X6")
#' autoplot(pred.spline, group.by = "X6", alpha = 0.5)
#' 
#' grid2 <- data.frame(X1 = factor(1,0:1), X2 = factor(0,0:1),
#'                     X6 = seq(min(d$X6),max(d$X6), length.out = 100))
#' pred.spline <- predictCox(m.cox.spline, newdata = rbind(grid,grid2), keep.newdata = TRUE,
#'                           se = TRUE, band = TRUE, centered = TRUE, type = "lp")
#' autoplot(pred.spline, group.by = c("X6","X1"), alpha = 0.5, plot = FALSE)$plot + facet_wrap(~X1)
#' }

## * autoplot.predictCox (code)
#' @rdname autoplot.predictCox
#' @method autoplot predictCox
#' @export
autoplot.predictCox <- function(object,
                                type = NULL,
                                ci = object$se,
                                band = object$band,
                                plot = TRUE,
                                smooth = NULL,
                                digits = 2,
                                alpha = NA,
                                group.by = "row",
                                reduce.data = FALSE,
                                ylab = NULL,
                                first.derivative = FALSE,
                                 ...){
  
    ## initialize and check    
    possibleType <- c("cumhazard","survival","lp")
    possibleType <- possibleType[possibleType %in% names(object)]
    nVar.lp <- length(object$var.lp)
        
    if(is.null(type)){
        if(length(possibleType) == 1){
            type <- possibleType
        }else{
            stop("argument \'type\' must be specified to choose between ",paste(possibleType, collapse = " "),"\n")
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
                           "cumhazard" = if(object$baseline && nVar.lp>0){"instantaneous baseline hazard"}else{"instantaneous hazard"},
                           "survival" = if(object$baseline && nVar.lp>0){"derivative of the baseline survival"}else{"derivative of the survival"})
        }else{
            ylab <- switch(type,
                           "lp" = "linear predictor",
                           "cumhazard" = if(object$baseline && nVar.lp>0){"cumulative baseline hazard"}else{"cumulative hazard"},
                           "survival" = if(object$baseline && nVar.lp>0){"baseline survival"}else{"survival"})
        }
    }

    if(type=="lp"){
        group.by <- match.arg(group.by, object$var.lp, several.ok = TRUE)
    }else{
        if(length(group.by)>1){
            stop("Argument \'group.by\' must have length 1.\n")
        }
        group.by <- match.arg(group.by, c("row","covariates","strata",object$var.lp,object$var.strata))
    }
    

    if(group.by[[1]] == "covariates" && ("newdata" %in% names(object)) == FALSE){
        stop("argument \'group.by\' cannot be \"covariates\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling the predictCox function \n")
    }

    if(group.by[[1]] %in% c("strata") && ("strata" %in% names(object) == FALSE)){
        stop("argument \'group.by\' cannot be \"strata\" when strata is missing in the object \n",
             "set argment \'keep.strata\' to TRUE when calling the predictCox function \n")
    }
    if(group.by[[1]] %in% c(object$var.lp,object$var.strata) && ("newdata" %in% names(object) == FALSE)){
        stop("argument \'group.by\' cannot be \"",group.by[[1]],"\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling the predictCox function \n")
    }
    if(!is.null(object$newdata$strata) && !is.null(object$var.strata)){
        if(any(object$var.strata %in% colnames(object$newdata))){
            stop("cannot handle covariate \"",paste(object$var.strata[object$var.strata %in% colnames(object$newdata)],collapse ="\" \""),"\" as this name is used internally.\n",
                 "consider refit the model using another name for the covariate. \n")
        }
        splitStrata <- do.call(rbind,strsplit(split = ",", x = as.character(object$newdata$strata), fixed = TRUE))
        colnames(splitStrata) <- object$infoVar$stratavars.original
        object$newdata <- cbind(object$newdata,splitStrata)
    }
        
    if(ci[[1]]==TRUE && (object$se[[1]]==FALSE || is.null(object$conf.level))){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set arguments \'se\' and \'confint\' to TRUE when calling the predictCox function \n")
    }
    if(band[[1]] && (object$band[[1]]==FALSE  || is.null(object$conf.level))){
        stop("argument \'band\' cannot be TRUE when the quantiles for the confidence bands have not been computed \n",
             "set arguments \'band\' and \'confint\' to TRUE when calling the predictCox function \n")
    }
    if(object$nTimes!= 0 && any(rank(object$times) != 1:length(object$times))){
        stop("Invalid object. The prediction times must be strictly increasing \n")
    }
    ## dots <- list(...)
    ## if(length(dots)>0){
    ##     txt <- names(dots)
    ##     txt.s <- if(length(txt)>1){"s"}else{""}
    ##     stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    ## }

    ## ** reshape data
    if(type == "lp"){
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

        if(!is.matrix(object[[type]])){
        
            ## baseline hazard/survival
            if(is.null(object[["strata"]])){
                object[[type]] <- rbind(object[[type]])
                if(object$nTimes==0){
                    if(0 %in% object$times == FALSE){
                        if(type=="cumhazard"){
                            object[[type]] <- cbind(0,object[[type]])
                        }else if(type=="survival"){
                            object[[type]] <- cbind(1,object[[type]])
                        }                   
                        object$times <- c(0,object$times)
                        if(!is.null(object$newdata)){
                            object$newdata <- rbind(data.table(start = 0, stop = 0, status = NA, strata = 1, strata.num = 0, eXb = NA, statusM1 = NA, XXXindexXXX = NA),
                                                    object$newdata)
                        }
                    }
                    if(object$lastEventTime %in% object$times == FALSE){
                        object[[type]] <- cbind(object[[type]],object[[type]][length(object[[type]])])
                        object$times <- c(object$times,pmin(object$lastEventTime,max(object$times)+1e-12))
                        if(!is.null(object$newdata)){
                            object$newdata <- rbind(data.table(start = 0, stop = object$lastEventTime, status = 0, strata = 1, strata.num = 0, eXb = NA, statusM1 = NA, XXXindexXXX = NA),
                                                    object$newdata)
                        }
                    }
                }

            }else{
                index.unique <- !duplicated(object$strata)
                strata <- object$strata[index.unique]
                if(!is.null(attr(object$strata,"covariates"))){
                    attr(strata,"covariates") <- attr(object$strata,"covariates")[index.unique]
                }
                n.strata <- length(strata)
                time <- unique(sort(object[["times"]])) 
                n.time <- length(time)
                type.tempo <- matrix(NA, nrow = n.strata, ncol = n.time)

            init <- switch(type,
                           "cumhazard" = 0,
                           "survival" = 1)

                for(iStrata in 1:n.strata){ ## iStrata <- 1
                    index.strata <- which(object[["strata"]]==strata[iStrata])
                    type.tempo[iStrata,]  <- stats::approx(x = object[["times"]][index.strata],
                                                           y = object[[type]][index.strata],
                                                           yleft = init,
                                                           yright = NA,
                                                           xout = time,
                                                           method = "constant")$y
                
            }
            object[[type]] <- type.tempo
            object[["strata"]] <- strata
            object[["times"]] <- time
        }

            newdata <- NULL
            if(object$nTimes==0){
                status <- object$newdata
            }else{
                status <- NULL
            }
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
                              status = status,
                              strata = object$strata,
                              times = object$times,
                              name.outcome = type,
                              group.by = group.by,
                              digits = digits
                              )
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
                           ...
                           )
  
    if(plot){
        print(gg.res$plot)
    }
  
    return(invisible(gg.res))
}

## * predict2melt
predict2melt <- function(outcome, name.outcome,
                         ci, outcome.lower, outcome.upper,
                         band, outcome.lowerBand, outcome.upperBand,
                         newdata, status, strata, times, 
                         group.by, digits){

    patterns <- NULL ## [:CRANtest:] data.table
    n.time <- NCOL(outcome)
    if(!is.null(time)){
        time.names <- times 
    }else{
        time.names <- 1:n.time
    }
    colnames(outcome) <- paste0(name.outcome,"_",time.names)

    ## merge outcome with CI and band ####
    pattern <- paste0(name.outcome,"_")
    if(!is.null(status)){
        Ustrata <- unique(status$strata)
        M.status <- matrix(as.numeric(NA), nrow = NROW(outcome), ncol = NCOL(outcome),
                           dimnames = list(NULL, paste0("status_",time.names)))
        status[,c("index") := match(.SD$stop, times)]

        for(iS in 1:length(Ustrata)){ ## iS <- 1
            iIndex <- which(status$strata == Ustrata[iS])
            iStatus <- status[iIndex, list("nevent" = sum(.SD$status), "index" = unique(.SD$index)),by="stop"]
            M.status[iS, iStatus$index] <- (iStatus$nevent>0)
        }
        ## M.status[1,times[]] <- status[strata==Ustrata[1]]
        ## times
        outcome <- cbind(outcome, M.status)
        pattern <- c(pattern,"status")
    }
    if(ci){
        pattern <- c(pattern,"lowerCI_","upperCI_")
    
        colnames(outcome.lower) <- paste0("lowerCI_",time.names)
        colnames(outcome.upper) <- paste0("upperCI_",time.names)
    }
    if(band){
        pattern <- c(pattern,"lowerBand_","upperBand_")
        
        colnames(outcome.lowerBand) <- paste0("lowerBand_",time.names)
        colnames(outcome.upperBand) <- paste0("upperBand_",time.names)

    }
    outcome <- data.table::data.table(
                               cbind(outcome,
                                     outcome.lower, outcome.upper,
                                     outcome.lowerBand,outcome.upperBand)
                           )

    ## merge with covariates ####
    outcome[, row := 1:.N]
    keep.col <- "row"
    if(!is.null(newdata)){ ## predicted survival/risk/average treatment effect
        cov.names <- names(newdata)

        ## used individually
        newdata <- newdata[, (cov.names) := lapply(cov.names,function(col){
            if (is.numeric(.SD[[col]]))
                round(.SD[[col]],digits) else .SD[[col]]})]
        outcome <- cbind(outcome,newdata)

        ## keep name of the covariate for when merging all covariates together
        newdata <- newdata[, (cov.names) := lapply(cov.names,function(col){
            paste0(col,"=",.SD[[col]])})]
        outcome[, c("covariates") := interaction(newdata,sep = " ")]
        
        keep.col <- c(keep.col,if(!is.null(cov.names)){"covariates"},cov.names)
    }

    if(!is.null(strata)){ ##  baseline survival/risk
        outcome[, strata := strata]
        if(!is.null(attr(strata,"covariates"))){
            outcome[,c(names(attr(strata,"covariates"))) := attr(strata,"covariates") ]
        }
        keep.col <- c(keep.col,"strata",names(attr(strata,"covariates")))
    }

    ## reshape to long format ####
    dataL <- data.table::melt(outcome, id.vars = keep.col,
                              measure = patterns(pattern),
                              variable.name = "time", value.name = gsub("_","",pattern))

    dataL[, time := as.numeric(as.character(factor(time, labels = time.names)))]
    dataL <- dataL[!is.na(dataL[[name.outcome]])]
    return(dataL)    
}

## * predict2plot
predict2plot <- function(dataL, name.outcome,
                         ci, band, group.by, smooth,                        
                         conf.level, alpha, xlab, ylab,
                         smoother = NULL, formula.smoother = NULL, first.derivative = FALSE,
                         size.estimate = 1.5, size.point = 3, size.ci = 1.1, size.band = 1.1, shape.point = c(3,18), n.sim = 250){
    
    .GRP <- NULL ## [:: for CRAN CHECK::]
    if(first.derivative && (smooth==FALSE)){
        stop("Set argument \'smooth\' to TRUE when \'first.derivative\' is TRUE. \n")
    }
    dataL <- data.table::copy(dataL)

    ## set at t- the value of t-1
    vec.outcome <- name.outcome
    if(ci){
        vec.outcome <- c(vec.outcome,"lowerCI","upperCI")
    }
    if(band){
        vec.outcome <- c(vec.outcome,"lowerBand","upperBand")
    }
    group.by2 <- unique(c(group.by,"row"))
    dataL[,c("timeRight") := c(.SD$time[2:.N]-1e-12,.SD$time[.N]+1e-12), by = group.by2] 
    dataL[,c(group.by) := as.factor(.SD[[group.by]])]
    
    ## smooth ####
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
            dataL[, c(paste0(name.outcome,".smooth")) := do.call(smoother, args = list(formula = ff,data = .SD))$fitted, by = "row"]
            if(ci){
                ff <- update(as.formula("lowerCI~."),formula.smoother)
                dataL[, c("lowerCI.smooth") := do.call(smoother, args = list(formula = ff,data = .SD))$fitted, by = "row"]
                ff <- update(as.formula("upperCI~."),formula.smoother)
                dataL[, c("upperCI.smooth") := do.call(smoother, args = list(formula = ff,data = .SD))$fitted, by = "row"]
            }
            if(band){
                ff <- update(as.formula("lowerBand~."),formula.smoother)
                dataL[, c("lowerBand.smooth") := do.call(smoother, args = list(formula = ff,data = .SD))$fitted, by = "row"]
                ff <- update(as.formula("upperBand~."),formula.smoother)
                dataL[, c("upperBand.smooth") := do.call(smoother, args = list(formula = ff,data = .SD))$fitted, by = "row"]
            }
        }
        
    }

    ## display ####
    labelCI <- paste0(conf.level*100,"% pointwise \n confidence interval")
    labelBand <- paste0(conf.level*100,"% simulaneous \n confidence interval \n")

    gg.base <- ggplot2::ggplot(data = dataL, mapping = ggplot2::aes(group = row))
    if(band){ ## confidence band
        if(smooth>0){
            if(!is.na(alpha)){
                gg.base <- gg.base + ggplot2::geom_ribbon(eval(parse(text = paste0(
                                                                         "ggplot2::aes(x = time, ymin = lowerBand.smooth, ymax = upperBand.smooth, group = ",group.by,")"))),
                                                          alpha = alpha)
            }else{
                gg.base <- gg.base + ggplot2::geom_line(eval(parse(text = paste0(
                                                                       "ggplot2::aes(x = time, y = lowerBand.smooth, group = ",group.by,", color = ",group.by,", linetype = \"band\")"))),
                                                        linewidth = size.band)
                gg.base <- gg.base + ggplot2::geom_line(eval(parse(text = paste0(
                                                                       "ggplot2::aes(x = time, y = upperBand.smooth, group = ",group.by, ", color = ",group.by,", linetype = \"band\")"))),
                                                        linewidth = size.band)
            }
        }else{
            if(!is.na(alpha)){
                gg.base <- gg.base + ggplot2::geom_rect(ggplot2::aes_string(xmin = "time", xmax = "timeRight", ymin = "lowerBand", ymax = "upperBand",
                                                                            fill = "labelBand"), linetype = 0, alpha = alpha)
                gg.base <- gg.base + scale_fill_manual("", values="grey12")        
            }else{
                gg.base <- gg.base + ggplot2::geom_segment(ggplot2::aes_string(x = "time", y = "lowerBand", xend = "timeRight", yend = "lowerBand", color = "\"band\""),
                                                           linewidth = size.band)
                gg.base <- gg.base + ggplot2::geom_segment(ggplot2::aes_string(x = "time", y = "upperBand", xend = "timeRight", yend = "upperBand", color = "\"band\""),
                                                           linewidth = size.band)
            }
        }
    }
    if(ci){ ## confidence interval
        if(smooth>0){
            if(!is.na(alpha)){
                gg.base <- gg.base + ggplot2::geom_errorbar(ggplot2::aes_string(x = "time", ymin = "lowerCI.smooth", ymax = "upperCI.smooth", linetype = "labelCI"),
                                                            width = size.ci)
                gg.base <- gg.base + ggplot2::scale_linetype_manual("",values=setNames(1,labelCI))
            }else{
                gg.base <- gg.base + ggplot2::geom_line(eval(parse(text = paste0(
                                                                       "ggplot2::aes(x = time, y = lowerCI.smooth, group = ",group.by,", color = ",group.by,", linetype = \"ci\")"))),
                                                        linewidth = size.ci)
                gg.base <- gg.base + ggplot2::geom_line(eval(parse(text = paste0(
                                                                       "ggplot2::aes(x = time, y = upperCI.smooth, group = ",group.by,", color = ",group.by,", linetype = \"ci\")"))),
                                                        linewidth = size.ci)

            }
        }else{
            if(!is.na(alpha)){
                gg.base <- gg.base + ggplot2::geom_errorbar(ggplot2::aes_string(x = "time", ymin = "lowerCI", ymax = "upperCI", linetype = "labelCI"),
                                                            width = size.ci)
                gg.base <- gg.base + ggplot2::scale_linetype_manual("",values=setNames(1,labelCI))

            }else{
                gg.base <- gg.base + ggplot2::geom_segment(ggplot2::aes_string(x = "time", y = "lowerCI", xend = "timeRight", yend = "lowerCI", color = "\"ci\""),
                                                           linewidth = size.ci)
                gg.base <- gg.base + ggplot2::geom_segment(ggplot2::aes_string(x = "time", y = "upperCI", xend = "timeRight", yend = "upperCI", color = "\"ci\""),
                                                           linewidth = size.ci)
            }
        }
    }
    ## estimate
    if(smooth>0){
        gg.base <- gg.base + ggplot2::geom_line(mapping = ggplot2::aes_string(x = "time", y = paste0(name.outcome,".smooth"), group = group.by, color = group.by),
                                                linewidth = size.estimate)
    }else{
        gg.base <- gg.base + ggplot2::geom_segment(mapping = ggplot2::aes_string(x = "timeRight", y = name.outcome, xend = "time", yend = name.outcome, color = group.by),
                                                   linewidth = size.estimate)
        if("status" %in% names(dataL)){
            dataL$status <- as.character(dataL$status)
            gg.base <- gg.base + ggplot2::geom_point(data = na.omit(dataL),
                                                     mapping = ggplot2::aes_string(x = "time", y = name.outcome, color = group.by, shape = "status"), size = size.point)
            gg.base <- gg.base + ggplot2::scale_shape_manual(breaks = c(0,1), values = shape.point, labels = c("censoring","event"))
        }else{
            gg.base <- gg.base + ggplot2::geom_point(data = dataL,
                                                     mapping = ggplot2::aes_string(x = "time", y = name.outcome, color = group.by), size = size.point)
        }
    }
    
    if(group.by=="row"){
        gg.base <- gg.base + ggplot2::labs(color="observation") + ggplot2::theme(legend.key.height=unit(0.1,"npc"),
                                                                                 legend.key.width=unit(0.08,"npc"))
        
        # display only integer values
        uniqueObs <- unique(dataL$row)

        if(length(uniqueObs)==1){
            gg.base <- gg.base + ggplot2::scale_color_discrete(guide="none")
        }
    }

    if(is.na(alpha)[[1]] && (band[[1]] || ci[[1]])){
        indexTempo <- which(c(ci,band)==1)
        if(smooth == FALSE){
            levels.group.by <- levels(dataL[[group.by]])
            n.levels.group.by <- length(levels.group.by)
            gg.base <- gg.base + ggplot2::scale_color_manual("", breaks = c(c("ci","band")[indexTempo],levels.group.by),
                                                             labels = c(c(labelCI,labelBand)[indexTempo],paste0(group.by," ",levels.group.by)),
                                                             values = c(c("grey","black")[indexTempo],
                                                                        grDevices::hcl(h = seq(15, 375, length = n.levels.group.by + 1), l = 65, c = 100)[1:n.levels.group.by]))
        }else{
            gg.base <- gg.base + ggplot2::scale_linetype_manual("", breaks = c("ci","band")[indexTempo],
                                                                labels = c(labelCI,labelBand)[indexTempo],
                                                                values = c("dotdash","longdash")[indexTempo])
        }
    }else if(ci[[1]] && band[[1]]){
        gg.base <- gg.base + ggplot2::guides(linetype = ggplot2::guide_legend(order = 1),
                                             fill = ggplot2::guide_legend(order = 2),
                                             group = ggplot2::guide_legend(order = 3)
                                             )
    }
    gg.base <- gg.base + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    if(name.outcome != "lp"){
        gg.base <- gg.base + ggplot2::coord_cartesian(xlim = c(0,max(dataL$timeRight)))
    }
    ## export
    ls.export <- list(plot = gg.base,
                      data = dataL)
    
    return(ls.export)
}



#----------------------------------------------------------------------
### autoplot.predictCox.R ends here
