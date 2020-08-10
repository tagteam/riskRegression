### autoplot.predictCox.R --- 
##----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 17 2017 (10:06) 
## Version: 
## last-updated: aug 10 2020 (12:05) 
##           By: Brice Ozenne
##     Update #: 604
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

                                        # {{{ autoplot.predictCox
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
#' @param group.by [character] The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}.
#' @param reduce.data [logical] If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param plot [logical] Should the graphic be plotted.
#' @param digits [integer] Number of decimal places.
#' @param ylab [character] Label for the y axis.
#' @param alpha [numeric, 0-1] Transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param ... Not used. Only for compatibility with the plot method.

## * autoplot.predictCox (examples)
#' @examples
#' library(survival)
#' library(ggplot2)
#'
#' #### simulate data ####
#' set.seed(10)
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' #### Cox model ####
#' m.cox <- coxph(Surv(time,event)~ X1 + X2 + X3,
#'                 data = d, x = TRUE, y = TRUE)
#'
#' ## display baseline hazard
#' e.basehaz <- predictCox(m.cox)
#' autoplot(e.basehaz, type = "cumhazard")
#'
#' ## display baseline hazard with type of event 
#' e.basehaz <- predictCox(m.cox, keep.newdata = TRUE)
#' autoplot(e.basehaz, type = "cumhazard")
#'
#' ## display predicted survival
#' pred.cox <- predictCox(m.cox, newdata = d[1:4,],
#'   times = 1:5, type = "survival", keep.newdata = TRUE)
#' autoplot(pred.cox)
#' autoplot(pred.cox, group.by = "covariates")
#' autoplot(pred.cox, group.by = "covariates", reduce.data = TRUE)
#' 
#' ## predictions with confidence interval/bands
#' pred.cox <- predictCox(m.cox, newdata = d[1,,drop=FALSE],
#'   times = 1:5, type = "survival", band = TRUE, se = TRUE, keep.newdata = TRUE)
#' autoplot(pred.cox, ci = TRUE, band = TRUE)
#' autoplot(pred.cox, ci = TRUE, band = TRUE, alpha = 0.1)
#'
#' #### Stratified Cox model ####
#' m.cox.strata <- coxph(Surv(time,event)~ strata(X1) + strata(X2) + X3 + X4,
#'                       data = d, x = TRUE, y = TRUE)
#'
#' ## baseline hazard
#' pred.baseline <- predictCox(m.cox.strata, keep.newdata = TRUE, type = "survival")
#' res <- autoplot(pred.baseline)
#' res$plot + facet_wrap(~strata, labeller = label_both)
#' 
#' ## predictions
#' pred.cox.strata <- predictCox(m.cox.strata, newdata = d[1:5,,drop=FALSE],
#'                               time = 1:5, keep.newdata = TRUE)
#'
#' autoplot(pred.cox.strata, type = "survival", group.by = "strata")

## * autoplot.predictCox (code)
#' @rdname autoplot.predictCox
#' @method autoplot predictCox
#' @export
autoplot.predictCox <- function(object,
                                type = NULL,
                                ci = FALSE,
                                band = FALSE,
                                group.by = "row",
                                reduce.data = FALSE,
                                plot = TRUE,
                                ylab = NULL,
                                digits = 2, alpha = NA, ...){
  
    ## initialize and check    
    possibleType <- c("cumhazard","survival")
    possibleType <- possibleType[possibleType %in% names(object)]

    if(is.null(type)){
        if(length(possibleType) == 1){
            type <- possibleType
        }else{
            stop("argument \'type\' must be specified to choose between ",paste(possibleType, collapse = " "),"\n")
        }
    }else{
        type <- match.arg(type, possibleType)  
    }
    if(is.null(ylab)){
        ylab <- switch(type,
                        "cumhazard" = "cumulative hazard",
                        "survival" = "survival")
    }

    group.by <- match.arg(group.by, c("row","covariates","strata"))

    if((group.by[[1]] == "covariates") && ("newdata" %in% names(object)) == FALSE){
        stop("argument \'group.by\' cannot be \"covariates\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling the predictCox function \n")
    }

    if(group.by[[1]] == "strata" && ("strata" %in% names(object) == FALSE)){
        stop("argument \'group.by\' cannot be \"strata\" when strata is missing in the object \n",
             "set argment \'keep.strata\' to TRUE when calling the predictCox function \n")
    }
  
    if(ci[[1]]==TRUE && (object$se[[1]]==FALSE || is.null(object$conf.level))){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set arguments \'se\' and \'confint\' to TRUE when calling the predictCox function \n")
    }
    if(band[[1]] && (object$band[[1]]==FALSE  || is.null(object$conf.level))){
        stop("argument \'band\' cannot be TRUE when the quantiles for the confidence bands have not been computed \n",
             "set arguments \'band\' and \'confint\' to TRUE when calling the predictCox function \n")
    }

    dots <- list(...)
    if(length(dots)>0){
        txt <- names(dots)
        txt.s <- if(length(txt)>1){"s"}else{""}
        stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    }

    ## reshape data
    if(!is.matrix(object[[type]])){
        if(is.null(object[["strata"]])){
             object[[type]] <- rbind(object[[type]])
        }else{
            strata <- unique(object[["strata"]])
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
            group.by <- "strata"
        }
        newdata <- NULL
        status <- object$status
        
    }else{
        newdata <- data.table::copy(object$newdata) ## can be NULL
        if(!is.null(newdata) && reduce.data[[1]]==TRUE){
            test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
            if(any(test)){
                newdata[, (names(test)[test]):=NULL]
            }        
        }
        status <- NULL
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

    ## display

    gg.res <- predict2plot(dataL = dataL,
                           name.outcome = type,
                           ci = ci,
                           band = band,
                           group.by = group.by,
                           conf.level = object$conf.level,
                           alpha = alpha,
                           xlab = if(is.null(object$infoVar)){"time"}else{object$infoVar$time},
                           ylab = ylab
                           )
  
  if(plot){
    print(gg.res$plot)
  }
  
    return(invisible(gg.res))
}
                                        # }}}

                                        # {{{ predict2melt
## * predict2melt
predict2melt <- function(outcome, name.outcome,
                         ci, outcome.lower, outcome.upper,
                         band, outcome.lowerBand, outcome.upperBand,
                         newdata, status, strata, times, group.by, digits){

    patterns <- NULL ## [:CRANtest:] data.table

    n.time <- NCOL(outcome)
    if(!is.null(time)){
        time.names <- times 
    }else{
        time.names <- 1:n.time
    }    
    colnames(outcome) <- paste0(name.outcome,"_",time.names)
    keep.cols <- unique(c("time",name.outcome,"row",group.by))

    ## add initial values ####
    first.dt <- switch(name.outcome,
                       "cumhazard" = data.table::data.table(time = 0, cumhazard = 0, origin = TRUE),
                       "survival" = data.table::data.table(time = 0, survival = 1, origin = TRUE),
                       "absRisk" = data.table::data.table(time = 0, absRisk = 0, origin = TRUE))
    
    ## merge outcome with CI and band ####
    pattern <- paste0(name.outcome,"_")
    if(!is.null(status)){
        Ustrata <- unique(status$strata)

        M.status <- matrix(as.numeric(NA), nrow = NROW(outcome), ncol = NCOL(outcome),
                           dimnames = list(NULL, paste0("status_",time.names)))
        status[,c("index") := match(.SD$time, times)]

        for(iS in 1:length(Ustrata)){ ## iS <- 1
            iStatus <- status[status$strata == Ustrata[iS]]
            M.status[iS, iStatus$index] <- (iStatus$nevent>0)
        }
        ## M.status[1,times[]] <- status[strata==Ustrata[1]]
        ## times
        outcome <- cbind(outcome, M.status)
        pattern <- c(pattern,"status")
        first.dt[, status := as.character(NA)]
    }

    if(ci){
        pattern <- c(pattern,"lowerCI_","upperCI_")
    
        colnames(outcome.lower) <- paste0("lowerCI_",time.names)
        colnames(outcome.upper) <- paste0("upperCI_",time.names)
        first.dt[, c("lowerCI") := .SD[[1]], .SDcols = name.outcome]
        first.dt[, c("upperCI") := .SD[[1]], .SDcols = name.outcome]
    }
    if(band){
        pattern <- c(pattern,"lowerBand_","upperBand_")
        keep.cols <- c(keep.cols,"lowerBand","upperBand")
        
        colnames(outcome.lowerBand) <- paste0("lowerBand_",time.names)
        colnames(outcome.upperBand) <- paste0("upperBand_",time.names)

        first.dt[, c("lowerBand") :=  .SD[[1]], .SDcols = name.outcome]
        first.dt[, c("upperBand") :=  .SD[[1]], .SDcols = name.outcome]
    }

    outcome <- data.table::data.table(
                               cbind(outcome,
                                     outcome.lower, outcome.upper,
                                     outcome.lowerBand,outcome.upperBand)
                           )

    ## merge with convariates ####
    outcome[, row := 1:.N]
    if(group.by == "covariates"){
        cov.names <- names(newdata)
        newdata <- newdata[, (cov.names) := lapply(cov.names,function(col){
            if (is.numeric(.SD[[col]]))
                paste0(col,"=",round(.SD[[col]],digits)) else paste0(col,"=",.SD[[col]])})]
        outcome[, ("covariates") := interaction(newdata,sep = " ")]
    }else if(group.by == "strata"){
        outcome[, strata := strata]
    }
    
    ## reshape to long format ####
    dataL <- melt(outcome, id.vars = union("row",group.by),
                  measure = patterns(pattern),
                  variable.name = "time", value.name = gsub("_","",pattern))

    dataL[, time := as.numeric(as.character(factor(time, labels = time.names)))]
    dataL[, c("origin") := FALSE]
    dataL <- dataL[!is.na(dataL[[name.outcome]])]
    dataL <- dataL[, rbind(first.dt,.SD), by = c(union("row",group.by))]
    return(dataL)    
}

                                        # }}}
                                        # {{{ predict2plot
## * predict2plot
predict2plot <- function(dataL, name.outcome,
                         ci, band, group.by,                         
                         conf.level, alpha, xlab, ylab){
                                        # for CRAN tests
    original <- lowerCI <- upperCI <- lowerBand <- upperBand <- timeLeft <- NULL
    ## duplicate observations to obtain step curves ####
    keep.cols <- unique(c("time",name.outcome,"row",group.by,"original"))
    if(ci){
        keep.cols <- c(keep.cols,"lowerCI","upperCI")
    }
    if(band){
        keep.cols <- c(keep.cols,"lowerBand","upperBand")
    }

    ## set at t- the value of t-1
    dtTempo <- copy(dataL)

    vec.outcome <- name.outcome
    if(ci){
        vec.outcome <- c(vec.outcome,"lowerCI","upperCI")
    }
    if(band){
        vec.outcome <- c(vec.outcome,"lowerBand","upperBand")
    }
    dataL[,c("timeLeft") := .SD[c(1,1:(.N-1))], .SDcols = "time", by = "row"] 
    
    ## display ####
    labelCI <- paste0(conf.level*100,"% confidence \n interval")
    labelBand <- paste0(conf.level*100,"% confidence \n band")

    gg.base <- ggplot(data = dataL, mapping = aes_string(group = "row"))

    if(band){ ## confidence band
        if(!is.na(alpha)){
            gg.base <- gg.base + ggplot2::geom_rect(aes(xmin = timeLeft, xmax = time, ymin = lowerBand, ymax = upperBand,
                                                        fill = labelBand), linetype = 0, alpha = alpha)
            gg.base <- gg.base + scale_fill_manual("", values="grey12")        
        }else{
            gg.base <- gg.base + ggplot2::geom_segment(aes(x = timeLeft, y = lowerBand, xend = time, yend = lowerBand, linetype = "band"),
                                                       size = 1.2, color = "black")
            gg.base <- gg.base + ggplot2::geom_segment(aes(x = timeLeft, y = upperBand, xend = time, yend = upperBand, linetype = "band"),
                                                       size = 1.2, color = "black")
        }
    }
    
    if(ci){ ## confidence interval
        if(!is.na(alpha)){
            gg.base <- gg.base + ggplot2::geom_errorbar(aes(x = time, ymin = lowerCI, ymax = upperCI, linetype = labelCI))
            gg.base <- gg.base + scale_linetype_manual("",values=setNames(1,labelCI))

        }else{
            gg.base <- gg.base + ggplot2::geom_segment(aes(x = timeLeft, y = lowerCI, xend = time, yend = lowerCI, linetype = "ci"),
                                                       size = 1.2, color = "black")
            gg.base <- gg.base + ggplot2::geom_segment(aes(x = timeLeft, y = upperCI, xend = time, yend = upperCI, linetype = "ci"),
                                                       size = 1.2, color = "black")
        }
    }

    ## estimate
    gg.base <- gg.base + ggplot2::geom_segment(mapping = aes_string(x = "timeLeft", y = name.outcome, xend = "time", yend = name.outcome, color = group.by), size = 1.5)
    if("status" %in% names(dataL)){
        gg.base <- gg.base + ggplot2::geom_point(data = dataL[dataL$origin==FALSE], mapping = aes_string(x = "timeLeft", y = name.outcome, color = group.by, shape = "status"), size = 3)
        gg.base <- gg.base + ggplot2::scale_shape_manual(breaks = c(0,1), values=c(3,18), labels = c("censoring","event"))
    }else{
        gg.base <- gg.base + ggplot2::geom_point(data = dataL[dataL$origin==FALSE], mapping = aes_string(x = "timeLeft", y = name.outcome, color = group.by), size = 2)
    }
    
    if(group.by=="row"){
        gg.base <- gg.base + ggplot2::labs(color="observation") + theme(legend.key.height=unit(0.1,"npc"),
                                                                        legend.key.width=unit(0.08,"npc"))
        
        # display only integer values
        uniqueObs <- unique(dataL$row)

        if(length(uniqueObs)==1){
            gg.base <- gg.base + scale_color_continuous(guide=FALSE)
        }else{
            gg.base <- gg.base + scale_color_continuous(breaks = uniqueObs[seq(1,length(uniqueObs), length.out = min(10,length(uniqueObs)))],
                                                        limits = c(0.5, length(uniqueObs) + 0.5))
        }
    }

    if(is.na(alpha)[[1]] && (band[[1]] || ci[[1]])){
        indexTempo <- which(c(ci,band)==1)
        if(band && ci){
            value <- c(1,2)
        }else{
            value <- 1
        }
        gg.base <- gg.base + scale_linetype_manual("", breaks = c("ci","band")[indexTempo],
                                                   labels = c(labelCI,labelBand)[indexTempo],
                                                   values = value)
    }else if(ci[[1]] && band[[1]]){
        gg.base <- gg.base + ggplot2::guides(linetype = ggplot2::guide_legend(order = 1),
                                             fill = ggplot2::guide_legend(order = 2),
                                             group = ggplot2::guide_legend(order = 3)
                                             )
    }
    gg.base <- gg.base + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
    
    ## export
    ls.export <- list(plot = gg.base,
                      data = dataL)
    
    return(ls.export)
}
                                        # }}}

#----------------------------------------------------------------------
### autoplot.predictCox.R ends here
