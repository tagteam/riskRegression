### plot.predictCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 17 2017 (10:06) 
## Version: 
## last-updated: Mar  1 2017 (10:01) 
##           By: Thomas Alexander Gerds
##     Update #: 112
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Plot predictions from a Cox model
#' @description Plot predictions from a Cox model
#' 
#' @param x object obtained with the function \code{predictCox}.
#' @inheritParams predictCox
#' @param ci Logical. If \code{TRUE} display the confidence intervals for the predictions.
#' @param groupBy The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}.
#' @param reduce.data Logical. If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param plot Logical. Should the graphic be plotted.
#' @param conf.level confidence level of the interval.
#' @param digit integer indicating the number of decimal places
#' @param ... not used. Only for compatibility with the plot method.
#' 
#' @examples
#' library(survival)
#' library(ggplot2)
#' 
#' d <- sampleData(1e2, outcome = "survival")
#' m.cox <- coxph(Surv(time,event)~ X1 + X2 + X3,
#'                 data = d, x = TRUE, y = TRUE)
#' dt.basehaz <- predictCox(m.cox)
#' ggplot(dt.basehaz, aes(x = time, y = survival)) + geom_point() + geom_line()
#'
#' pred.cox <- predictCox(m.cox, newdata = d[1:4,],
#'   times = 1:5, type = "survival", se = TRUE, keep.newdata = TRUE)
#' plot(pred.cox,ci=TRUE)
#' plot(pred.cox, groupBy = "covariates")
#' plot(pred.cox, groupBy = "covariates", reduce.data = TRUE)
#' 
#' 
#' m.cox.strata <- coxph(Surv(time,event)~ strata(X1) + strata(X2) + X3 + X6,
#' data = d, x = TRUE, y = TRUE)
#' pred.cox.strata <- predictCox(m.cox.strata, newdata = d[c(1:5,10,50),],
#' time = 1:5, keep.newdata = TRUE)
#' plot(pred.cox.strata, type = "survival")
#' plot(pred.cox.strata, type = "survival", groupBy = "strata")
#' res <- plot(pred.cox.strata, type = "survival",
#'             groupBy = "covariates")
#'
#' # customize display
#' res$plot + geom_point(size = 3)
#'
#' @method plot predictCox
#' @export
plot.predictCox <- function(x,
                            type = NULL,
                            ci = FALSE,
                            groupBy = "row",
                            reduce.data = FALSE,
                            plot = TRUE,
                            conf.level = 0.95,
                            digit = 2, ...){

    ## initialize and check    
    possibleType <- c("hazard","cumhazard","survival")
    possibleType <- possibleType[possibleType %in% names(x)]

    if(is.null(type)){
        if(length(possibleType) == 1){
            type <- possibleType
        }else{
            stop("argument \'type\' must be specified to choose between ",paste(possibleType, collapse = " "),"\n")
        }
    }else if(length(type)>1){
        stop("argument \'type\' must have length 1 \n")        
    }else if(type %in% possibleType == FALSE){
        stop("argument \'type\' can only be ",paste(possibleType, collapse = " or ")," \n")        
    }
    typename <- switch(type,
                       hazard = "hazard",
                       cumhazard = "cumulative hazard",
                       survival = "survival")
    
    possibleGroupBy <- c("row","covariates","strata")
    if(groupBy %in% possibleGroupBy == FALSE){
        stop("argument \"groupBy\" must be in \"",paste(possibleGroupBy, collapse = "\" \""),"\"\n")
    }
    
    if(groupBy == "covariates" && ("newdata" %in% names(x) == FALSE)){
        stop("argument \'groupBy\' cannot be \"covariates\" when newdata is missing in the object \n",
             "set argment \'keep.newdata\' to TRUE when calling predictCox \n")
    }
    if(groupBy == "strata" && ("strata" %in% names(x) == FALSE)){
        stop("argument \'groupBy\' cannot be \"strata\" when strata is missing in the object \n",
             "set argment \'keep.strata\' to TRUE when calling predictCox \n")
    }

    if(ci && (paste0(type,".se") %in% names(x) == FALSE)){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set argment \'se\' to TRUE when calling predictCox \n")
    }

    ## display
    newdata <- copy(x$newdata)
    if(!is.null(newdata) && reduce.data){
        test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
        if(any(test)){
            newdata[, (names(test)[test]):=NULL]
        }        
    }
    
    gg.res <- predict2plot(outcome = x[[type]],
                           outcome.se = if(ci){x[[paste0(type,".se")]]}else{NULL},
                           newdata = newdata,
                           strata = x$strata,
                           times = x$times,
                           digit = digit,
                           name.outcome = typename,
                           conf.level = conf.level,
                           groupBy = groupBy)

    if(plot){
        print(gg.res)
    }
    
    return(invisible(list(plot = gg.res,
                          data = newdata)))
}

predict2plot <- function(outcome, outcome.se, newdata, strata, times,
                         digit, name.outcome, conf.level, groupBy){

    n.obs <- NROW(outcome)
    n.time <- NCOL(outcome)
    if(!is.null(time)){
        time.names <- times
    }else{
        time.names <- 1:n.time
    }    
    colnames(outcome) <- paste0(name.outcome,"_",time.names)
    

    if(!is.null(outcome.se)){
        pattern <- c(paste0(name.outcome,"_"),"lower_","upper_")

        lower <- upper <- matrix(NA, nrow = n.obs, ncol = n.time)
        lower[] <- pmax(0,outcome + qnorm((1-conf.level)/2) * outcome.se)
        upper[] <- pmin(1,outcome + qnorm(1-(1-conf.level)/2) * outcome.se)
        colnames(lower) <- paste0("lower_",time.names)
        colnames(upper) <- paste0("upper_",time.names)
        outcome <- data.table::as.data.table(cbind(outcome,lower,upper))
    }else{
        outcome <- data.table::as.data.table(outcome)
        pattern <- paste0(name.outcome,"_")
    }

    outcome[, row := 1:.N]
    if(groupBy == "covariates"){
        cov.names <- names(newdata)
        newdata <- newdata[, (cov.names) := lapply(cov.names, function(col){paste0(col,"=",round(.SD[[col]],digit))})]
        outcome[, ("covariates") := interaction(newdata,sep = " ")]
    }else if(groupBy == "strata"){
        outcome[, strata := strata]
    }
    gg.dtL <- melt(outcome, id.vars = union("row",groupBy),
                   variable.name = "time", value.name = gsub("_","",pattern))
    gg.dtL[, time := gsub(pattern[1],"",time)]
    if(!is.null(times)){
        gg.dtL[, time := as.numeric(time)]
    }
    gg.base <- ggplot(data = gg.dtL, aes_string(x = "time",y = name.outcome, group = "row", color = groupBy))

    gg.base <- ggplot(data = gg.dtL, aes_string(x = "time", y = name.outcome, group = "row", color = groupBy))
    gg.base <- gg.base + geom_point() + geom_line()
    if(!is.null(outcome.se)){
        gg.base <- gg.base + geom_errorbar(aes(ymin = lower, ymax = upper))
    }
    
    return(gg.base)
}


#----------------------------------------------------------------------
### plot.predictCox.R ends here
