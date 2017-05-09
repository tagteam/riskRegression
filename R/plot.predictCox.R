### plot.predictCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: feb 17 2017 (10:06) 
## Version: 
## last-updated: apr 28 2017 (15:36) 
##           By: Brice Ozenne
##     Update #: 258
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ plot.predictCox
#' @title Plot predictions from a Cox model
#' @description Plot predictions from a Cox model
#' 
#' @param x object obtained with the function \code{predictCox}.
#' @inheritParams predictCox
#' @param ci Logical. If \code{TRUE} display the confidence intervals for the predictions.
#' @param band Logical. If \code{TRUE} display the confidence bands for the predictions.
#' @param groupBy The grouping factor used to color the prediction curves. Can be \code{"row"}, \code{"strata"}, or \code{"covariates"}.
#' @param reduce.data Logical. If \code{TRUE} only the covariates that does take indentical values for all observations are displayed.
#' @param plot Logical. Should the graphic be plotted.
#' @param digits integer indicating the number of decimal places
#' @param alpha transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param ... not used. Only for compatibility with the plot method.
#' 
#' @examples
#' library(survival)
#' library(ggplot2)
#'
#' ## predictions ##
#' d <- sampleData(1e2, outcome = "survival")
#' m.cox <- coxph(Surv(time,event)~ X1 + X2 + X3,
#'                 data = d, x = TRUE, y = TRUE)
#' dt.basehaz <- predictCox(m.cox)
#' ggplot(as.data.table(dt.basehaz), aes(x = time, y = survival)) + geom_point() + geom_line()
#'
#' pred.cox <- predictCox(m.cox, newdata = d[1:4,],
#'   times = 1:5, type = "survival", keep.newdata = TRUE)
#' plot(pred.cox)
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
#' res$plot + geom_point(data = res$dataDuplicated, size = 3)
#'
#' ## predictions with confidence interval
#' pred.cox <- predictCox(m.cox, newdata = d[1:4,],
#'   times = 1:5, type = "survival", se = TRUE, keep.newdata = TRUE)
#' plot(pred.cox, ci = TRUE)
#'
#' ## predictions with confidence bands
#' pred.cox <- predictCox(m.cox, newdata = d[1,,drop=FALSE],
#'   times = 1:5, type = "survival", nSim.band = 500, keep.newdata = TRUE)
#' plot(pred.cox, band = TRUE)
#'
#' 
#' @method plot predictCox
#' @export
plot.predictCox <- function(x,
                            type = NULL,
                            ci = FALSE,
                            band = FALSE,
                            groupBy = "row",
                            reduce.data = FALSE,
                            plot = TRUE,
                            digits = 2, alpha = 0.1, ...){
  
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

    if(ci && (paste0(type,".se") %in% names(x) == FALSE)){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set argment \'se\' to TRUE when calling predictCox \n")
    }

    if(band && ("quantile.band" %in% names(x) == FALSE)){
        stop("argument \'band\' cannot be TRUE when no quantiles for the confidence bands have not been computed \n",
             "set argment \'nSim.band\' to a positive integer when calling predictCox \n")
    }
    
    ## display
    newdata <- copy(x$newdata)
    if(!is.null(newdata) && reduce.data){
        test <- unlist(newdata[,lapply(.SD, function(col){length(unique(col))==1})])
        if(any(test)){
            newdata[, (names(test)[test]):=NULL]
        }        
    }

    dataL <- predict2melt(outcome = x[[type]], ci = ci, band = band,
                          outcome.lower = if(ci){x[[paste0(type,".lower")]]}else{NULL},
                          outcome.upper = if(ci){x[[paste0(type,".upper")]]}else{NULL},
                          outcome.lowerBand = if(band){x[[paste0(type,".lowerBand")]]}else{NULL},
                          outcome.upperBand = if(band){x[[paste0(type,".upperBand")]]}else{NULL},
                          newdata = newdata,
                          strata = x$strata,
                          times = x$times,
                          name.outcome = typename,
                          groupBy = groupBy,
                          digits = digits
                          )

    gg.res <- predict2plot(dataL = dataL,
                           name.outcome = typename,
                           ci = ci,
                           band = band,
                           groupBy = groupBy,
                           conf.level = x$conf.level,
                           alpha = alpha
                           )
  
  if(plot){
    print(gg.res$plot)
  }
  
  return(invisible(gg.res))
}
# }}}

# {{{ predict2melt
predict2melt <- function(outcome, name.outcome,
                         ci, outcome.lower, outcome.upper,
                         band, outcome.lowerBand, outcome.upperBand,
                         newdata, strata, times, groupBy, digits){

    ## for CRAN tests
    patterns <- function(){}
    
    n.time <- NCOL(outcome)
    if(!is.null(time)){
        time.names <- times 
    }else{
        time.names <- 1:n.time
    }    
    colnames(outcome) <- paste0(name.outcome,"_",time.names)
    keep.cols <- unique(c("time",name.outcome,"row",groupBy))
    
    #### merge outcome with CI and band ####
    pattern <- paste0(name.outcome,"_")
    if(ci){
        pattern <- c(pattern,"lowerCI_","upperCI_")
    
        colnames(outcome.lower) <- paste0("lowerCI_",time.names)
        colnames(outcome.upper) <- paste0("upperCI_",time.names)
    }
    if(band){
        pattern <- c(pattern,"lowerBand_","upperBand_")
        keep.cols <- c(keep.cols,"lowerBand","upperBand")
        
        colnames(outcome.lowerBand) <- paste0("lowerBand_",time.names)
        colnames(outcome.upperBand) <- paste0("upperBand_",time.names)
    }

    outcome <- data.table::as.data.table(
                               cbind(outcome,
                                     outcome.lower, outcome.upper,
                                     outcome.lowerBand,outcome.upperBand)
                           )

    #### merge with convariates ####
    outcome[, row := 1:.N]
    if(groupBy == "covariates"){
        cov.names <- names(newdata)
        newdata <- newdata[, (cov.names) := lapply(cov.names, function(col){paste0(col,"=",round(.SD[[col]],digits))})]
        outcome[, ("covariates") := interaction(newdata,sep = " ")]
    }else if(groupBy == "strata"){
        outcome[, strata := strata]
    }
    
    #### reshape to long format ####
    dataL <- melt(outcome, id.vars = union("row",groupBy),
                   measure= patterns(pattern),
                   variable.name = "time", value.name = gsub("_","",pattern))
    dataL[, time := as.numeric(as.character(factor(time, labels = time.names)))]

    return(dataL)    
}

# }}}
# {{{ predict2plot
predict2plot <- function(dataL, name.outcome,
                         ci, band, groupBy,                         
                         conf.level, alpha){

    # for CRAN tests
    original <- lowerCI <- upperCI <- lowerBand <- upperBand <- NULL
    #### duplicate observations to obtain step curves ####
    keep.cols <- unique(c("time",name.outcome,"row",groupBy))
    if(band){
        keep.cols <- c(keep.cols,"lowerBand","upperBand")
    }
    dataL[, original := TRUE]

    dtTempo <- copy(dataL)
    dtTempo[, (c("time","original")) := list(time = c(0,.SD$time[-.N] + .Machine$double.eps*100),
                                             original = FALSE),
            by = row]

    dataL_duplicated <- rbind(dataL[,unique(keep.cols), with = FALSE],
                              dtTempo[,unique(keep.cols), with = FALSE])

    #### display ####

    gg.base <- ggplot(mapping = aes_string(x = "time", y = name.outcome, group = "row", color = groupBy))
    gg.base <- gg.base + geom_point(data = dataL_duplicated) + geom_line(data = dataL_duplicated)
    if(groupBy=="row"){
        gg.base <- gg.base + ggplot2::labs(color="observation") + theme(legend.key.height=unit(0.1,"npc"),
                                                                        legend.key.width=unit(0.08,"npc"))
        
        # display only integer values
        uniqueObs <- unique(dataL$row)
        gg.base <- gg.base + scale_color_continuous(breaks = uniqueObs[seq(1,length(uniqueObs), length.out = min(10,length(uniqueObs)))],
                                                    limits = c(0.5, length(uniqueObs) + 0.5))
    }
    if(ci){
        labelCI <- paste0(conf.level*100,"% confidence interval")
        gg.base <- gg.base + geom_errorbar(data = dataL, aes(ymin = lowerCI, ymax = upperCI, linetype = labelCI))
        gg.base <- gg.base + scale_linetype_manual("",values=setNames(1,labelCI))
    }
    if(band){
        labelBand <- paste0(conf.level*100,"% confidence band")
        gg.base <- gg.base + geom_ribbon(data = dataL_duplicated, aes(ymin = lowerBand, ymax = upperBand, fill = labelBand), alpha = alpha)
        gg.base <- gg.base + scale_fill_manual("", values="grey12")        
    }
    
    if(ci && band){
        gg.base <- gg.base + ggplot2::guides(group = ggplot2::guide_legend(order = 1),
                                             se = ggplot2::guide_legend(order = 2),
                                             fill = ggplot2::guide_legend(order = 3))
    }
    
    ## export
    ls.export <- list(plot = gg.base,
                      data = dataL,
                      dataDuplicated = dataL_duplicated)
    
    return(ls.export)
}
# }}}

#----------------------------------------------------------------------
### plot.predictCox.R ends here
