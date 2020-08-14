### autoplot.ate.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 28 2017 (14:19) 
## Version: 
## last-updated: aug 14 2020 (11:47) 
##           By: Brice Ozenne
##     Update #: 118
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * autoplot.ate (documentation)
#' @title Plot Average Risks
#' @description Plot average risks.
#' @name autoplot.ate
#' 
#' @param object Object obtained with the function \code{ate}.
#' @param ci [logical] If \code{TRUE} display the confidence intervals for the average risks.
#' @param band [logical] If \code{TRUE} display the confidence bands for the average risks.
#' @param plot [logical] Should the graphic be plotted.
#' @param digits [integer, >0] Number of decimal places.
#' @param alpha [numeric, 0-1] Transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param smooth [logical] Should a smooth version of the risk function be plotted instead of a simple function?
#' @param estimator [character] The type of estimator relative to which the risks should be displayed. 
#' @param ... not used. Only for compatibility with the plot method.
#' 
#' @seealso
#' \code{\link{ate}} to compute average risks.

## * autoplot.ate (examples)
#' @examples
#' \dontrun{
#' library(survival)
#' library(rms)
#' library(ggplot2)
#' #### simulate data ####
#' n <- 1e2
#' set.seed(10)
#' dtS <- sampleData(n,outcome="survival")
#' 
#' 
#' #### Cox model ####
#' fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
#'
#' #### Average treatment effect ####
#' seqTimes <- sort(unique(fit$y[,1]))
#' seqTimes5 <- seqTimes[seqTimes>5 & seqTimes<10]
#' ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
#'               times = seqTimes, B = 0, band = TRUE, nsim.band = 500, y = TRUE,
#'               mc.cores=1)
#'
#' #### display #### 
#' ggplot2::autoplot(ateFit)
#' 
#' outGG <- autoplot(ateFit, band = TRUE, ci = TRUE, alpha = 0.1)
#' dd <- as.data.frame(outGG$data[treatment == 0])
#' outGG$plot + facet_wrap(~treatment, labeller = label_both)
#' }

## * autoplot.ate (code)
#' @rdname autoplot.ate
#' @method autoplot ate
#' @export
autoplot.ate <- function(object,
                         estimator = object$estimator[1],
                         ci = FALSE,
                         band = FALSE,
                         plot = TRUE,
                         smooth = FALSE,
                         digits = 2,
                         alpha = NA,
                         ...){

    ## initialize and check
    estimator <- match.arg(estimator, choices =  object$estimator, several.ok = FALSE)

    if(ci[[1]]==TRUE && (object$se[[1]]==FALSE || is.null(object$conf.level))){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set arguments \'se\' and \'confint\' to TRUE when calling the ate function \n")
    }
    if(band[[1]] && (object$band[[1]]==FALSE  || is.null(object$conf.level))){
        stop("argument \'band\' cannot be TRUE when the quantiles for the confidence bands have not been computed \n",
             "set arguments \'band\' and \'confint\' to TRUE when calling the ate function \n")
    }
    
    dots <- list(...)
    if(length(dots)>0){
        txt <- names(dots)
        txt.s <- if(length(txt)>1){"s"}else{""}
        stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    }

    ## display
    dataL <- copy(object$meanRisk)
    dataL[,row := as.numeric(as.factor(.SD$treatment))]
    dataL[,c("origin") := FALSE]
    data.table::setnames(dataL, old = paste0("meanRisk.",estimator), new = "meanRisk")
    if(min(dataL$time)>1e-12){
        first.dt <- data.table::data.table(time = c(0,min(dataL$time)-1e-12), meanRisk = 0, origin = TRUE)
    }else{
        first.dt <- NULL
    }
    
    if(ci){
        data.table::setnames(dataL,
                             old = c(paste0("meanRisk.",estimator,".lower"),paste0("meanRisk.",estimator,".upper")),
                             new = c("lowerCI","upperCI"))
        if(!is.null(first.dt)){
            first.dt[, c("lowerCI","upperCI") := 0]
        }
    }
    if(band){
        data.table::setnames(dataL,
                             old = c(paste0("meanRisk.",estimator,".lowerBand"),paste0("meanRisk.",estimator,".upperBand")),
                             new = c("lowerBand","upperBand"))
        if(!is.null(first.dt)){
            first.dt[, c("lowerBand","upperBand") := 0]
        }
    }
    if(!is.null(first.dt)){
        dataL <- dataL[, .SD, .SDcol = c(names(first.dt),c("row","treatment"))]
        dataL <- dataL[, rbind(first.dt,.SD), by = c("row","treatment")]
    }else{
        dataL <- dataL[, .SD, .SDcol = c("time","meanRisk","row","treatment", if(ci){c("lowerCI","upperCI")}, if(band){c("lowerBand","upperBand")})]
        dataL[,c("origin") := TRUE]
    }
    gg.res <- predict2plot(dataL = dataL,
                           name.outcome = "meanRisk", # must not contain space to avoid error in ggplot2
                           ci = ci, band = band,
                           group.by = "treatment",
                           conf.level = object$conf.level,
                           smooth = smooth,
                           alpha = alpha,
                           xlab = "time",
                           ylab = "Average absolute risk")
  
    if(plot){
        print(gg.res$plot)
    }
  
    return(invisible(gg.res))
}



#----------------------------------------------------------------------
### autoplot.ate.R ends here
