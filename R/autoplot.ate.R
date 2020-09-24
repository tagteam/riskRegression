### autoplot.ate.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 28 2017 (14:19) 
## Version: 
## last-updated: sep 24 2020 (13:45) 
##           By: Brice Ozenne
##     Update #: 139
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
#' @param ... Additional parameters to cutomize the display.
#' 
#' @return Invisible. A list containing:
#' \itemize{
#' \item plot: the ggplot object.
#' \item data: the data used to create the plot.
#' }
#'
#' @seealso
#' \code{\link{ate}} to compute average risks.

## * autoplot.ate (examples)
#' @examples
#' library(survival)
#' library(rms)
#' library(ggplot2)
#' 
#' #### simulate data ####
#' n <- 1e2
#' set.seed(10)
#' dtS <- sampleData(n,outcome="survival")
#' seqTimes <- c(0,sort(dtS$time[dtS$event==1]),max(dtS$time))
#' 
#' #### Cox model ####
#' fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
#'
#' #### Average treatment effect ####
#' ateFit <- ate(fit, data = dtS, treatment = "X1",
#'               times = seqTimes, se = TRUE, band = TRUE)
#'
#' #### display #### 
#' ggplot2::autoplot(ateFit)
#' ggplot2::autoplot(ateFit, smooth = TRUE)
#' outGG <- ggplot2::autoplot(ateFit, alpha = 0.1)
#' 
#' dd <- as.data.frame(outGG$data[treatment == 0])
#' outGG$plot + facet_wrap(~treatment, labeller = label_both)

## * autoplot.ate (code)
#' @rdname autoplot.ate
#' @method autoplot ate
#' @export
autoplot.ate <- function(object,
                         estimator = object$estimator[1],
                         ci = object$inference$ci,
                         band = object$inference$band,
                         plot = TRUE,
                         smooth = FALSE,
                         digits = 2,
                         alpha = NA,
                         ...){

    ## initialize and check
    estimator <- match.arg(estimator, choices =  object$estimator, several.ok = FALSE)

    if(ci[[1]]==TRUE && object$inference$ci[[1]]==FALSE){
        stop("argument \'ci\' cannot be TRUE when no standard error have been computed \n",
             "set arguments \'se\' and \'confint\' to TRUE when calling the ate function \n")
    }
    if(band[[1]] && object$inference$band[[1]]==FALSE){
        stop("argument \'band\' cannot be TRUE when the quantiles for the confidence bands have not been computed \n",
             "set arguments \'band\' and \'confint\' to TRUE when calling the ate function \n")
    }
    if(any(rank(object$times) != 1:length(object$times))){
        stop("Invalid object. The prediction times must be strictly increasing \n")
    }
    
    ## dots <- list(...)
    ## if(length(dots)>0){
    ##     txt <- names(dots)
    ##     txt.s <- if(length(txt)>1){"s"}else{""}
    ##     stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    ## }

    ## display
    dataL <- object$meanRisk[estimator,.SD,on="estimator"]
    dataL[,row := as.numeric(as.factor(.SD$treatment))]
    if(ci){
        setnames(dataL, old = c("lower","upper"), new = c("lowerCI","upperCI"))
    }
    if(attr(object$estimator,"TD")){
        dataL$time <- dataL$time + dataL$landmark
    }
    gg.res <- predict2plot(dataL = dataL,
                           name.outcome = "estimate", # must not contain space to avoid error in ggplot2
                           ci = ci, band = band,
                           group.by = "treatment",
                           conf.level = object$conf.level,
                           smooth = smooth,
                           alpha = alpha,
                           xlab = "time",
                           ylab = "Average absolute risk",
                           ...)
  
    if(plot){
        print(gg.res$plot)
    }
  
    return(invisible(gg.res))
}



#----------------------------------------------------------------------
### autoplot.ate.R ends here
