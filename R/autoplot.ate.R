### autoplot.ate.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 28 2017 (14:19) 
## Version: 
## last-updated: okt 26 2020 (09:57) 
##           By: Brice Ozenne
##     Update #: 188
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
#' @param type [character vector] what to displayed.
#' Can be \code{"meanRisk"} to display the risks specific to each treatment group,
#' \code{"diffRisk"} to display the difference in risks between treatment groups,
#' or \code{"ratioRisk"} to display the ratio of risks between treatment groups,.
#' @param ci [logical] If \code{TRUE} display the confidence intervals for the average risks.
#' @param band [logical] If \code{TRUE} display the confidence bands for the average risks.
#' @param plot [logical] Should the graphic be plotted.
#' @param plot.type [character] Type of plot to be used.
#' \code{plot.type="2"} is useful when looking simulateneous at all eventtimes.
#' Otherwise use \code{plot.type="1"}.
#' @param digits [integer, >0] Number of decimal places.
#' @param alpha [numeric, 0-1] Transparency of the confidence bands. Argument passed to \code{ggplot2::geom_ribbon}.
#' @param ylab [character] Label for the y axis.
#' @param first.derivative [logical] If \code{TRUE}, display the first derivative over time of the risks/risk differences/risk ratios.
#' (confidence intervals are obtained via simulation).
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
#' #### plot.type = 1: for few timepoints ####
#' ateFit <- ate(fit, data = dtS, treatment = "X1",
#'               times = c(1,2,5,10), se = TRUE, band = TRUE)
#' ggplot2::autoplot(ateFit)
#' ggplot2::autoplot(ateFit, band = FALSE)
#' ggplot2::autoplot(ateFit, type = "diffRisk")
#' ggplot2::autoplot(ateFit, type = "ratioRisk")
#'
#' #### plot.type = 2: when looking at all jump times ####
#' ateFit <- ate(fit, data = dtS, treatment = "X1",
#'               times = seqTimes, se = TRUE, band = TRUE)
#'
#' ggplot2::autoplot(ateFit, plot.type = "2")
#' 
#' ## customize plot
#' outGG <- ggplot2::autoplot(ateFit, plot.type = "2", alpha = 0.25)
#' outGG$plot + facet_wrap(~X1, labeller = label_both)
#'
#' 
#' ## Looking at the difference after smoothing
#' \dontrun{
#' outGGS <- ggplot2::autoplot(ateFit, plot.type = "2", alpha = NA, smooth = TRUE)
#' outGGS$plot + facet_wrap(~X1, labeller = label_both)
#' }
#' 
#' ## first derivative
#' ## (computation of the confidence intervals takes time)
#' ## (based on simulation - n.sim parameter)
#' \dontrun{ 
#' ggplot2::autoplot(ateFit, plot.type = "2", smooth = TRUE,
#'                   band = FALSE, type = "diffRisk")
#' ggplot2::autoplot(ateFit, plot.type = "2", smooth = TRUE, first.derivative = TRUE,
#'                   band = FALSE, type = "diffRisk")
#' }

## * autoplot.ate (code)
#' @rdname autoplot.ate
#' @method autoplot ate
#' @export
autoplot.ate <- function(object,
                         type = "meanRisk",
                         first.derivative = FALSE,
                         estimator = object$estimator[1],
                         ci = object$inference$ci,
                         band = object$inference$band,
                         plot.type = "1",
                         plot = TRUE,
                         smooth = FALSE,
                         digits = 2,
                         alpha = NA,
                         ylab = NULL,
                         ...){

    ## initialize and check
    estimator <- match.arg(estimator, choices =  object$estimator, several.ok = FALSE)
    plot.type <- match.arg(as.character(plot.type), choices = c("1","2"), several.ok = FALSE)
    name.treatment <- object$variable["treatment"]
    
    type <- match.arg(type, c("meanRisk","diffRisk","ratioRisk"))
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
    if(first.derivative && plot.type == "1"){
        stop("Argument \'first.derivative\' currently only implemented for plot.type=\"1\" \n")
    }
    ## dots <- list(...)
    ## if(length(dots)>0){
    ##     txt <- names(dots)
    ##     txt.s <- if(length(txt)>1){"s"}else{""}
    ##     stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    ## }

    ## display
    test.ylab <- is.null(ylab)
    dataL <- object[[type]][estimator,.SD,on="estimator"]
    if(type %in% c("meanRisk")){
        data.table::setnames(dataL, old = "treatment", new = name.treatment)
        if(first.derivative && test.ylab){
            ylab <- "Derivative of the average absolute risk"
        }else if(is.null(ylab)){
            ylab <- "Average absolute risk"
        }
    }else if(type %in% c("diffRisk")){
        dataL[, c(name.treatment) :=  paste0(.SD$B,"-",.SD$A)]
        if(first.derivative && test.ylab){
            ylab <- "Derivative of the difference in average absolute risks"
        }else if(is.null(ylab)){
            ylab <- "Difference in average absolute risks"
        }
    }else if(type %in% c("ratioRisk")){
        dataL[, c(name.treatment) :=  paste0(.SD$B,"/",.SD$A)]
        if(first.derivative && test.ylab){
            ylab <- "Derivative of the ratio between average absolute risks"
        }else if(is.null(ylab)){
            ylab <- "Ratio between average absolute risks"
        }
    }
    
    dataL[,row := as.numeric(as.factor(.SD[[name.treatment]]))]
    if(ci){
        setnames(dataL, old = c("lower","upper"), new = c("lowerCI","upperCI"))
    }
    if(attr(object$estimator,"TD")){
        dataL$time <- dataL$time + dataL$landmark
    }
    if(first.derivative && ci){
        attr(first.derivative,"vcov") <- attr(object[[type]],"vcov")[[estimator]]
    }

    if(plot.type=="1"){
        dataL$time <- as.factor(dataL$time)
        gg.res <- list(data = dataL,
                       plot = NULL)
        if(band>0){
            if(test.ylab){
                if(band==1){ylab <- paste0(ylab, " (simultaneous ci over time)")}
                if(band==2){ylab <- paste0(ylab, " (simultaneous ci over time and treatment)")}
            }
            gg.res$plot <- ggplot2::ggplot(data = dataL, aes_string(x = "time", y = "estimate", ymin = "lowerBand", ymax = "upperBand", color = name.treatment))
            gg.res$plot <- gg.res$plot + ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.5))
        }else if(ci){
            if(test.ylab){ylab <- paste0(ylab, " (pointwise ci)")}
            gg.res$plot <- ggplot2::ggplot(data = dataL, ggplot2::aes_string(x = "time", y = "estimate", ymin = "lowerCI", ymax = "upperCI", color = name.treatment))
            gg.res$plot <- gg.res$plot + ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.5))
        }else{
            gg.res$plot <- ggplot(data = dataL, ggplot2::aes_string(x = "time", y = "estimate", color = name.treatment))
            gg.res$plot <- gg.res$plot + ggplot2::geom_point()
        }
        gg.res$plot <- gg.res$plot + ggplot2::labs(y = ylab)
    }else if(plot.type=="2"){
        gg.res <- predict2plot(dataL = dataL,
                               name.outcome = "estimate", # must not contain space to avoid error in ggplot2
                               ci = ci, band = band,
                               group.by = name.treatment,
                               conf.level = object$inference$conf.level,
                               smooth = smooth,
                               alpha = alpha,
                               xlab = "time",
                               ylab = ylab,
                               first.derivative = first.derivative,
                               ...)
    }
    
    if(plot){
        print(gg.res$plot)
    }
  
    return(invisible(gg.res))
}



#----------------------------------------------------------------------
### autoplot.ate.R ends here
