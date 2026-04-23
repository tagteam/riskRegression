 ### autoplot.ate.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 28 2017 (14:19) 
## Version: 
## last-updated: Apr 23 2026 (11:56) 
##           By: Brice Ozenne
##     Update #: 264
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
#' \code{"ratioRisk"} to display the ratio of risks between treatment groups,
#' \code{"IPTW"} to display the treatment weights (if any) or \code{"probaT"} to display for each observation the probability of being assigned each treatment level,
#' and \code{"IPCW"} to display the censoring weights (if any).
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
#' library(ggplot2)
#' 
#' #### simulate data ####
#' n <- 1e2
#' set.seed(10)
#' dtS <- sampleData(n,outcome="survival")
#' seqTimes <- c(0,sort(dtS$time[dtS$event==1]),max(dtS$time))
#' 
#' #### Cox model ####
#' fit <- coxph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
#'
#' #### plot.type = 1: for few timepoints ####
#' ateFit <- ate(fit, data = dtS, treatment = "X1",
#'               times = c(1,2,5,10), se = TRUE, band = TRUE)
#' ggplot2::autoplot(ateFit)
#' \dontrun{
#' ggplot2::autoplot(ateFit, band = FALSE)
#' ggplot2::autoplot(ateFit, type = "diffRisk")
#' ggplot2::autoplot(ateFit, type = "ratioRisk")
#' }
#' 
#' #### plot.type = 2: when looking at all jump times ####
#' \dontrun{
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
#' outGGS <- ggplot2::autoplot(ateFit, plot.type = "2", alpha = NA, smooth = TRUE)
#' outGGS$plot + facet_wrap(~X1, labeller = label_both)
#' 
#' ## first derivative
#' ## (computation of the confidence intervals takes time)
#' ## (based on simulation - n.sim parameter)
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
    .data <- NULL ## [:: for CRAN CHECK::]

    ## ** initialize and check
    estimator <- match.arg(estimator, choices =  object$estimator, several.ok = FALSE)
    plot.type <- match.arg(as.character(plot.type), choices = c("1","2"), several.ok = FALSE)
    name.treatment <- object$variable["treatment"]

    type <- match.arg(type, c("meanRisk","diffRisk","ratioRisk","probaT","IPTW","IPCW"))
    if(type %in% c("IPTW", "probaT")){
        if(all(c("IPTW","AIPTW") %in% object$estimator == FALSE)){
            stop("No IPTW to display when the estimator is only G-formula \n",
                 "Consider providing a treatment model or setting argument \'estimator\' to \"IPTW\" or \"AIPTW\" when calling ate. \n")
        }
        if("probaT" %in% names(object$weights) == FALSE){
            stop("No IPTW values stored in the object  \n",
                 "Consider setting the argument \'store\' to c(weights = TRUE) to store the weights. \n")
        }
    }
    if(type == "IPCW"){
        if(all(c("IPTW","AIPTW") %in% object$estimator == FALSE)){
            stop("No IPCW to display when the estimator is only G-formula \n",
                 "Consider providing a treatment and a censoring model or setting argument \'estimator\' to \"IPTW\" or \"AIPTW\" when calling ate. \n")
        }
        if("probaT" %in% names(object$weights) == FALSE){
            stop("No IPCW values stored in the object  \n",
                 "Consider setting the argument \'store\' to c(weights = TRUE) to store the weights. \n")
        }
        if("probaC" %in% names(object$weights) == FALSE){
            stop("No IPCW values stored in the object \n",
                 "Likely because there was no censored observations. \n")
        }
    }
    

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

    ## ** create graph
    if(type == "IPTW"){

        IPTW.wide <- data.frame(id = 1:NROW(object$weights$probaT),
                                group = object$contrasts[max.col(object$weights$indicatorT)],
                                weight = rowSums(weights(object, type = "IPTW")))
        IPTW.wide$group <- as.factor(paste0(object$variable["treatment"],"=",IPTW.wide$group))

        gg <- ggplot2::ggplot(IPTW.wide, ggplot2::aes(x = .data$weight, fill = .data$group)) + ggplot2::geom_boxplot()
        gg <- gg + ggplot2::labs(fill = "Observations from", x = "Inverse Probability of Treatment Weight (IPTW)")

        gg.res <- list(data = IPTW.wide,
                       plot = gg)

    }else if(type == "IPCW"){

        IPCW.wide <- data.frame(id = 1:NROW(object$weights$probaC), group = object$contrasts[max.col(object$weights$indicatorT)], weights(object, type = "IPCW"))
        IPCW.long <- stats::reshape(IPCW.wide, direction = "long", idvar = c("id","group"),
                                    varying = 3:NCOL(IPCW.wide), v.names = "values", times = object$eval.times)
        IPCW.long$time <- as.factor(IPCW.long$time)
        IPCW.long$group <- as.factor(paste0(object$variable["treatment"],"=",IPCW.long$group))
        ## percentage of 0 at each timepoint
        IPCW.long0 <- as.data.frame(data.table::as.data.table(IPCW.long)[, list(mean=mean(.SD$values),
                                                                                label=paste0(round(100*mean(.SD$values>0),2),"%")),                           
                                                                         by = c("time","group")])
        IPCW.long0$y <- 0.9
        
        gg <- ggplot2::ggplot(mapping  = ggplot2::aes(.data$time, .data$values, fill = .data$group))
        gg <- gg + ggplot2::geom_boxplot(data = IPCW.long[IPCW.long$values>0,,drop=FALSE],
                                         position = ggplot2::position_dodge(.9))

        gg <- gg + ggplot2::geom_text(data = IPCW.long0, 
                             mapping = ggplot2::aes(x = .data$time, y = .data$y, color = .data$group, group = .data$group, label = .data$label),
                             position = ggplot2::position_dodge(.9))
        gg <- gg + ggplot2::labs(color = "Observations from", fill = "Observations from", y = "Inverse Probability of Censoring Weight (IPCW) \n (proportion of non-zero weights)")

        gg.res <- list(data = IPCW.long,
                       plot = gg)
        attr(gg.res$data,"zero") <- IPCW.long0

    }else if(type == "probaT"){

        IPTW.wide <- data.frame(1:NROW(object$weights$probaT), object$contrasts[max.col(object$weights$indicatorT)], weights(object, type = "probaT"))
        names(IPTW.wide) <- c("id","group",paste0(object$variable["treatment"],".",object$contrasts))
        IPTW.long <- stats::reshape(IPTW.wide, direction = "long",
                                    idvar = c("id","group"), varying = paste0(object$variable["treatment"],".",object$contrasts),
                                    v.names = "value", timevar = "target", times = object$contrasts)
        IPTW.long$target <- as.factor(IPTW.long$target)
        IPTW.long$group <- as.factor(paste0(object$variable["treatment"],"=",IPTW.long$group))

        gg <- ggplot2::ggplot(IPTW.long, ggplot2::aes(.data$value, .data$target, fill = .data$group)) + ggplot2::geom_boxplot()
        gg <- gg + ggplot2::labs(fill = "Observations from", y = paste0("Targeted ",object$variable["treatment"]), x = "Propensity score (1/IPTW)")

        gg.res <- list(data = IPTW.long,
                       plot = gg)

    }else{
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
                gg.res$plot <- ggplot2::ggplot(data = dataL, ggplot2::aes(x = .data$time, y = .data$estimate,
                                                                          ymin = .data$lowerBand, ymax = .data$upperBand, color = .data[[name.treatment]]))
                gg.res$plot <- gg.res$plot + ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.5))
            }else if(ci){
                if(test.ylab){ylab <- paste0(ylab, " (pointwise ci)")}
                gg.res$plot <- ggplot2::ggplot(data = dataL, ggplot2::aes(x = .data$time, y = .data$estimate,
                                                                          ymin = .data$lowerCI, ymax = .data$upperCI, color = .data[[name.treatment]]))
                gg.res$plot <- gg.res$plot + ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.5))
            }else{
                gg.res$plot <- ggplot(data = dataL, ggplot2::aes(x = .data$time, y = .data$estimate, color = .data[[name.treatment]]))
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
    }

    ## ** display
    if(plot){
        print(gg.res$plot)
    }


    ## ** export
    return(invisible(gg.res))
}



#----------------------------------------------------------------------
### autoplot.ate.R ends here
