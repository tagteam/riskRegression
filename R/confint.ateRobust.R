### confint.influenceTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun  1 2018 (12:15) 
## Version: 
## Last-Updated: sep 19 2018 (08:08) 
##           By: Brice Ozenne
##     Update #: 57
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * confint.ateRobust (documentation)
##' @title Confidence Intervals for the Average Treatment Effect
##' @description Confidence intervals for the average treatment effect
##' @name confint.ateRobust
##' 
##' @param object A \code{ateRobust} object, i.e. output of the \code{ateRobust} function.
##' @param parm not used. For compatibility with the generic method.
##' @param level [numeric, 0-1] Level of confidence.
##' @param meanRisk.transform [character] the transformation used
##' to improve coverage of the confidence intervals for the mean risk in small samples.
##' Can be \code{"none"}, \code{"log"}, \code{"loglog"}, \code{"cloglog"}.
##' @param diffRisk.transform [character] the transformation used to improve coverage
##' of the confidence intervals for the risk difference in small samples.
##' Can be \code{"none"}, \code{"atanh"}.
##' @param ... not used.
##'
##' @author Brice Ozenne
##'

## * confint.ateRobust (code)
##' @rdname confint.ateRobust
##' @method confint ateRobust
##' @export
confint.ateRobust <- function(object,
                              parm = NULL,
                              level = 0.95,
                              meanRisk.transform = "none",
                              diffRisk.transform = "none",
                              ...){


    ## ** check arguments
    dots <- list(...)
    if(length(dots)>0){
        txt <- names(dots)
        txt.s <- if(length(txt)>1){"s"}else{""}
        stop("unknown argument",txt.s,": \"",paste0(txt,collapse="\" \""),"\" \n")
    }

    object$meanRisk.transform <- match.arg(meanRisk.transform, c("none","log","loglog","cloglog"))
    object$diffRisk.transform <- match.arg(diffRisk.transform, c("none","atanh"))

    ## ** compute se, CI/CB
    ate.se <- object$ate.se
    if(object$se==FALSE){
         ate.se[,grep("Gformula|^IPTW",colnames(ate.se))] <- NA
    }
   
    outCIBP.risk <- transformCIBP(estimate = object$ate.value[c("risk.0","risk.1"),],
                                  se = ate.se[c("risk.0","risk.1"),],
                                  null = 0,
                                  conf.level = level,
                                  type = object$meanRisk.transform,
                                  min.value = switch(object$meanRisk.transform,
                                                     "none" = 0,
                                                     "log" = NULL,
                                                     "loglog" = NULL,
                                                     "cloglog" = NULL),
                                  max.value = switch(object$meanRisk.transform,
                                                     "none" = 1,
                                                     "log" = 1,
                                                     "loglog" = NULL,
                                                     "cloglog" = NULL),
                                  ci = TRUE,
                                  band = FALSE,
                                  p.value = FALSE)

    outCIBP.ate <- transformCIBP(estimate = object$ate.value["ate.diff",,drop=FALSE],
                                 se = ate.se["ate.diff",,drop=FALSE],
                                 null = 0,
                                 conf.level = level,
                                 type = diffRisk.transform,
                                 min.value = switch(diffRisk.transform,
                                                    "none" = 0,
                                                    "atanh" = NULL),
                                 max.value = switch(diffRisk.transform,
                                                    "none" = 1,
                                                    "atanh" = NULL),
                                 ci = TRUE,
                                 band = FALSE,
                                 p.value = TRUE)

    ## export
    object$ate.lower <- rbind(outCIBP.risk$lower,outCIBP.ate$lower)
    object$ate.upper <- rbind(outCIBP.risk$upper,outCIBP.ate$upper)
    object$ate.p.value <- rbind(risk.0 = NA, risk.1 = NA,outCIBP.ate$p.value)
    object$conf.level <- level
    return(object)
}



######################################################################
### confint.influenceTest.R ends here
