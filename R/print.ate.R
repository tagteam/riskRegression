### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: sep 23 2020 (17:22) 
##           By: Brice Ozenne
##     Update #: 498
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * print.ate (documentation)
#' @title Print Average Treatment Effects
#' @description Print average treatment effects.
#' @name print.ate
#' 
#' @param x object obtained with function \code{ate}
#' @param ... for internal use
#'
#' @seealso
#' \code{\link{summary.ate}} to obtained a more detailed output
#' \code{\link{confint.ate}} to compute confidence intervals/bands.
#' \code{\link{ate}} to compute the average treatment effects.

## * print.ate (code)
#' @rdname print.ate
#' @method print ate
#' @export
print.ate <- function(x, estimator = x$estimator, ...){

    dots <- list(...)

    estimator.print <- sapply(estimator,
                              switch,
                              "GFORMULA" = "G-formula",
                              "GFORMULATD" = "G-formula with time-dependent covariates)",
                              "IPTW,IPCW" = "Inverse probability of treatment weighting",
                              "AIPTW,AIPCW" = "Augmented estimator")


    ## ** prepare
    post.calculation <- !is.null(x$meanRisk)
    if(!is.null(dots$display.results)){
        display.results <- FALSE
    }else{
        display.results <- post.calculation
    }
    ## First line
    cause <- attr(x$theCause,"cause")
    n.cause <- length(cause)
    if(n.cause>1){
        txt.cause <- paste0(" for cause ",x$theCause)
    }else{
        txt.cause <- ""
    }
    ## event variable
    txt.eventvar <- paste0("cause: ",cause)
    if(length(setdiff(cause,x$theCause))>0){
        txt.eventvar <- paste0(txt.eventvar,", competing risk(s): ",paste(setdiff(cause,x$theCause),collapse = " "))
    }
    if(any(attr(x$eval.times,"n.censored")>0)){
        txt.eventvar <- paste0(txt.eventvar, ", censoring: ",attr(x$theCause,"level.censoring"))
    }
    ## number at risk at the evaluation time
    M.at.risk <- cbind("- Eval. time       "="number at risk",attr(x$eval.times,"n.at.risk"))
    colnames(M.at.risk)[2] <- " :"
    M.at.risk[,2] <- paste(as.character(M.at.risk[[2]]),"  ",sep="")

    ## statistical inference
    if(!is.na(x$inference$B) && !identical(x$inference$B,0)){
        bootci.method <- switch(x$inference$bootci.method,
                                "norm" = "Normal",
                                "basic" = "Basic",
                                "stud" = "Studentized",
                                "perc" = "Percentile",
                                "bca" = "BCa",
                                "wald" = "Wald",
                                "quantile" = "Percentile")
        txt.inference <- paste(bootci.method," bootstrap based on ",x$inference$B," bootstrap samples\n",
                               "                          that were drawn with replacement from the original data.\n",sep="")
        test.inference <- TRUE
    }else if(any(x$inference[,c("se","band","ci","p.value")]>0)){
        txt.inference <- "decomposition in iid terms \n"
        test.inference <- TRUE
    }else{
        test.inference <- FALSE
    }
 

        
    ## ** input variables
    if(!is.na(x$variable["strata"])){
        if(post.calculation){
            cat("     Average risk",txt.cause," by strata \n\n", sep = "")
        }else{
            cat(" Input variables \n")
        }
        
        cat(" - Strata               : ",x$variable["strata"]," (",length(x$contrasts)," levels: \"",paste0(x$contrasts,collapse ="\" \""),"\")\n",sep="")
    }else if(!is.na(x$variable["treatment"])){
        if(post.calculation){
            cat("     Average treatment effect",txt.cause," \n\n", sep = "")
        }else{
            cat(" Input variables \n")
        }
        cat(" - Treatment            : ",x$var["treatment"]," (",length(x$contrasts)," levels: \"",paste0(x$contrasts,collapse ="\" \""),"\")\n",sep="")
    }

    cat(" - Event                : ",x$var["event"]," (",txt.eventvar,")\n",sep="")
    cat(" - Time  [min;max]      : ",x$var["time"]," [",signif(attr(x$var,"range.time")[1],3),";",signif(attr(x$var,"range.time")[2], 3),"]\n",sep="")
    print(M.at.risk, row.names = FALSE)

    cat("\n Estimation procedure \n")
    cat(" - Estimator",if(length(estimator.print)>1){"s"}else{" "},"           : ",paste(estimator.print, collapse = " "),"\n",sep="")
    if(test.inference){
        cat(" - Statistical inference: ", txt.inference,sep="")
    }
    
    if(display.results){
        cat("\n Results \n")
        dt.res <- cbind(x$diffRisk[,list(risk.A = formatCI(lower = min(estimate.A),upper = max(estimate.A)),
                                         risk.B = formatCI(lower = min(estimate.B),upper = max(estimate.B)),
                                         "difference (B-A)" = formatCI(lower = min(estimate),upper = max(estimate))),
                                   by=c("A","B")],
                        "ratio (B/A)" = x$ratioRisk[,list(ratio = formatCI(lower = min(estimate),upper = max(estimate))),
                                                    by=c("A","B")]$ratio
                        )
        if(!is.na(x$variable["treatment"])){
            setnames(dt.res, old = c("A","B"), new = paste(x$variable["treatment"], c("A","B"),sep="="))
        }else if(!is.na(x$variable["strata"])){
            setnames(dt.res, old = c("A","B"), new = paste(x$variable["strata"], c("A","B"),sep="="))
        }
        cat(" - Standardized risks   :    [min;max]   \n")
        print(dt.res, row.names = FALSE)
        cat("\n")
            
        cat(" - Computation time     : ",x$computation.time$point," ",attr(x$computation.time$point,"units")," (point estimate)\n", sep="")
        if(!is.na(x$inference$B) && !identical(x$inference$B,0)){
            cat("                          ",x$computation.time$bootstrap," ",attr(x$computation.time$bootstrap,"units")," (bootstrap)\n", sep="")
        }else if(!is.null(x$iid) ||x$band || x$se){
            cat("                          ",x$computation.time$iid," ",attr(x$computation.time$iid,"units")," (iid)\n", sep="")
        }
        cat("\n")
    }
    
    return(NULL)
}


#----------------------------------------------------------------------
### print.ate.R ends here
