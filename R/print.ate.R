### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: okt  7 2021 (10:18) 
##           By: Brice Ozenne
##     Update #: 552
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
#' @param estimator [character] The type of estimator relative to which the risks should be output. 
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

    estimator.print <- sapply(attr(estimator,"full"),
                              switch,
                              "GFORMULA" = "G-formula",
                              "GFORMULATD" = "G-formula with time-dependent covariates)",
                              "IPTW,IPCW" = "Inverse probability of treatment weighting",
                              "IPTW" = "Inverse probability of treatment weighting",
                              "AIPTW,AIPCW" = "Augmented estimator",
                              "AIPTW" = "Augmented estimator")


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
    txt.eventvar <- paste0("cause: ",x$theCause)
    if(length(setdiff(cause,x$theCause))>0){
        txt.eventvar <- paste0(txt.eventvar,", competing risk(s): ",paste(setdiff(cause,x$theCause),collapse = " "))
    }
    if(any(attr(x$eval.times,"n.censored")>0)){
        txt.eventvar <- paste0(txt.eventvar, ", censoring: ",attr(x$theCause,"level.censoring"))
    }
    ## number at risk at the evaluation time
    if(!is.na(x$variable["time"])){
        M.at.risk <- cbind("- Eval. time       "="number at risk",attr(x$eval.times,"n.at.risk"))
        colnames(M.at.risk)[2] <- " :"
        M.at.risk[,2] <- paste(as.character(M.at.risk[[2]]),"  ",sep="")
    }
    
    ## statistical inference
    if(x$inference$bootstrap){
        bootci.method <- switch(x$inference$bootci.method,
                                "norm" = "Normal",
                                "basic" = "Basic",
                                "stud" = "Studentized",
                                "perc" = "Percentile",
                                "bca" = "BCa",
                                "wald" = "Wald",
                                "quantile" = "Percentile")
        txt.inference <- paste(bootci.method," bootstrap based on ",x$inference$n.bootstrap," bootstrap samples\n",
                               "                that were drawn with replacement from the original data.\n",sep="")
        test.inference <- TRUE
    }else if(any(x$inference[,c("se","band","ci","p.value")]>0)){
        txt.inference <- "Gaussian approximation \n                where the variance is estimated via the influence function \n"
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
    if(!is.na(x$var["time"])){
        cat(" - Time  [min;max]      : ",x$var["time"]," [",signif(attr(x$var,"range.time")[1],3),";",signif(attr(x$var,"range.time")[2], 3),"]\n",sep="")
        print(M.at.risk, row.names = FALSE)        
    }
    
    cat("\n Estimation procedure \n")
    cat(" - Estimator",if(length(estimator.print)>1){"s"}else{" "}," : ",paste(estimator.print, collapse = " "),"\n",sep="")
    if(test.inference){
        cat(" - Uncertainty: ", txt.inference,sep="")
    }
    
    if(display.results){
        cat("\n Results \n")
        if(all(x$diffRisk[,.N,by=c("A","B")]$N==1)){
            test.onerow <- TRUE
            dt.res <- cbind(x$diffRisk[,list(risk.A = .SD$estimate.A,
                                             risk.B = .SD$estimate.B,
                                             "difference (B-A)" = .SD$estimate),
                                       by=c("A","B")],
                            "ratio (B/A)" = x$ratioRisk[,list(ratio = .SD$estimate),
                                                        by=c("A","B")]$ratio
                            )
        }else{
            test.onerow <- FALSE
            dt.res <- cbind(x$diffRisk[,list(risk.A = Publish::formatCI(lower = min(.SD$estimate.A),upper = max(.SD$estimate.A)),
                                             risk.B = Publish::formatCI(lower = min(.SD$estimate.B),upper = max(.SD$estimate.B)),
                                             "difference (B-A)" = Publish::formatCI(lower = min(.SD$estimate),upper = max(.SD$estimate))),
                                       by=c("A","B")],
                            "ratio (B/A)" = x$ratioRisk[,list(ratio = Publish::formatCI(lower = min(.SD$estimate),upper = max(.SD$estimate))),
                                                        by=c("A","B")]$ratio
                            )
        }
        if(!is.na(x$variable["treatment"])){
            setnames(dt.res, old = c("A","B"), new = paste(x$variable["treatment"], c("A","B"),sep="="))
        }else if(!is.na(x$variable["strata"])){
            setnames(dt.res, old = c("A","B"), new = paste(x$variable["strata"], c("A","B"),sep="="))
        }
        if(test.onerow){
            cat(" - Standardized risks   :                \n")
        }else{
            cat(" - Standardized risks   :    [min;max]   \n")
        }        
        print(dt.res, row.names = FALSE)

        cat("\n")
        
        cat(" - Computation time     : ",x$computation.time$point," ",attr(x$computation.time$point,"units")," (point estimate)\n", sep="")
        if(x$inference$bootstrap){
            cat("                          ",x$computation.time$bootstrap," ",attr(x$computation.time$bootstrap,"units")," (bootstrap)\n", sep="")            
        }else if(!is.null(x$iid) ||x$inference$band || x$inference$se){
            cat("                          ",x$computation.time$iid," ",attr(x$computation.time$iid,"units")," (iid)\n", sep="")
        }
        cat("\n")
    }
    
    return(invisible(NULL))
}


#----------------------------------------------------------------------
### print.ate.R ends here
