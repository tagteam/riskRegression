### print.Score.R ---
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: May 31 2016 (11:32)
## Version:
## last-updated: Sep  7 2019 (10:13) 
##           By: Thomas Alexander Gerds
##     Update #: 44
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

##' Print method for risk prediction scores
##'
##' @title Print Score object
##' @param x Object obtained with \code{Score.list}
##' @param digits Number of digits
##' @param ... passed to print
#'
#' @method print Score
#' @export
print.Score <- function(x,digits=3,...){
    for (s in c(x$summary)){
        cat(paste0("\nSummary statistics ",s,":\n"))
        print(x[[s]],digits=digits,response.type=x$response.type,...)
    }
    for (m in c(x$metrics)){
        cat(paste0("\nMetric ",m,":\n"))
        print(x[[m]],digits=digits,response.type=x$response.type,...)
    }
    for (p in c(x$plots)){
        cat(paste0("\nData for ",p," plot are stored in the object as x[[\"",p,"\"]].\n"))
        ## print(x[[m]],digits=digits, ...)
    }
    switch(x$split.method$name,"BootCv"={
        cat("\nBootstrap cross-validation based on ",
            x$split.method$B,
            " ",
            ifelse(x$split.method$N==x$split.method$M,
                   "bootstrap samples (drawn with replacement)",
                   "bootstrap subsamples (drawn without replacement)"),
            " each of size ",
            x$split.method$M,
            ".\n",
            ifelse(x$call$se.fit,paste0("The level of significance is set at ",x$alpha,"\nThe 'confidence intervals' are bootstrap quantiles"),""),
            ifelse(x$call$multi.split.test,"\nThe 'p-values' are median p-values across the splits",""),
            "\n",
            sep="")
    },"LeaveOneOutBoot"={
        cat("\nEfron's leave-one-out-bootstrap based on ",
            x$split.method$B,
            " ",
            ifelse(x$split.method$N==x$split.method$M,
                   "bootstrap samples (drawn with replacement)",
                   "bootstrap subsamples (drawn without replacement)"),
            " each of size ",
            x$split.method$M,
            ".\n",
            "The 'confidence intervals' and 'p-values' are obtained with the delta method after bootstrap.\n",
            sep="")
    })
    if (x$split.method$internal.name == "crossval"){
        cat("\n",x$split.method$name, " repeated ",
            x$split.method$B,
            " times.\n",
            sep="")
    }
}

#' @method print scoreAUC
#' @export
print.scoreAUC <- function(x,B,digits=3,response.type,...){
    cat("\nResults by model:\n\n")
    print(x$score,digits=digits,...)
    if (length(x$contrasts)>0){
        cat("\nResults of model comparisons:\n\n")
        print(x$contrasts,digits=digits,...)
    }
}
#' @method print scoreBrier
#' @export
print.scoreBrier <- function(x,B,digits=3,response.type,...){
    cat("\nResults by model:\n\n")
    print(x$score,digits=digits,...)
    if (length(x$contrasts)>0){
        cat("\nResults of model comparisons:\n\n")
        print(x$contrasts,digits=digits,...)
    }
}

#' @method print scorerisks
#' @export
print.scorerisks <- function(x,B,digits=3,response.type,...){
    if (response.type=="binary")
        data.table::setkeyv(x$score,c("risk"))
    else
        data.table::setkeyv(x$score,c("times","risk"))
    print(x$score,digits=digits)
    if ("times"%in%names(x$score))
        print(data.table::dcast(x$score,times+ID~model,value.var="risk"))
    else
        print(data.table::dcast(x$score,ID~model,value.var="risk"))
}


#----------------------------------------------------------------------
### print.Score.R ends here
