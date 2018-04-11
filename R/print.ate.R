### print.ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jun  6 2016 (06:48) 
## Version: 
## last-updated: apr 11 2018 (20:50) 
##           By: Brice Ozenne
##     Update #: 90
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Print average treatment effects
#'
#' Print average treatment effects
#' @param x object obtained with function \code{ate}
#' @param type Character. Method for constructing bootstrap confidence intervals.
#' Either "perc" (the default), "norm", "basic", "stud", or "bca".
#' Argument passed to \code{boot::boot.ci}.
#' @param digits Number of digits
#' @param print should something be displayed in the console?
#' @param ... passed to print
#'
#' @details When using bootstrap resampling the p-values are computing using a test-inversion method,
#' i.e. find the critical confidence level such that one side of the confidence interval overlap the null hypothesis.
#' The p-value is 1 minus the critical confidence level.
#' 
#' @method print ate
#' @export
print.ate <- function(x, type = "perc", digits = 3, print = TRUE, ...){

                                        # {{{ compute confidence intervals and p-values when using the bootstrap
    if(!is.null(x$boot)){
        name.estimate <- names(x$boot$t0)
        n.estimate <- length(name.estimate)
        index <- 1:n.estimate

        slot.boot.ci <- switch(type,
                               "norm" = "normal",
                               "basic" = "basic",
                               "stud" = "student",
                               "perc" = "percent",
                               "bca" = "bca")
        index.lowerCI <- switch(type,
                                "norm" = 2,
                                "basic" = 4,
                                "stud" = 4,
                                "perc" = 4,
                                "bca" = 4)
        index.upperCI <- switch(type,
                                "norm" = 3,
                                "basic" = 5,
                                "stud" = 5,
                                "perc" = 5,
                                "bca" = 5)

        test.NA <- !is.na(x$boot$t)
        test.Inf <- !is.infinite(x$boot$t)
        n.boot <- colSums(test.NA*test.Inf)

        ## boot estimate
        boot.estimate <- apply(x$boot$t, 2, mean, na.rm = TRUE)

        ## standard error
        boot.se <- sqrt(apply(x$boot$t, 2, var, na.rm = TRUE))

        ## confidence interval
        ls.CI <- lapply(index, function(iP){ # iP <- 9
            if(n.boot[iP]==0){
                return(c(lower = NA, upper = NA))
            }else{
                out <- boot::boot.ci(x$boot,
                                     conf = x$conf,
                                     type = type,
                                     index = iP)[[slot.boot.ci]][index.lowerCI:index.upperCI]
                return(setNames(out,c("lower","upper")))
            }
        })
        boot.CI <- do.call(rbind,ls.CI)
        
        ## pvalue
        null <- setNames(rep(0,length(name.estimate)),name.estimate)
        null[grep("^compRisk:ratio", name.estimate)] <- 1

        boot.p <- sapply(index, function(iP){ # iP <- 1
            if(n.boot[iP]==0){
                return(NA)
            }else{
                ## search confidence level such that quantile of CI which is close to 0
                p.value <- boot2pvalue(x = x$boot$t[,iP],
                                       null = null[iP],
                                       estimate = x$boot$t0[iP],
                                       alternative = "two.sided",
                                       FUN.ci = function(p.value, sign.estimate, ...){ ## p.value <- 0.4
                                           side.CI <- c(index.lowerCI,index.upperCI)[2-sign.estimate]
                                           boot::boot.ci(x$boot,
                                                         conf = 1-p.value,
                                                         type = type,
                                                         index = iP)[[slot.boot.ci]][side.CI]
                                       })
                return(p.value)
            }
        })
        
        ## export
            
        index.mrisk <- grep("^meanRisk:",name.estimate)
        boot.mrisks <- cbind(x$meanRisk[,.(Treatment,time)],
                             meanRiskBoot = boot.estimate[index.mrisk],
                             se = boot.se[index.mrisk],
                             lower = boot.CI[index.mrisk,"lower"],
                             upper = boot.CI[index.mrisk,"upper"],
                             n.boot = n.boot[index.mrisk])

        index.crisk.diff <- grep("^compRisk:diff:",name.estimate)
        index.crisk.ratio <- grep("^compRisk:ratio:",name.estimate)
        boot.crisks <- cbind(x$riskComparison[,.(Treatment.A,Treatment.B,time)],
                             diffMeanBoot = boot.estimate[index.crisk.diff],
                             diff.se = boot.se[index.crisk.diff],
                             diff.lower = boot.CI[index.crisk.diff,"lower"],
                             diff.upper = boot.CI[index.crisk.diff,"upper"],
                             diff.p.value = boot.p[index.crisk.diff],
                             diff.n.boot = n.boot[index.crisk.diff],
                             ratioMeanBoot = boot.estimate[index.crisk.ratio],
                             ratio.se = boot.se[index.crisk.ratio],
                             ratio.lower = boot.CI[index.crisk.ratio,"lower"],
                             ratio.upper = boot.CI[index.crisk.ratio,"upper"],
                             ratio.p.value = boot.p[index.crisk.ratio],
                             ratio.n.boot = n.boot[index.crisk.ratio]
                             )
                                    
        if (x$TD){
            key1 <- c("Treatment","landmark")
            key2 <- c("Treatment.A","Treatment.B","landmark")
        }
        else{
            key1 <- c("Treatment","time")
            key2 <- c("Treatment.A","Treatment.B","time")
        }

        
        x$meanRisk <- merge(x$meanRisk[,.(Treatment,time,meanRisk)],
                            boot.mrisks,by=key1)
        x$riskComparison <- merge(x$riskComparison[,.(Treatment.A,Treatment.B,time,diff,ratio)],
                                  boot.crisks,by=key2)            

    }
                                        # }}}

                                        # {{{ display
    if(print){
        cat("The treatment variable ",x$treatment," has the following options:\n",sep="")
        cat(paste(x$contrasts,collapse=", "),"\n")
        cat("\nMean risks on probability scale [0,1] in hypothetical worlds\nin which all subjects are treated with one of the treatment options:\n\n")
        print(x$meanRisk,digits=digits,...)
        cat("\nComparison of risks on probability scale [0,1] between\nhypothetical worlds are interpretated as if the treatment was randomized:\n\n")    
        print(x$riskComparison,digits=digits,...)
        ##
        if(x$se && (x$conf.level > 0 && (x$conf.level < 1))){
            if(x$B==0){
                cat("\nWald confidence intervals are based on asymptotic standard errors.",sep="")
            }else {
                type <- switch(type,
                               "norm" = "Normal",
                               "basic" = "Basic",
                               "stud" = "Studentized",
                               "perc" = "Percentile",
                               "bca" = "BCa")
                cat("\n",type," bootstrap confidence intervals based on ",x$B," bootstrap samples\nthat were drawn with replacement from the original data.",sep="")
            }
            if(x$nsim.band>0){
                cat("\nConfidence bands are based on ",x$nsim.band," simulations",sep="")
            }
            cat("\nConfidence level:",x$conf.level,"\n")
        }
    }
                                        # }}}
    
    ## export
    return(invisible(x))
}


#----------------------------------------------------------------------
### print.ate.R ends here
