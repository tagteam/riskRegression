### calcCIbootATE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 23 2018 (14:07) 
## Version: 
## Last-Updated: maj 23 2018 (14:08) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


calcCIbootATE <- function(boot, meanRisk, riskComparison, type, conf, TD){

    type <- tolower(type) ## convert to lower case
    name.estimate <- names(boot$t0)
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

    test.NA <- !is.na(boot$t)
    test.Inf <- !is.infinite(boot$t)
    n.boot <- colSums(test.NA*test.Inf)

    ## boot estimate
    boot.estimate <- apply(boot$t, 2, mean, na.rm = TRUE)

    ## standard error
    boot.se <- sqrt(apply(boot$t, 2, var, na.rm = TRUE))

    ## confidence interval
    alpha <- 1-conf
    ls.CI <- lapply(index, function(iP){ # iP <- 1
        if(n.boot[iP]==0){
            return(c(lower = NA, upper = NA))
        }else if(type == "Wald"){
            return(c(lower = as.double(boot$t0[iP] + qnorm(alpha/2) * boot.se[iP]),
                     upper = as.double(boot$t0[iP] - qnorm(alpha/2) * boot.se[iP])
                     ))
        }else if(type == "quantile"){
            return(c(lower = as.double(quantile(boot$t[,iP], probs = alpha/2, na.rm = TRUE)),
                     upper = as.double(quantile(boot$t[,iP], probs = 1-(alpha/2), na.rm = TRUE))
                     ))
        }else{
            out <- boot::boot.ci(boot,
                                 conf = conf,
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
            p.value <- boot2pvalue(x = boot$t[,iP],
                                   null = null[iP],
                                   estimate = boot$t0[iP],
                                   alternative = "two.sided",
                                   FUN.ci = function(p.value, sign.estimate, ...){ ## p.value <- 0.4
                                       side.CI <- c(index.lowerCI,index.upperCI)[2-sign.estimate]
                                       boot::boot.ci(boot,
                                                     conf = 1-p.value,
                                                     type = type,
                                                     index = iP)[[slot.boot.ci]][side.CI]
                                   })
            return(p.value)
        }
    })
        
    ## merge
    index.mrisk <- grep("^meanRisk:",name.estimate)
    boot.mrisks <- cbind(meanRisk[,.SD, .SDcols = c("Treatment","time")],
                         meanRiskBoot = boot.estimate[index.mrisk],
                         se = boot.se[index.mrisk],
                         lower = boot.CI[index.mrisk,"lower"],
                         upper = boot.CI[index.mrisk,"upper"],
                         n.boot = n.boot[index.mrisk])

    index.crisk.diff <- grep("^compRisk:diff:",name.estimate)
    index.crisk.ratio <- grep("^compRisk:ratio:",name.estimate)
    boot.crisks <- cbind(riskComparison[, .SD, .SDcols = c("Treatment.A","Treatment.B","time")],
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
                                    
    if (TD){
        key1 <- c("Treatment","landmark")
        key2 <- c("Treatment.A","Treatment.B","landmark")
    }
    else{
        key1 <- c("Treatment","time")
        key2 <- c("Treatment.A","Treatment.B","time")
    }


    ## export
    return(list(meanRisk = merge(meanRisk[,.SD,.SDcols = c("Treatment","time","meanRisk")],
                                 boot.mrisks,by=key1),
                riskComparison = merge(riskComparison[,.SD,.SDcols = c("Treatment.A","Treatment.B","time","diff","ratio")],
                                       boot.crisks,by=key2)
                ))
}

######################################################################
### calcCIbootATE.R ends here
