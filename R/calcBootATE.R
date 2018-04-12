### calcBootATE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2018 (17:05) 
## Version: 
## Last-Updated: apr 12 2018 (13:02) 
##           By: Brice Ozenne
##     Update #: 45
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

# {{{ calcBootATE
## generate a boot object for the ate function that will be used to compute CI and p.values
calcBootATE <- function(object, pointEstimate, Gformula, data, 
                        treatment, contrasts, times, cause, landmark, n.contrasts, levels,
                        dots, n.obs,
                        handler, B, seed, mc.cores,
                        verbose){

    name.estimate <- names(pointEstimate)
    
    ## cores
    x.cores <- parallel::detectCores()
    if(mc.cores > x.cores){
        warning("Not enough available cores \n","available: ",parallel::detectCores()," | requested: ",mc.cores,"\n")
        mc.cores <- x.cores
    }

    ## package to be exported to cluster
    if(handler %in% c("snow","foreach") ){
        pp <- find(as.character(object$call[[1]]))
        index.package <- grep("package:",pp)
        addPackage <- if(length(index.package)>0){gsub("package:","",pp[index.package])}else{NULL}
        addPackage <- unique(c("riskRegression","parallel","survival",addPackage))
    }
    
    ## seed
    if (!missing(seed)) set.seed(seed)
    bootseeds <- sample(1:1000000,size=B,replace=FALSE)

    ## bootstrap
    if(handler[[1]] %in% c("snow","parallel")) {
                                        # {{{ use boot package

        if(handler=="snow"){
            ## initialize CPU
            cl <- parallel::makeCluster(mc.cores)
            ## load packages
            parallel::clusterCall(cl, function(x){sapply(x, library, character.only = TRUE)}, addPackage)
            ## set seeds
            parallel::clusterApply(cl, bootseeds, function(x){set.seed(x)})
            ## check
            ## clusterCall(cl, function(x){rnorm(5)})
        }else{
            ## set seeds
            bootseeds <- sum(bootseeds)
            set.seed(bootseeds)
            ##
            cl <- NULL
        }

        ## run bootstrap
        boot.object <- boot::boot(data = data, R = B, statistic = function(data, index, ...){
            dataBoot <- data[index]
            object$call$data <- dataBoot
            objectBoot <- try(eval(object$call),silent=TRUE)
            if ("try-error" %in% class(objectBoot)){
                stop(paste0("Failed to fit model ",class(object)))
            }
            iBoot <- tryCatch(Gformula(object=objectBoot,
                                       data=dataBoot,
                                       treatment=treatment,
                                       contrasts=contrasts,
                                       times=times,
                                       cause=cause,
                                       landmark=landmark,
                                       n.contrasts = n.contrasts,
                                       levels = levels,
                                       dots),
                              error = function(x){return(NULL)})

                if(is.null(iBoot)){ ## error handling
                    out <- setNames(rep(NA, length(name.estimate), name.estimate))
                }else{
                    out <- setNames(c(iBoot$meanRisk$meanRisk,
                                      iBoot$riskComparison$diff,
                                      iBoot$riskComparison$ratio), name.estimate)
                }
                return(out)
            }, sim = "ordinary", stpe = "indices", strata = rep(1, n.obs),
            parallel = handler[[1]], ncpus = mc.cores, cl = cl)

                                        # }}}
    }else {
        if (handler[[1]]=="foreach" && mc.cores>1){
                                        # {{{ foreach
            if(verbose){
                cl <- parallel::makeCluster(mc.cores, outfile = "")
                pb <- txtProgressBar(max = B, style = 3)          
            }else{
                cl <- parallel::makeCluster(mc.cores)
            }
            doParallel::registerDoParallel(cl)
            b <- NULL ## [:forCRANcheck:] foreach
            boots <- foreach::`%dopar%`(foreach::foreach(b = 1:B, .packages = addPackage, .export = NULL), {
                set.seed(bootseeds[[b]])
                if(verbose){setTxtProgressBar(pb, b)}
                dataBoot <- data[sample(1:n.obs, size = n.obs, replace = TRUE),]
                object$call$data <- dataBoot
                objectBoot <- try(eval(object$call),silent=TRUE)
                if ("try-error" %in% class(objectBoot)){
                    stop(paste0("Failed to fit model ",class(object),ifelse(try(b>0,silent=TRUE),paste0(" in bootstrap step ",b,"."))))
                }
                tryCatch(Gformula(object=objectBoot,
                                  data=dataBoot,
                                  treatment=treatment,
                                  contrasts=contrasts,
                                  times=times,
                                  cause=cause,
                                  landmark=landmark,
                                  n.contrasts = n.contrasts,
                                  levels = levels,
                                  dots),
                         error = function(x){return(NULL)})
            })
            if(verbose){close(pb)}
            parallel::stopCluster(cl)
                                        # }}}
        }else{
                                        # {{{ mcapply
            if(Sys.info()["sysname"] == "Windows" && mc.cores>1){
                message("mclapply cannot perform parallel computations on Windows \n",
                        "consider setting argument handler to \"foreach\" \n")
                mc.cores <- 1
            }
            boots <- parallel::mclapply(1:B, function(b){
                set.seed(bootseeds[[b]])
                dataBoot <- data[sample(1:n.obs, size = n.obs, replace = TRUE),]
                object$call$data <- dataBoot
                objectBoot <- try(eval(object$call),silent=TRUE)
                if ("try-error" %in% class(objectBoot)){
                    stop(paste0("Failed to fit model",ifelse(try(b>0,silent=TRUE),paste0(" in bootstrap step ",b,"."))))
                }
                tryCatch(Gformula(object=objectBoot,
                                  data=dataBoot,
                                  treatment=treatment,
                                  contrasts=contrasts,
                                  times=times,
                                  cause=cause,
                                  landmark=landmark,
                                  n.contrasts = n.contrasts,
                                  levels = levels,
                                  dots),
                         error = function(x){return(NULL)})
            }, mc.cores = mc.cores)
                                        # }}}
        }
    
                                        # {{{ convert to boot object
        M.bootEstimate <- do.call(rbind,lapply(boots,function(iBoot){
            c(iBoot$meanRisk$meanRisk, iBoot$riskComparison$diff, iBoot$riskComparison$ratio)
        }))
        if(NROW(M.bootEstimate)==0){
            stop("Error in all bootstrap samples.")
        }
        colnames(M.bootEstimate) <- name.estimate

        boot.object <- list(t0 = pointEstimate,
                            t = M.bootEstimate,
                            R = B,
                            data = data,
                            seed = bootseeds,
                            statistic = NULL,
                            sim = "ordinary",
                            call = quote(boot(data = XX, statistic = XX, R = XX)),
                            stype = "i",
                            strata = rep(1,n.obs),
                            weights = rep(1/n.obs,n.obs),
                            pred.i = NULL,  ## Omitted if m is 0 or sim is not "ordinary" (from doc of boot::boot)
                            L = NULL, ## only used when sim is "antithetic" (from doc of boot::boot)
                            ran.gen = NULL, ## only used when sim is "parametric" (from doc of boot::boot)
                            mle = NULL ## only used when sim is "parametric" (from doc of boot::boot)
                            )
        class(boot.object) <- "boot"
                                        # }}}
    }

    return(list(boot = boot.object,
                bootseeds = bootseeds))
    
}
                                        # }}}

# {{{ calcCIboot
calcCIboot <- function(boot, meanRisk, riskComparison, type, conf, TD){

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
# }}}

######################################################################
### calcBootATE.R ends here
