### calcBootATE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2018 (17:05) 
## Version: 
## Last-Updated: Oct  7 2018 (14:59) 
##           By: Thomas Alexander Gerds
##     Update #: 70
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
calcBootATE <- function(object, pointEstimate, Gformula, data, formula, TD,
                        treatment, contrasts, times, cause, landmark, n.contrasts, levels,
                        dots, n.obs,
                        handler, B, seed, mc.cores, cl,
                        verbose){
    name.estimate <- names(pointEstimate)
    no.cl <- is.null(cl)
    if( (no.cl == FALSE) && (mc.cores == 1) ){ ## i.e. the user has not initialized the number of cores
        mc.cores <- length(cl)
    }
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
        if(handler=="snow" && no.cl){
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
        }
        ## run bootstrap
        boot.object <- boot::boot(data = data, R = B, statistic = function(data, index, ...){
            dataBoot <- data[index]
            object$call$data <- dataBoot
            objectBoot <- try(eval(object$call),silent=TRUE)
            if ("try-error" %in% class(objectBoot)){
                stop(paste0("Failed to fit model ",class(object)))
            }
            Gargs <- list(object=objectBoot,data=dataBoot,treatment=treatment,contrasts=contrasts,times=times,cause=cause,landmark=landmark,n.contrasts = n.contrasts,levels = levels,dots)
            if (TD) Gargs <- c(Gargs,list(formula=formula))
            iBoot <- tryCatch(do.call(Gformula, Gargs),
                              error = function(x){return(NULL)})
            if(is.null(iBoot)){ ## error handling
                out <- setNames(rep(NA, length(name.estimate), name.estimate))
            }else{
                out <- setNames(c(iBoot$meanRisk$meanRisk,iBoot$riskComparison$diff,iBoot$riskComparison$ratio), name.estimate)
            }
            return(out)
        }, sim = "ordinary", stpe = "indices", strata = rep(1, n.obs),
        parallel = handler[[1]], ncpus = mc.cores, cl = cl)
        # }}}
    }else {
        if (handler[[1]]=="foreach" && mc.cores>1){
            # {{{ foreach
            if(no.cl){
                if(verbose){
                    cl <- parallel::makeCluster(mc.cores, outfile = "")
                }else{
                    cl <- parallel::makeCluster(mc.cores)
                }                
            }           
            doParallel::registerDoParallel(cl)
            ## progress bar 
            if(verbose){pb <- txtProgressBar(max = B, style = 3)}
            b <- NULL ## [:forCRANcheck:] foreach
            boots <- foreach::`%dopar%`(foreach::foreach(b = 1:B, .packages = addPackage, .export = c("formula")), { ## b <- 1
                set.seed(bootseeds[[b]])
                if(verbose){setTxtProgressBar(pb, b)}
                dataBoot <- data[sample(1:n.obs, size = n.obs, replace = TRUE),]
                object$call$data <- dataBoot
                objectBoot <- try(eval(object$call),silent=TRUE)
                if ("try-error" %in% class(objectBoot)){
                    stop(paste0("Failed to fit model ",class(object),ifelse(try(b>0,silent=TRUE),paste0(" in bootstrap step ",b,"."))))
                }
                Gargs <- list(object=objectBoot,data=dataBoot,treatment=treatment,contrasts=contrasts,times=times,cause=cause,landmark=landmark,n.contrasts = n.contrasts,levels = levels,dots)
                if (TD) Gargs <- c(Gargs,list(formula=formula))
                iBoot <- tryCatch(do.call(Gformula, Gargs),
                                  error = function(x){return(NULL)})
                ## error = function(x){return(x)})
            })
            if(verbose){close(pb)}
            if(no.cl){parallel::stopCluster(cl)}
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

                Gargs <- list(object=objectBoot,data=dataBoot,treatment=treatment,contrasts=contrasts,times=times,cause=cause,landmark=landmark,n.contrasts = n.contrasts,levels = levels,dots)
                if (TD) Gargs <- c(Gargs,list(formula=formula))
                iBoot <- tryCatch(do.call(Gformula, Gargs),
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


######################################################################
### calcBootATE.R ends here
