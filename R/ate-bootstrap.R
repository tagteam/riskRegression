### ate-bootstrap.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2018 (17:05) 
## Version: 
## Last-Updated: jun 27 2019 (12:08) 
##           By: Brice Ozenne
##     Update #: 144
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calcBootATE
## generate a boot object for the ate function that will be used to compute CI and p.values
calcBootATE <- function(object, pointEstimate, fct.pointEstimate, data, formula, TD,
                        treatment, contrasts, times, cause, landmark, n.contrasts, levels,
                        dots, n.obs,
                        handler, B, seed, mc.cores, cl,
                        verbose){

                                        # {{{ prepare arguments
    name.estimate <- names(pointEstimate)
    
    no.cl <- is.null(cl)
    if( (no.cl[[1]] == FALSE) && (mc.cores[[1]] == 1) ){ ## i.e. the user has not initialized the number of cores
        mc.cores <- length(cl)
    }
    ## package to be exported to cluster
    if(handler %in% c("snow","foreach") ){
        pp <- find(as.character(object$call[[1]]))
        index.package <- grep("package:",pp)
        addPackage <- if(length(index.package)>0){gsub("package:","",pp[index.package])}else{NULL}
        addPackage <- unique(c("riskRegression","data.table","parallel","survival",addPackage))
    }
    ## seed
    if (!missing(seed)){
        set.seed(seed)
    }
    bootseeds <- sample(1:1000000,size=B,replace=FALSE)

    ## allArgs <- c("warperBootATE","data","n.obs","fct.pointEstimate",
                 ## "object","treatment","contrasts","times","cause","landmark",
                 ## "n.contrasts","levels","TD","name.estimate","formula","dots")

                                        # }}}

                                        # {{{ warper
    warperBootATE <- function(dataBoot, fct.pointEstimate,
                              object, treatment, contrasts, times, cause, landmark, n.contrasts, levels, TD, name.estimate, formula, dots){
        ## update dataset
        object$call$data <- dataBoot

        ## refit models for the conditional mean
        objectBoot <- try(eval(object$call),silent=TRUE)
        if ("try-error" %in% class(objectBoot)){
            iBoot <- paste0("Failed to fit model ",class(object))
            class(iBoot) <- "try-error"
            return(iBoot)
        }
        ## gather information
        args.bootstrap <- list(object=objectBoot,
                               data=dataBoot,
                               treatment=treatment,
                               contrasts=contrasts,
                               times=times,
                               cause=cause,
                               landmark=landmark,
                               n.contrasts = n.contrasts,
                               levels = levels,
                               dots)
        if (TD){
            args.bootstrap <- c(args.bootstrap, list(formula=formula))
        }
        ## compute ate
        iBoot <- try(do.call(fct.pointEstimate, args.bootstrap), silent = TRUE)

        ## export
        if(inherits(iBoot,"try-error")){ ## error handling
            out <- setNames(rep(NA, length(name.estimate)), name.estimate)
            attr(out,"error") <- iBoot
            return(out)        
        }else{
            return(setNames(c(iBoot$meanRisk$meanRisk,iBoot$riskComparison$diff,iBoot$riskComparison$ratio), name.estimate))
        }
    }

                                        # }}}
    
    ## bootstrap
    if(handler %in% c("snow","multicore")) {
                                        # {{{ use boot package
        if(handler=="snow" && no.cl[[1]]==TRUE){ 
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
        boot.object <- boot::boot(data = data,
                                  R = B,
                                  sim = "ordinary",
                                  stpe = "indices",
                                  strata = rep(1, n.obs),
                                  parallel = handler,
                                  ncpus = mc.cores,
                                  cl = cl,
                                  statistic = function(data, index, ...){
                                      warperBootATE(data = data[index], fct.pointEstimate = fct.pointEstimate,
                                                    object = object, treatment = treatment, contrasts = contrasts, times = times, cause = cause, landmark = landmark,
                                                    n.contrasts = n.contrasts, levels = levels, TD = TD, name.estimate = name.estimate, formula = formula, dots = dots)                                      
                                  })
                                        # }}}
    }else{

        if (handler=="foreach"){ # [@TAG: removed mc.cores > 1]
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
            boots <- foreach::`%dopar%`(foreach::foreach(b = 1:B, .packages = addPackage), { ## b <- 1
                if(verbose){setTxtProgressBar(pb, b)}
                set.seed(bootseeds[[b]])
                warperBootATE(dataBoot = data[sample(1:n.obs, size = n.obs, replace = TRUE)], fct.pointEstimate = fct.pointEstimate,
                              object = object, treatment = treatment, contrasts = contrasts, times = times, cause = cause, landmark = landmark,
                              n.contrasts = n.contrasts, levels = levels, TD = TD, name.estimate = name.estimate, formula = formula, dots = dots)
            })            
            if(verbose){close(pb)}
            if(no.cl){parallel::stopCluster(cl)}
                                        # }}}
        }else if(handler=="mclapply"){
                                        # {{{ mclapply
            boots <- parallel::mclapply(1:B, mc.cores = mc.cores, FUN = function(b){
                set.seed(bootseeds[[b]])
                warperBootATE(dataBoot = data[sample(1:n.obs, size = n.obs, replace = TRUE)], fct.pointEstimate = fct.pointEstimate,
                              object = object, treatment = treatment, contrasts = contrasts, times = times, cause = cause, landmark = landmark,
                              n.contrasts = n.contrasts, levels = levels, TD = TD, name.estimate = name.estimate, formula = formula, dots = dots)
            })
                                        # }}}
        }

                                        # {{{ convert to boot object
        M.bootEstimate <- do.call(rbind,boots)
        
        if(all(is.na(M.bootEstimate))){
            stop(paste0("Error in all bootstrap samples: ", attr(boots[[1]],"error")[1]))
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
                                        
    ## output
    return(list(boot = boot.object,
                bootseeds = bootseeds))
}

######################################################################
### ate-bootstrap.R ends here
