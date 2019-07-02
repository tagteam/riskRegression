### ate-bootstrap.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 11 2018 (17:05) 
## Version: 
## Last-Updated: jul  2 2019 (16:17) 
##           By: Brice Ozenne
##     Update #: 181
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
calcBootATE <- function(args, name.estimate, n.obs, fct.pointEstimate,
                        handler, B, seed, mc.cores, cl,
                        verbose){

                                        # {{{ prepare arguments
    n.estimate <- length(name.estimate)
    
    ## hard copy of the dataset before bootstrap
    ls.data <- list(object.event = NULL,
                    object.treatment = NULL,
                    object.censor = NULL)
    for(iModel in c("object.event","object.treatment","object.censor")){
        ls.data[[iModel]] <- data.table::as.data.table(eval(args[[iModel]]$call$data))
    }

    ## package to be exported to cluster
    vec.fitter <- unique(c(as.character(args$object.event$call[[1]]),
                           as.character(args$object.treatment$call[[1]]),
                           as.character(args$object.censor$call[[1]])))
    ls.package <- lapply(vec.fitter,function(iFitter){
        iSource <- utils::find(iFitter)
        if(grepl("package:",iSource)){gsub("package:","",iSource)}else{NULL}
    })
    add.Package <- unique(c("riskRegression","data.table","parallel","survival",unlist(ls.package)))

    ## if cluster already defined by the user
    no.cl <- is.null(cl)
    if( (no.cl[[1]] == FALSE) && (mc.cores[[1]] == 1) ){ ## i.e. the user has not initialized the number of cores
        mc.cores <- length(cl)
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
    warperBootATE <- function(index, args, fct.pointEstimate, name.estimate, n.estimate){
        ## models for the conditional mean
        for(iModel in c("object.event","object.treatment","object.censor")){
            if(!is.null(args[[iModel]])){
                args[[iModel]]$call$data <- ls.data[[iModel]][index] ## resample dataset
                args[[iModel]] <- try(eval(args[[iModel]]$call),silent=TRUE) ## refit  model
                if ("try-error" %in% class(args[[iModel]])){
                    iBoot <- paste0("Failed to fit model ",iModel," on the bootstrap sample", sep = "")
                    class(iBoot) <- "try-error"
                    return(iBoot)
                }
            }
        }

        ## compute ate
        iBoot <- try(do.call(fct.pointEstimate, args), silent = TRUE)

        ## export
        if(inherits(iBoot,"try-error")){ ## error handling
            out <- setNames(rep(NA, n.estimate), name.estimate)
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
            parallel::clusterCall(cl, function(x){sapply(x, library, character.only = TRUE)}, add.Package)
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
                                      warperBootATE(index = index,
                                                    args = args,
                                                    fct.pointEstimate = fct.pointEstimate,
                                                    name.estimate = name.estimate,
                                                    n.estimate = n.estimate)                                      
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
            boots <- foreach::`%dopar%`(foreach::foreach(b = 1:B, .packages = add.Package, .export = c(".calcLterm","SurvResponseVar")), { ## b <- 1
                if(verbose){setTxtProgressBar(pb, b)}
                set.seed(bootseeds[[b]])
                warperBootATE(index = sample(1:n.obs, size = n.obs, replace = TRUE),
                              args = args,
                              fct.pointEstimate = fct.pointEstimate,
                              name.estimate = name.estimate,
                              n.estimate = n.estimate)                                      
            })            
            if(verbose){close(pb)}
            if(no.cl){parallel::stopCluster(cl)}
                                        # }}}
        }else if(handler=="mclapply"){
                                        # {{{ mclapply
            boots <- parallel::mclapply(1:B, mc.cores = mc.cores, FUN = function(b){
                set.seed(bootseeds[[b]])
                warperBootATE(index = sample(1:n.obs, size = n.obs, replace = TRUE),
                              args = args,
                              fct.pointEstimate = fct.pointEstimate,
                              name.estimate = name.estimate,
                              n.estimate = n.estimate)                                      
            })
                                        # }}}
        }

                                        # {{{ convert to boot object
        M.bootEstimate <- do.call(rbind,boots)
        
        if(all(is.na(M.bootEstimate))){
            stop(paste0("Error in all bootstrap samples: ", attr(boots[[1]],"error")[1]))
        }
        colnames(M.bootEstimate) <- name.estimate
                                        # }}}
    }
                                    
    ## output
    return(list(boot = M.bootEstimate,
                bootseeds = bootseeds))
}

######################################################################
### ate-bootstrap.R ends here
