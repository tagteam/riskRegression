### ate-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: Oct 17 2024 (11:39) 
##           By: Brice Ozenne
##     Update #: 1360
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * iidATE
iidATE <- function(estimator,
                   object.event,
                   object.treatment,
                   object.censor,
                   mydata,
                   contrasts,
                   times,
                   cause,
                   iid.GFORMULA,
                   iid.IPTW,
                   iid.AIPTW,
                   Y.tau,
                   F1.ctf.tau,
                   iW.IPTW,
                   iW.IPTW2,
                   iW.IPCW,
                   iW.IPCW2,
                   augTerm,
                   F1.tau,
                   F1.jump,
                   S.jump,
                   G.jump,
                   dM.jump,
                   index.obsSINDEXjumpC,
                   eventVar.time,
                   time.jumpC,
                   beforeEvent.jumpC,
                   beforeTau.nJumpC,
                   n.obs,
                   n.times,
                   product.limit,
                   store,
                   ...){

    ## ** prepare output
    n.contrasts <- length(contrasts)
    grid <- expand.grid(tau = 1:n.times, contrast = 1:n.contrasts)
    n.grid <- NROW(grid)

    ## ** Precompute quantities
    tol <- 1e-12
    if(attr(estimator,"integral")){
        ls.F1tau_F1t <- lapply(1:n.times, function(iT){-colCenter_cpp(F1.jump, center = F1.tau[,iT])})

        SG <- S.jump*G.jump
        SGG <- SG*G.jump
        SSG <- SG*S.jump

        iSG <- matrix(0, nrow = NROW(SG), ncol = NCOL(SG))
        iSGG <- matrix(0, nrow = NROW(SG), ncol = NCOL(SG))
        iSSG <- matrix(0, nrow = NROW(SG), ncol = NCOL(SG))
        index.beforeEvent.jumpC <- which(beforeEvent.jumpC)
        if(length(index.beforeEvent.jumpC)>0){
            iSG[index.beforeEvent.jumpC] <- 1/SG[index.beforeEvent.jumpC]
            iSGG[index.beforeEvent.jumpC] <- 1/SGG[index.beforeEvent.jumpC]
            iSSG[index.beforeEvent.jumpC] <- 1/SSG[index.beforeEvent.jumpC]
        }
        
        dM_SG <- dM.jump * iSG
        dM_SGG <- dM.jump * iSGG
        dM_SSG <- dM.jump * iSSG
        
        ls.F1tau_F1t_SG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]*iSG})
        ls.F1tau_F1t_dM_SGG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]*dM_SGG})
        ls.F1tau_F1t_dM_SSG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]*dM_SSG})
    }

    test.IPTW <- attr(estimator,"IPTW")
    test.IPCW <- attr(estimator,"IPCW")
    any.IPTW <- "IPTW" %in% attr(estimator,"full")
    any.IPTW.IPCW <- "IPTW,IPCW" %in% attr(estimator,"full")
    any.AIPTW <- "AIPTW" %in% attr(estimator,"full")
    any.AIPTW.AIPCW <- "AIPTW,AIPCW" %in% attr(estimator,"full")
    
    ## ** Compute influence function relative to the prediction
    ## *** outcome model
    ## already done in ATE_TI (ate-pointEstimate.R)

    ## cat("Outcome (method=1) \n")
    ## print(sapply(lapply(iid.AIPTW,abs),colSums))

    ## *** treatment model
    if(test.IPTW){
        for(iC in 1:n.contrasts){ ## iC <- 1
            factor <- TRUE
            attr(factor,"factor") <- list()

            if(any.IPTW.IPCW){
                attr(factor,"factor") <- c(attr(factor,"factor"),
                                           list(IPTW = colMultiply_cpp(Y.tau * iW.IPCW, scale = -iW.IPTW2[,iC]))
                                           )
            }else if(any.IPTW){
                attr(factor,"factor") <- c(attr(factor,"factor"),
                                           list(IPTW = colMultiply_cpp(Y.tau, scale = -iW.IPTW2[,iC]))
                                           )
            }
            
            if(any.AIPTW.AIPCW){
                attr(factor,"factor") <- c(attr(factor,"factor"),
                                           list(AIPTW = colMultiply_cpp(Y.tau * iW.IPCW - F1.ctf.tau[[iC]], scale = -iW.IPTW2[,iC]))
                                           )
            }else if(any.AIPTW){
                attr(factor,"factor") <- c(attr(factor,"factor"),
                                           list(AIPTW = colMultiply_cpp(Y.tau - F1.ctf.tau[[iC]], scale = -iW.IPTW2[,iC]))
                                           )
            }
            term.treatment <- attr(predictRisk(object.treatment,
                                               newdata = mydata,
                                               average.iid = factor,
                                               level = contrasts[iC]), "average.iid")

            if(any.IPTW || any.IPTW.IPCW){
                iid.IPTW[[iC]] <- iid.IPTW[[iC]] + term.treatment[["IPTW"]]
            }
            if(any.AIPTW || any.AIPTW.AIPCW){
                iid.AIPTW[[iC]] <- iid.AIPTW[[iC]] + term.treatment[["AIPTW"]]
            }
        }
    }
    
    ## cat("Treatment (method=1) \n")
    ## print(sapply(lapply(iid.AIPTW,abs),colSums))

    ## *** censoring model
    if(test.IPCW){
        factor <- TRUE        
        for(iTime in 1:n.times){ ## iTime <- 1
            attr(factor, "factor") <- lapply(1:n.contrasts, function(iC){cbind(-iW.IPTW[,iC]*iW.IPCW2[,iTime]*Y.tau[,iTime])})
## browser()
            term.censoring <- attr(predictRisk(object.censor, newdata = mydata, times = c(0,time.jumpC)[index.obsSINDEXjumpC[,iTime]+1],
                                               diag = TRUE, product.limit = product.limit, average.iid = factor, store = store),"average.iid")

            for(iC in 1:n.contrasts){ ## iGrid <- 1
                ## - because predictRisk outputs the risk instead of the survival                 
                if(any.IPTW || any.IPTW.IPCW){
                    iid.IPTW[[iC]][,iTime] <- iid.IPTW[[iC]][,iTime] - term.censoring[[iC]]
                }
                if(any.AIPTW || any.AIPTW.AIPCW){
                    iid.AIPTW[[iC]][,iTime] <- iid.AIPTW[[iC]][,iTime] - term.censoring[[iC]]
                }
            }
        }
        ## colSums(abs(do.call(cbind,term.censoring)))
    }

    ## cat("Censoring (method=1) \n")
    ## print(sapply(lapply(iid.AIPTW,abs),colSums))
    ## *** augmentation censoring term
    if(attr(estimator,"integral")){

        ## **** outcome term
        ## at tau
            
        ## compute integral
        int.IFF1_tau <- cbind(0,rowCumSum(dM_SG))

        ## extract integral over the right time spand
        factor <- TRUE
        attr(factor,"factor") <- lapply(1:n.contrasts, function(iC){
            matrix(NA, nrow = n.obs, ncol = n.times)
        })

        for(iTau in 1:n.times){ ## iTau <- 1
            index.col <- prodlim::sindex(jump.times = c(0,time.jumpC), eval.times = pmin(mydata[[eventVar.time]],times[iTau]))
            for(iC in 1:n.contrasts){
                attr(factor,"factor")[[iC]][,iTau] <- cbind(iW.IPTW[,iC] * int.IFF1_tau[(1:n.obs) + (index.col-1) * n.obs])
            }
        }

        term.intF1_tau <- attr(predictRisk(object.event, newdata = mydata, times = times, cause = cause,
                                           average.iid = factor, product.limit = product.limit, store = store),"average.iid")

        for(iC in 1:n.contrasts){ ## iC <- 1
            iid.AIPTW[[iC]] <- iid.AIPTW[[iC]] + term.intF1_tau[[iC]]
        }
          
        ## ## at t
        factor <- TRUE
        attr(factor, "factor") <- lapply(1:n.contrasts, function(iC){
            -colMultiply_cpp(dM_SG*beforeEvent.jumpC, scale = iW.IPTW[,iC])
        })
        
        integrand.F1t <- attr(predictRisk(object.event, newdata = mydata, times = time.jumpC, cause = cause,
                                          average.iid = factor, product.limit = product.limit, store = store), "average.iid")

        for(iC in 1:n.contrasts){ ## iC <- 1
            iid.AIPTW[[iC]] <- iid.AIPTW[[iC]] + subsetIndex(rowCumSum(integrand.F1t[[iC]]),
                                                             index = beforeTau.nJumpC,
                                                             default = 0, col = TRUE)
        }
        ## cat("Augmentation outcome (method=1) \n")
        ## print(sapply(lapply(iid.AIPTW,abs),colSums))

        ## **** treatment term
        for(iC in 1:n.contrasts){ ## iC <- 1
            factor <- TRUE
            attr(factor,"factor") <- list(AIPTW = colMultiply_cpp(augTerm, scale = -iW.IPTW2[,iC]))

            term.treatment <- attr(predictRisk(object.treatment,
                                               newdata = mydata,
                                               average.iid = factor,
                                               level = contrasts[iC]), "average.iid")

            ## print(colSums(term.treatment[["AIPTW"]]^2))
            iid.AIPTW[[iC]] <- iid.AIPTW[[iC]] + term.treatment[["AIPTW"]]
        }
        ## cat("Augmentation treatment (method=1) \n")
        ## print(sapply(lapply(iid.AIPTW,abs),colSums))

        ## **** survival term
        factor <- TRUE
        attr(factor,"factor") <- lapply(1:n.grid, function(iGrid){ ## iGrid <- 1
            iTau <- grid[iGrid,"tau"]
            iC <- grid[iGrid,"contrast"]
            return(-colMultiply_cpp(ls.F1tau_F1t_dM_SSG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
        })
        integrand.St <- attr(predictRisk(object.event, type = "survival", newdata = mydata, times = time.jumpC-tol, cause = cause,
                                         average.iid = factor, product.limit = product.limit, store = store), "average.iid")
        for(iGrid in 1:n.grid){ ## iGrid <- 1
            iTau <- grid[iGrid,"tau"]
            iC <- grid[iGrid,"contrast"]
                        
            if(beforeTau.nJumpC[iTau]>0){
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowSums(integrand.St[[iGrid]][,1:beforeTau.nJumpC[iTau],drop=FALSE])
            }
        }
        ## cat("Augmentation survival (method=1) \n")
        ## print(sapply(lapply(iid.AIPTW,abs),colSums))
        
        ## **** censoring term            
        ## ## integral censoring denominator
        factor <- TRUE
        attr(factor,"factor") <- lapply(1:n.grid, function(iGrid){ ## iGrid <- 1
            iTau <- grid[iGrid,"tau"]
            iC <- grid[iGrid,"contrast"]
            return(-colMultiply_cpp(ls.F1tau_F1t_dM_SGG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
        })

        integrand.G1 <- predictCox(object.censor, newdata = mydata, times = time.jumpC - tol, 
                                   average.iid = factor, store = store)$survival.average.iid

        ## ## integral censoring martingale
        factor <- TRUE
        attr(factor,"factor") <- lapply(1:n.grid, function(iGrid){ ## iGrid <- 1
            iTau <- grid[iGrid,"tau"]
            iC <- grid[iGrid,"contrast"]
            return(-colMultiply_cpp(ls.F1tau_F1t_SG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
        })
        integrand.G2 <- predictCox(object.censor, newdata = mydata, times = time.jumpC, type = "hazard",
                                   average.iid = factor, store = store)$hazard.average.iid

        for(iGrid in 1:n.grid){ ## iGrid <- 1
            iTau <- grid[iGrid,"tau"]
            iC <- grid[iGrid,"contrast"]
            if(beforeTau.nJumpC[iTau]>0){
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowSums(integrand.G1[[iGrid]][,1:beforeTau.nJumpC[iTau],drop=FALSE] + integrand.G2[[iGrid]][,1:beforeTau.nJumpC[iTau],drop=FALSE])
            }
        }
    } ## end attr(estimator,"integral")
    ## cat("Augmentation (method=1) \n")
    ## print(sapply(lapply(iid.AIPTW,abs),colSums))
    
    ## ** export
    out <- list()

    if(attr(estimator,"export.GFORMULA")){
        out <- c(out, list(GFORMULA = iid.GFORMULA))
    }
    if(attr(estimator,"export.IPTW")){
        out <- c(out, list(IPTW = iid.IPTW))
    }
    if(attr(estimator,"export.AIPTW")){
        out <- c(out, list(AIPTW = iid.AIPTW))
    }
    return(out)            
}

## * iidATE2 (same as iidATE but perform the average from the iid, less efficient but easier to understand)
iidATE2 <- function(estimator,
                    contrasts,
                    iid.GFORMULA,
                    iid.IPTW,
                    iid.AIPTW,
                    Y.tau,
                    F1.ctf.tau,
                    iW.IPTW,
                    iW.IPTW2,
                    iW.IPCW,
                    iW.IPCW2,
                    augTerm,
                    F1.tau,
                    F1.jump,
                    S.jump,
                    G.jump,
                    dM.jump,
                    iid.nuisance.outcome,
                    iid.nuisance.treatment,
                    iid.nuisance.censoring,
                    iid.nuisance.censoring.diag,
                    iid.nuisance.survival,
                    iid.nuisance.martingale,
                    index.obsSINDEXjumpC,
                    index.obsSINDEXjumpC.int,
                    n.obs,
                    n.times,
                    n.jumps,
                    method.iid,
                    ...){

    ## ** prepare output
    n.contrasts <- length(contrasts)
    
    ## ** Precompute quantities
    tol <- 1e-12
    if(attr(estimator,"integral")){
        SG <- S.jump*G.jump
        dM_SG <- dM.jump/SG   
        ls.F1tau_F1t <- lapply(1:n.times, function(iT){-colCenter_cpp(F1.jump, center = F1.tau[,iT])})
    }
    
    test.IPTW <- attr(estimator,"IPTW")
    test.IPCW <- attr(estimator,"IPCW")
    any.IPTW <- "IPTW" %in% attr(estimator,"full")
    any.IPTW.IPCW <- "IPTW,IPCW" %in% attr(estimator,"full")
    any.AIPTW <- "AIPTW" %in% attr(estimator,"full")
    any.AIPTW.AIPCW <- "AIPTW,AIPCW" %in% attr(estimator,"full")

    ## ** Compute influence function relative to the predictions

    ## *** outcome model
    ## already done in ATE_TI (ate-pointEstimate.R)

    ## cat("Outcome (method=2) \n")
    ## print(sapply(lapply(iid.AIPTW,abs),colSums))

    ## *** treatment model
    if(test.IPTW){
        
        for(iC in 1:n.contrasts){ ## iC <- 1
            for(iTau in 1:n.times){ ## iTau <- 1

                if(any.IPTW || any.IPTW.IPCW){
                    if(any.IPTW){
                        iFactor <- - Y.tau[,iTau] * iW.IPTW2[,iC]
                    }else if(any.IPTW.IPCW){
                        iFactor <- - Y.tau[,iTau] * iW.IPCW[,iTau] * iW.IPTW2[,iC]
                    }
                    iid.IPTW[[iC]][,iTau] <- iid.IPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(iid.nuisance.treatment[[iC]], scale = iFactor))
                }
                if(any.AIPTW || any.AIPTW.AIPCW){
                    if(any.AIPTW){
                        iFactor <- - (Y.tau[,iTau] - F1.ctf.tau[[iC]][,iTau]) * iW.IPTW2[,iC]
                    }else if(any.AIPTW.AIPCW){
                        iFactor <- - (Y.tau[,iTau] * iW.IPCW[,iTau] - F1.ctf.tau[[iC]][,iTau]) * iW.IPTW2[,iC]
                    }
                    iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(iid.nuisance.treatment[[iC]], scale = iFactor))
                }
                
            }
        }
    }
    ## cat("Treatment (method=2) \n")
    ## print(sapply(lapply(iid.AIPTW,abs),colSums))
    
    ## *** censoring model
    if(test.IPCW){

        for(iTau in 1:n.times){ ## iTau <- 1
            for(iC in 1:n.contrasts){ ## iC <- 1
                if(any.IPTW || any.IPTW.IPCW){
                    iid.IPTW[[iC]][,iTau] <- iid.IPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(iid.nuisance.censoring.diag[[iTau]][,1,], -iW.IPCW2[,iTau]*Y.tau[,iTau]*iW.IPTW[,iC]))
                }
                if(any.AIPTW || any.AIPTW.AIPCW){
                    iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(iid.nuisance.censoring.diag[[iTau]][,1,], -iW.IPCW2[,iTau]*Y.tau[,iTau]*iW.IPTW[,iC]))
                }
            }
        }
    }
    
    ## cat("Censoring (method=2) \n")
    ## print(sapply(lapply(iid.AIPTW,abs),colSums))

    ## *** augmentation censoring term
    if(attr(estimator,"integral")){

        
        ## **** outcome model
        for(iTau in 1:n.times){ ## iTau <- 1
            ## at tau
            integral.F1tau <- calcItermOut(factor = dM_SG,
                                           iid = iid.nuisance.outcome[,iTau,,drop=FALSE],
                                           indexJump = index.obsSINDEXjumpC.int[,iTau],
                                           n = n.obs)

            ## at t
            integral.F1t <- calcItermIn(factor = -dM_SG,
                                        iid = iid.nuisance.outcome[, n.times+(1:n.jumps),,drop=FALSE],
                                        indexJump = index.obsSINDEXjumpC.int[,iTau],
                                        n = n.obs)

            ## assemble
            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(integral.F1tau + integral.F1t, scale = iW.IPTW[,iC]))
            }

        }

        ## cat("Augmentation outcome (method=2) \n")
        ## print(sapply(lapply(iid.AIPTW,abs),colSums))

        ## **** treatment model
        for(iTau in 1:n.times){ ## iTau <- 1
            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(iid.nuisance.treatment[[iC]],
                                                                                            scale = - augTerm[,iTau] * iW.IPTW2[,iC]))
            }
        }
        ## cat("Augmentation treatment (method=2) \n")
        ## print(sapply(lapply(iid.AIPTW,abs),colSums))

        ## **** survival model
        for(iTau in 1:n.times){ ## iTau <- 1
            integral.Surv <- calcItermIn(factor = -ls.F1tau_F1t[[iTau]]*dM_SG/S.jump,
                                         iid = iid.nuisance.survival,
                                         indexJump = index.obsSINDEXjumpC.int[,iTau],
                                         n = n.obs)

            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(integral.Surv, scale = iW.IPTW[,iC]))
            }
        }
        ## cat("Augmentation survival (method=2) \n")
        ## print(sapply(lapply(iid.AIPTW,abs),colSums))
            
        ## **** censoring term
        for(iTau in 1:n.times){ ## iTau <- 1
            ## integral censoring denominator
            integral.G <- calcItermIn(factor = -ls.F1tau_F1t[[iTau]]*dM_SG/G.jump,
                                      iid = iid.nuisance.censoring,
                                      indexJump = index.obsSINDEXjumpC.int[,iTau],
                                      n = n.obs)

            ## integral censoring martingale
            integral.dLambda <- calcItermIn(factor = -ls.F1tau_F1t[[iTau]]/SG,
                                            iid = iid.nuisance.martingale,
                                            indexJump = index.obsSINDEXjumpC.int[,iTau],
                                            n = n.obs)

            ## collect
            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(integral.G + integral.dLambda, scale = iW.IPTW[,iC]))
            }
        }
    } ## end attr(estimator,"integral")
    ## cat("Augmentation (method=2) \n")
    ## print(sapply(lapply(iid.AIPTW,abs),colSums))

    ## ** export
    out <- list()
    if(attr(estimator,"export.GFORMULA")){
        out <- c(out, list(GFORMULA = iid.GFORMULA))
    }
    if(attr(estimator,"export.IPTW")){
        out <- c(out, list(IPTW = iid.IPTW))
    }
    if(attr(estimator,"export.AIPTW")){
        out <- c(out, list(AIPTW = iid.AIPTW))
    }
    
    return(out)            
}

## * calcIterm (for iidATE2)
calcItermOut <- function(factor, iid, indexJump, n){ ## compute iid int factor
    ls.I <- lapply(1:n, function(iObs){ ## iObs <- 2
        iIID <- iid[,1,iObs]
        iFactor <- factor[iObs,1:indexJump[iObs]]
        
        if(indexJump[iObs]==0){
            return(rep(0,n))
        }else {
            return(iIID*sum(iFactor))
        }
    })
    return(do.call(cbind,ls.I))
}

calcItermIn <- function(factor, iid, indexJump, n){ ## compute int iid * factor

    ls.I <- lapply(1:n, function(iObs){ ## iObs <- 1
        iIID <- iid[,1:indexJump[iObs],iObs]
        iFactor <- factor[iObs,1:indexJump[iObs]]
        
        if(indexJump[iObs]==0){
            return(rep(0,n))
        }else if(indexJump[iObs]==1){
            return(iIID*iFactor)
        }else{
            return(rowSums(rowMultiply_cpp(iIID, iFactor)))
        }
    })
    return(do.call(cbind,ls.I))
}

######################################################################
### ate-iid.R ends here
