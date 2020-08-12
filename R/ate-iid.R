### ate-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: aug 12 2020 (11:42) 
##           By: Brice Ozenne
##     Update #: 1177
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
                   iid.Gformula,
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
                   dLambda.jump,
                   eventVar.time,
                   time.jumpC,
                   beforeEvent.jumpC,
                   beforeTau.nJumpC,
                   n.obs,
                   n.times,
                   product.limit,
                   ...){

    ## ** prepare output
    n.contrasts <- length(contrasts)
    grid <- expand.grid(tau = 1:n.times, contrast = 1:n.contrasts)
    n.grid <- NROW(grid)
    
    ## ** Precompute quantities
    tol <- 1e-12
    if(attr(estimator,"integral")){
        SG <- S.jump*G.jump
        ls.F1tau_F1t <- lapply(1:n.times, function(iT){-colCenter_cpp(F1.jump, center = F1.tau[,iT])})
        dM_SG <- dM.jump/SG
        dM_SGG <- dM_SG/G.jump
        ls.F1tau_F1t_dM_SGG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]*dM_SG/G.jump})
        ls.F1tau_F1t_dM_SSG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]*dM_SG/S.jump})
        ls.F1tau_F1t_SG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]/SG})
    }

    ## ** Compute influence function relative to the prediction
    ## *** outcome model
    ## already done in ATE_TI (ate-pointEstimate.R)

    ## cat("Outcome (method=1) \n")
    ## print(head(iid.AIPTW[[1]]))

    ## ## *** treatment model
    if(attr(estimator,"IPTW")){
        for(iC in 1:n.contrasts){ ## iC <- 1

            factor <- TRUE
            attr(factor,"factor") <- list()

            if("IPTW,IPCW" %in% estimator){
                attr(factor,"factor") <- c(attr(factor,"factor"),
                                           list(IPTW = colMultiply_cpp(Y.tau * iW.IPCW, scale = -iW.IPTW2[,iC]))
                                           )
            }else if("IPTW" %in% estimator){
                attr(factor,"factor") <- c(attr(factor,"factor"),
                                           list(IPTW = colMultiply_cpp(Y.tau, scale = -iW.IPTW2[,iC]))
                                           )
            }
            
            if("AIPTW,AIPCW" %in% estimator){
                attr(factor,"factor") <- c(attr(factor,"factor"),
                                           list(AIPTW = colMultiply_cpp(Y.tau * iW.IPCW - F1.ctf.tau[[iC]], scale = -iW.IPTW2[,iC]))
                                           )
            }else if("AIPTW" %in% estimator){
                attr(factor,"factor") <- c(attr(factor,"factor"),
                                           list(AIPTW = colMultiply_cpp(Y.tau - F1.ctf.tau[[iC]], scale = -iW.IPTW2[,iC]))
                                           )
            }

            term.treatment <- attr(predictRisk(object.treatment,
                                               newdata = mydata,
                                               average.iid = factor,
                                               level = contrasts[iC]), "average.iid")

            if(attr(estimator,"export.IPTW")){
                iid.IPTW[[iC]] <- iid.IPTW[[iC]] + term.treatment[["IPTW"]]
            }
            if(attr(estimator,"export.AIPTW")){
                iid.AIPTW[[iC]] <- iid.AIPTW[[iC]] + term.treatment[["AIPTW"]]
            }
        }
    }
    ## cat("Treatment (method=1) \n")
    ## print(head(iid.AIPTW[[1]]))

    ## *** censoring model
    if(attr(estimator,"IPCW")){

        factor <- TRUE
        attr(factor,"factor") <- lapply(1:n.grid, function(iGrid){
            iTau <- grid[iGrid,"tau"]
            iC <- grid[iGrid,"contrast"]
            return(cbind(-iW.IPTW[,iC]*iW.IPCW2[,iTau]*Y.tau[,iTau]))
        })
        
        term.censoring <- predictCox(object.censor,
                                     newdata = mydata,
                                     times = mydata[[eventVar.time]] - tol, ## same as pmin(times[iTau], mydata[[eventVar.time]] - tol) because indicator in the weights
                                     diag = TRUE,
                                     average.iid = factor)$survival.average.iid

        for(iGrid in 1:n.grid){ ## iGrid <- 1
            iTau <- grid[iGrid,"tau"]
            iC <- grid[iGrid,"contrast"]
            if(attr(estimator,"export.IPTW")){
                iid.IPTW[[iC]][,iTau] <- iid.IPTW[[iC]][,iTau] + term.censoring[[iGrid]]
            }
            if(attr(estimator,"export.AIPTW")){
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + term.censoring[[iGrid]]
            }
        }
    }

    ## cat("Censoring (method=1) \n")
    ## print(head(iid.AIPTW[[1]]))

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
                                           average.iid = factor, product.limit = product.limit),"average.iid")
  
        for(iC in 1:n.contrasts){ ## iC <- 1
            iid.AIPTW[[iC]] <- iid.AIPTW[[iC]] + term.intF1_tau[[iC]]
        }
          
        ## ## at t
        factor <- TRUE
        attr(factor, "factor") <- lapply(1:n.contrasts, function(iC){
            -colMultiply_cpp(dM_SG*beforeEvent.jumpC, scale = iW.IPTW[,iC])
        })
        
        integrand.F1t <- attr(predictRisk(object.event, newdata = mydata, times = time.jumpC, cause = cause,
                                          average.iid = factor, product.limit = product.limit), "average.iid")

        for(iC in 1:n.contrasts){ ## iC <- 1
            iid.AIPTW[[iC]] <- iid.AIPTW[[iC]] + subsetIndex(rowCumSum(integrand.F1t[[iC]]),
                                                             index = beforeTau.nJumpC,
                                                             default = 0, col = TRUE)
        }
        ## cat("Augmentation outcome (method=1) \n")
        ## print(head(iid.AIPTW[[1]]))

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
        ## print(head(iid.AIPTW[[1]]))

            ## **** survival term
            factor <- TRUE
            attr(factor,"factor") <- lapply(1:n.grid, function(iGrid){ ## iGrid <- 1
                iTau <- grid[iGrid,"tau"]
                iC <- grid[iGrid,"contrast"]
                return(-colMultiply_cpp(ls.F1tau_F1t_dM_SSG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
            })
            integrand.St <- attr(predictRisk(object.event, type = "survival", newdata = mydata, times = time.jumpC-tol, cause = cause,
                                             average.iid = factor, product.limit = product.limit), "average.iid")
        
            for(iGrid in 1:n.grid){ ## iGrid <- 1
                iTau <- grid[iGrid,"tau"]
                iC <- grid[iGrid,"contrast"]
                        
                if(beforeTau.nJumpC[iTau]>0){
                    iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowSums(integrand.St[[iGrid]][,1:beforeTau.nJumpC[iTau],drop=FALSE])
                }
            }
        ## cat("Augmentation survival (method=1) \n")
        ## print(head(iid.AIPTW[[1]]))

            ## **** censoring term            
            ## ## integral censoring denominator
            factor <- TRUE
            attr(factor,"factor") <- lapply(1:n.grid, function(iGrid){ ## iGrid <- 1
                iTau <- grid[iGrid,"tau"]
                iC <- grid[iGrid,"contrast"]
                return(-colMultiply_cpp(ls.F1tau_F1t_dM_SGG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
            })

            integrand.G1 <- predictCox(object.censor, newdata = mydata, times = time.jumpC - tol, 
                                       average.iid = factor)$survival.average.iid

            ## ## integral censoring martingale
            factor <- TRUE
            attr(factor,"factor") <- lapply(1:n.grid, function(iGrid){ ## iGrid <- 1
                iTau <- grid[iGrid,"tau"]
                iC <- grid[iGrid,"contrast"]
                return(-colMultiply_cpp(ls.F1tau_F1t_SG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
            })
            integrand.G2 <- predictCox(object.censor, newdata = mydata, times = time.jumpC, type = "hazard",
                                       average.iid = factor)$hazard.average.iid

            for(iGrid in 1:n.grid){ ## iGrid <- 1
                iTau <- grid[iGrid,"tau"]
                iC <- grid[iGrid,"contrast"]
                if(beforeTau.nJumpC[iTau]>0){
                    iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowSums(integrand.G1[[iGrid]][,1:beforeTau.nJumpC[iTau],drop=FALSE] + integrand.G2[[iGrid]][,1:beforeTau.nJumpC[iTau],drop=FALSE])
                }
            }
    } ## end attr(estimator,"integral")
    ## cat("Augmentation (method=1) \n")
    ## print(head(iid.AIPTW[[1]]))

    ## ** export
    out <- list()
    if(attr(estimator,"export.Gformula")){
        out <- c(out, list(Gformula = iid.Gformula))
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
                    iid.Gformula,
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
                    iid.nuisance.survival,
                    iid.nuisance.martingale,
                    index.obsSINDEXjumpC,
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
    
    ## ** Compute influence function relative to the predictions

    ## *** outcome model
    ## already done in ATE_TI (ate-pointEstimate.R)

    ## cat("Outcome (method=2) \n")
    ## print(head(iid.AIPTW[[1]]))

    ## *** treatment model
    if(attr(estimator,"IPTW")){
        
        for(iC in 1:n.contrasts){ ## iC <- 1
            for(iTau in 1:n.times){ ## iTau <- 1

                if(attr(estimator,"export.IPTW")){
                    if("IPTW" %in% estimator){
                        iFactor <- - Y.tau[,iTau] * iW.IPTW2[,iC]
                    }else if("IPTW,IPCW" %in% estimator){
                        iFactor <- - Y.tau[,iTau] * iW.IPCW[,iTau] * iW.IPTW2[,iC]
                    }
                    iid.IPTW[[iC]][,iTau] <- iid.IPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(iid.nuisance.treatment[[iC]], scale = iFactor))
                }
                if(attr(estimator,"export.AIPTW")){
                    if("AIPTW" %in% estimator){
                        iFactor <- - (Y.tau[,iTau] - F1.ctf.tau[[iC]][,iTau]) * iW.IPTW2[,iC]
                    }else if("AIPTW,AIPCW" %in% estimator){
                        iFactor <- - (Y.tau[,iTau] * iW.IPCW[,iTau] - F1.ctf.tau[[iC]][,iTau]) * iW.IPTW2[,iC]
                    }
                    iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(iid.nuisance.treatment[[iC]], scale = iFactor))
                }
                
            }
        }
    }
    ## cat("Treatment (method=2) \n")
    ## print(head(iid.AIPTW[[1]]))
    
    ## *** censoring model
    if(attr(estimator,"IPCW")){

        for(iTau in 1:n.times){ ## iTau <- 1
            predIID.censoringSurv_obs <- do.call(cbind,lapply(1:n.obs, function(iObs){
                if(index.obsSINDEXjumpC[iObs,iTau]==0){
                    return(rep(0,n.obs))
                }else{
                    return(iid.nuisance.censoring[,index.obsSINDEXjumpC[iObs,iTau],iObs])
                }
            }))
            for(iC in 1:n.contrasts){ ## iC <- 1
                if(attr(estimator,"export.IPTW")){
                    iid.IPTW[[iC]][,iTau] <- iid.IPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(predIID.censoringSurv_obs, -iW.IPCW2[,iTau]*Y.tau[,iTau]*iW.IPTW[,iC]))
                }
                if(attr(estimator,"export.AIPTW")){
                    iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(predIID.censoringSurv_obs, -iW.IPCW2[,iTau]*Y.tau[,iTau]*iW.IPTW[,iC]))
                }
            }
        }
    }
    
    ## cat("Censoring (method=2) \n")
    ## print(head(iid.AIPTW[[1]]))

    ## *** augmentation censoring term
    if(attr(estimator,"integral")){

        ## **** outcome model
        for(iTau in 1:n.times){ ## iTau <- 1
            ## at tau
            integral.F1tau <- calcItermOut(factor = dM_SG,
                                           iid = iid.nuisance.outcome[,iTau,,drop=FALSE],
                                           indexJump = index.obsSINDEXjumpC[,iTau],
                                           n = n.obs)

            ## at t
            integral.F1t <- calcItermIn(factor = -dM_SG,
                                         iid = iid.nuisance.outcome[, n.times+(1:n.jumps),,drop=FALSE],
                                         indexJump = index.obsSINDEXjumpC[,iTau],
                                         n = n.obs)

            ## assemble
            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(integral.F1tau + integral.F1t, scale = iW.IPTW[,iC]))
            }

        }

        ## cat("Augmentation outcome (method=2) \n")
        ## print(head(iid.AIPTW[[1]]))

        ## **** treatment model
        for(iTau in 1:n.times){ ## iTau <- 1
            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(iid.nuisance.treatment[[iC]],
                                                                                            scale = - augTerm[,iTau] * iW.IPTW2[,iC]))
            }
        }
        ## cat("Augmentation treatment (method=2) \n")
        ## print(head(iid.AIPTW[[1]]))

        ## **** survival model
        for(iTau in 1:n.times){ ## iTau <- 1
            integral.Surv <- calcItermIn(factor = -ls.F1tau_F1t[[iTau]]*dM_SG/S.jump,
                                         iid = iid.nuisance.survival,
                                         indexJump = index.obsSINDEXjumpC[,iTau],
                                         n = n.obs)

            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(integral.Surv, scale = iW.IPTW[,iC]))
            }
        }
        ## cat("Augmentation survival (method=2) \n")
        ## print(head(iid.AIPTW[[1]]))
            
        ## **** censoring term
        ## set G iid at t-
        iid.nuisance.censoring_tminus <- array(0, dim = c(n.obs,n.jumps,n.obs))
        iid.nuisance.censoring_tminus[,2:n.jumps,] <- iid.nuisance.censoring[,1:(n.jumps-1),]
            
        for(iTau in 1:n.times){ ## iTau <- 1
            ## integral censoring denominator
            integral.G <- calcItermIn(factor = -ls.F1tau_F1t[[iTau]]*dM_SG/G.jump,
                                      iid = iid.nuisance.censoring_tminus,
                                      indexJump = index.obsSINDEXjumpC[,iTau],
                                      n = n.obs)

            ## integral censoring martingale
            integral.dLambda <- calcItermIn(factor = -ls.F1tau_F1t[[iTau]]/SG,
                                            iid = iid.nuisance.martingale,
                                            indexJump = index.obsSINDEXjumpC[,iTau],
                                            n = n.obs)

            ## collect
            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.AIPTW[[iC]][,iTau] <- iid.AIPTW[[iC]][,iTau] + rowMeans(rowMultiply_cpp(integral.G + integral.dLambda, scale = iW.IPTW[,iC]))
            }
        }
    } ## end attr(estimator,"integral")
    ## cat("Augmentation (method=2) \n")
    ## print(head(iid.AIPTW[[1]]))

    ## ** export
    out <- list()
    if(attr(estimator,"export.Gformula")){
        out <- c(out, list(Gformula = iid.Gformula))
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
