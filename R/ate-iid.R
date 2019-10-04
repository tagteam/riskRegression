### ate-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: okt  4 2019 (17:15) 
##           By: Brice Ozenne
##     Update #: 843
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
                   data,
                   contrasts,
                   times,
                   cause,
                   iid.ate,
                   iid.outcome,
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
                   eventVar.time,
                   time.jumpC,
                   beforeEvent.jumpC,
                   beforeTau.nJumpC,
                   n.obs,
                   n.times,
                   ...){

    ## ** prepare output
    n.contrasts <- length(contrasts)
    iid.treatment <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.treatment) <- contrasts
    iid.censoring <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.censoring) <- contrasts
    iid.survival <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.survival) <- contrasts              

    ## ** Precompute quantities
    tol <- 1e-12
    if(attr(estimator,"integral")){
        SG <- S.jump*G.jump
        dM_SG <- dM.jump/SG
        dM_SGG <- dM_SG/G.jump
        ls.F1tau_F1t <- lapply(1:n.times, function(iT){-colCenter_cpp(F1.jump, center = F1.tau[,iT])})
        ls.F1tau_F1t_dM_SGG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]*dM_SG/G.jump})
        ls.F1tau_F1t_dM_SSG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]*dM_SG/S.jump})
        ls.F1tau_F1t_SG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]/SG})
    }

    ## ** Compute influence function relative to the prediction
    ## *** outcome model
    if(attr(estimator,"integral")){

        ## **** integral outcome model at tau
        ## compute integral
        int.IFF1_tau <- cbind(0,rowCumSum(dM_SG))
        ## extract integral over the right time spand
        for(iTau in 1:n.times){ ## iTau <- 1
            index.col <- prodlim::sindex(jump.times = c(0,time.jumpC), eval.times = pmin(data[[eventVar.time]],times[iTau]))
            factor <- TRUE
            attr(factor,"factor") <- lapply(1:n.contrasts, function(iC){cbind(iW.IPTW[,iC] * int.IFF1_tau[(1:n.obs) + (index.col-1) * n.obs])})

            ## outputs -IF_S or IF_r : ok
            term.intF1_tau <- attr(predictRisk(object.event, newdata = data, times = times[iTau], cause = cause,
                                     average.iid = factor, diag = FALSE),"average.iid")

            for(iC in 1:n.contrasts){ ## iC <- 1 
                iid.outcome[[iC]][,iTau] <- iid.outcome[[iC]][,iTau] + term.intF1_tau[[iC]]
            }
        }
        
        ## **** integral outcome model at t
        for(iC in 1:n.contrasts){ ## iC <- 1
            factor <- TRUE
            attr(factor, "factor") <- list(-colMultiply_cpp(dM_SG*beforeEvent.jumpC, scale = iW.IPTW[,iC]))                                
            integrand.F1t <- attr(predictRisk(object.event, newdata = data, times = time.jumpC, cause = cause,
                                              average.iid = factor, diag = 2), "average.iid")[[1]] ## argument diag=2 is only used by CauseSpecificCox models
            iid.outcome[[iC]] <- iid.outcome[[iC]] + rowCumSum(integrand.F1t)[,beforeTau.nJumpC]
        }
    }

    ## *** treatment model
    if(attr(estimator,"IPTW")){
        for(iC in 1:n.contrasts){ ## iC <- 1

            factor <- TRUE
            if(estimator == "IPTW"){
                attr(factor,"factor") <- list(colMultiply_cpp(Y.tau, scale = -iW.IPTW2[,iC]))
            }else if(estimator == "IPTW,IPCW"){
                attr(factor,"factor") <- list(colMultiply_cpp(Y.tau * iW.IPCW, scale = -iW.IPTW2[,iC]))
            }else if(estimator == "AIPTW"){
                attr(factor,"factor") <- list(colMultiply_cpp(Y.tau - F1.ctf.tau[[iC]], scale = -iW.IPTW2[,iC]))
            }else if(estimator == "AIPTW,AIPCW"){
                attr(factor,"factor") <- list(colMultiply_cpp(Y.tau * iW.IPCW - F1.ctf.tau[[iC]] + augTerm, scale = -iW.IPTW2[,iC]))
            }

            term.treatment <- attr(predictRisk(object.treatment,
                                               newdata = data,
                                               average.iid = factor,
                                               level = contrasts[iC]), "average.iid")

            iid.treatment[[iC]] <- term.treatment[[1]]
        }
    }
     
    ## *** survival term
    if(attr(estimator,"integral")){
        
        for(iTau in 1:n.times){ ## iTau <- 1
            for(iC in 1:n.contrasts){ ## iC <- 1
                factor <- TRUE
                attr(factor,"factor") <- list(-colMultiply_cpp(ls.F1tau_F1t_dM_SSG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
                integrand.St <- attr(predictRisk(object.event, type = "survival", newdata = data, times = time.jumpC-tol, cause = cause,
                                                 average.iid = factor, diag = 2), "average.iid")[[1]]
                iid.survival[[iC]][,iTau] <- iid.survival[[iC]][,iTau] + rowSums(integrand.St[,1:beforeTau.nJumpC[iTau]])
            }
        }
        
    }
    
    ## *** censoring model
    if(attr(estimator,"IPCW")){

        ## **** IPCW
        for(iTau in 1:n.times){ ## iTau <- 1
            factor <- TRUE
            attr(factor,"factor") <- lapply(1:n.contrasts, function(iC){cbind(-iW.IPTW[,iC]*iW.IPCW2[,iTau]*Y.tau[,iTau])})
            
            term.censoring <- predictCox(object.censor,
                                         newdata = data,
                                         times = pmin(times[iTau], data[[eventVar.time]] - tol),
                                         diag = TRUE,
                                         average.iid = factor)$survival.average.iid

            for(iC in 1:n.contrasts){ ## iC <- 1 
                iid.censoring[[iC]][,iTau] <- term.censoring[[iC]]
            }
        }

        ## **** integral term
        if(attr(estimator,"integral")){
            for(iTau in 1:n.times){ ## iTau <- 1
                for(iC in 1:n.contrasts){ ## iC <- 1
                    ## integral censoring denominator
                    factor <- TRUE
                    attr(factor,"factor") <- list(-colMultiply_cpp(ls.F1tau_F1t_dM_SGG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
                    integrand.G1 <- predictCox(object.censor, newdata = data, times = time.jumpC - tol, 
                                               average.iid = factor)$survival.average.iid[[1]]

                    ## integral censoring martingale
                    factor <- TRUE
                    attr(factor,"factor") <- list(-colMultiply_cpp(ls.F1tau_F1t_SG[[iTau]]*beforeEvent.jumpC, scale = iW.IPTW[,iC]))
                    integrand.G2 <- predictCox(object.censor, newdata = data, times = time.jumpC, type = "hazard",
                                               average.iid = factor)$hazard.average.iid[[1]]

                    ## collect
                    iid.censoring[[iC]][,iTau] <- iid.censoring[[iC]][,iTau] + rowSums(integrand.G1[,1:beforeTau.nJumpC[iTau]] + integrand.G2[,1:beforeTau.nJumpC[iTau]])
                }
            }
        }
    }

    ## ** export
    iid.total <- vector(mode = "list", length = n.contrasts)
    for(iC in 1:n.contrasts){
        iid.total[[iC]] <- iid.ate[[iC]] + iid.outcome[[iC]] + iid.treatment[[iC]] + iid.survival[[iC]] + iid.censoring[[iC]]
    }
    names(iid.total) <- contrasts
    attr(iid.total,"ate") <- iid.ate
    attr(iid.total,"outcome") <- iid.outcome
    attr(iid.total,"treatment") <- iid.treatment
    attr(iid.total,"survival") <- iid.survival
    attr(iid.total,"censoring") <- iid.censoring

    return(iid.total)            
}

## * iidATE2
iidATE2 <- function(estimator,
                    contrasts,
                    iid.ate,
                    iid.outcome,
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
                    ...){

    ## ** prepare output
    n.contrasts <- length(contrasts)
    iid.treatment <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.treatment) <- contrasts
    iid.censoring <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.censoring) <- contrasts
    iid.survival <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.survival) <- contrasts
    
    ## ** Precompute quantities
    tol <- 1e-12
    if(attr(estimator,"integral")){
        SG <- S.jump*G.jump
        dM_SG <- dM.jump/SG   
        ls.F1tau_F1t <- lapply(1:n.times, function(iT){-colCenter_cpp(F1.jump, center = F1.tau[,iT])})
    }
    
    ## ** Compute influence function relative to the predictions

    ## *** outcome model
    if(attr(estimator,"integral")){

        ## **** integral outcome model at tau
        for(iTau in 1:n.times){ ## iTau <- 1
            integral.F1tau <- calcIterm(factor = dM_SG,
                                        iid = iid.nuisance.outcome[,iTau,,drop=FALSE],
                                        indexJump = index.obsSINDEXjumpC[,iTau],
                                        iid.outsideI = TRUE)

            for(iC in 1:n.contrasts){ ## iC <- 1 
                iid.outcome[[iC]][,iTau] <- iid.outcome[[iC]][,iTau] + colMeans(colMultiply_cpp(integral.F1tau, scale = iW.IPTW[,iC]))
            }
        }
        
        ## **** integral outcome model at t
        for(iTau in 1:n.times){ ## iTau <- 1
            integral.F1t <- calcIterm(factor = -dM_SG,
                                      iid = iid.nuisance.outcome[, n.times+(1:n.jumps),,drop=FALSE],
                                      indexJump = index.obsSINDEXjumpC[,iTau],
                                      iid.outsideI = FALSE)

            for(iC in 1:n.contrasts){ ## iC <- 1 
                iid.outcome[[iC]][,iTau] <- iid.outcome[[iC]][,iTau] + colMeans(colMultiply_cpp(integral.F1t, scale = iW.IPTW[,iC]))
            }
        }
    }

    
    ## *** treatment model
    if(attr(estimator,"IPTW")){
        
        for(iC in 1:n.contrasts){ ## iC <- 1
            
            for(iTau in 1:n.times){ ## iTau <- 1
                if(estimator == "IPTW"){
                    iFactor <- - Y.tau[,iTau] * iW.IPTW2[,iC]
                }else if(estimator == "IPTW,IPCW"){
                    iFactor <- - Y.tau[,iTau] * iW.IPCW[,iTau] * iW.IPTW2[,iC]
                }else if(estimator == "AIPTW"){
                    iFactor <- - (Y.tau[,iTau] - F1.ctf.tau[[iC]][,iTau]) * iW.IPTW2[,iC]
                }else if(estimator == "AIPTW,AIPCW"){
                    iFactor <- - (Y.tau[,iTau] * iW.IPCW[,iTau] - F1.ctf.tau[[iC]][,iTau] + augTerm[,iTau]) * iW.IPTW2[,iC]
                }

                iid.treatment[[iC]][,iTau] <- colMeans(colMultiply_cpp(iid.nuisance.treatment[[iC]], scale = iFactor))
            }
        }
    }

    ## *** survival model
    if(attr(estimator,"integral")){

        for(iTau in 1:n.times){ ## iTau <- 1
            integral.Surv <- calcIterm(factor = -ls.F1tau_F1t[[iTau]]*dM_SG/S.jump,
                                       iid = iid.nuisance.survival,
                                       indexJump = index.obsSINDEXjumpC[,iTau],
                                       iid.outsideI = FALSE)

            for(iC in 1:n.contrasts){ ## iC <- 1 
                iid.survival[[iC]][,iTau] <- iid.survival[[iC]][,iTau] + colMeans(colMultiply_cpp(integral.Surv, scale = iW.IPTW[,iC]))
            }
        }
    }
    
    ## *** censoring model
    if(attr(estimator,"IPCW")){

        ## **** IPCW
        for(iTau in 1:n.times){ ## iTau <- 1
            predIID.censoringSurv_obs <- do.call(rbind,lapply(1:n.obs, function(iObs){
                if(index.obsSINDEXjumpC[iObs,iTau]==0){
                    return(rep(0,n.obs))
                }else{
                    return(iid.nuisance.censoring[iObs,index.obsSINDEXjumpC[iObs,iTau],])
                }
            }))
            for(iC in 1:n.contrasts){ ## iC <- 1 
                iid.censoring[[iC]][,iTau] <- colMeans(colMultiply_cpp(predIID.censoringSurv_obs, -iW.IPCW2[,iTau]*Y.tau[,iTau]*iW.IPTW[,iC]))
            }
        }

        ## **** integral term
        if(attr(estimator,"integral")){

            ## set G iid at t-
            iid.nuisance.censoring_tminus <- array(0, dim = c(n.obs,n.jumps,n.obs))
            iid.nuisance.censoring_tminus[,2:n.jumps,] <- iid.nuisance.censoring[,1:(n.jumps-1),]
            
            for(iTau in 1:n.times){ ## iTau <- 1
                ## integral censoring denominator
                integral.G <- calcIterm(factor = -ls.F1tau_F1t[[iTau]]*dM_SG/G.jump,
                                        iid = iid.nuisance.censoring_tminus,
                                        indexJump = index.obsSINDEXjumpC[,iTau],
                                        iid.outsideI = FALSE)

                ## integral censoring martingale
                integral.dLambda <- calcIterm(factor = -ls.F1tau_F1t[[iTau]]/SG,
                                              iid = iid.nuisance.martingale,
                                              indexJump = index.obsSINDEXjumpC[,iTau],
                                              iid.outsideI = FALSE)

                ## collect
                for(iC in 1:n.contrasts){ ## iC <- 1
                    iid.censoring[[iC]][,iTau] <- iid.censoring[[iC]][,iTau] + colMeans(colMultiply_cpp(integral.G + integral.dLambda, scale = iW.IPTW[,iC]))
                }
            }
        }
    }

    ## ** export
    iid.total <- vector(mode = "list", length = n.contrasts)
    for(iC in 1:n.contrasts){
        iid.total[[iC]] <- iid.ate[[iC]] + iid.outcome[[iC]] + iid.treatment[[iC]] + iid.survival[[iC]] + iid.censoring[[iC]]
    }
    names(iid.total) <- contrasts
    attr(iid.total,"ate") <- iid.ate
    attr(iid.total,"outcome") <- iid.outcome
    attr(iid.total,"treatment") <- iid.treatment
    attr(iid.total,"survival") <- iid.survival
    attr(iid.total,"censoring") <- iid.censoring

    return(iid.total)            
}

## * calcIterm (for iidATE2)
calcIterm <- function(factor, iid, indexJump, iid.outsideI){
    n <- length(indexJump)

    if(iid.outsideI){ ## compute iid int factor
        ls.I <- lapply(1:n, function(iObs){ ## iObs <- 2
            iIID <- iid[iObs,1,]
            iFactor <- factor[iObs,1:indexJump[iObs]]
        
            if(indexJump[iObs]==0){
                return(rep(0,n))
            }else {
                return(iIID*sum(iFactor))
            }
        })
    } else { ## compute int iid * factor
        ls.I <- lapply(1:n, function(iObs){ ## iObs <- 2
            iIID <- iid[iObs,1:indexJump[iObs],]
            iFactor <- factor[iObs,1:indexJump[iObs]]
        
            if(indexJump[iObs]==0){
                return(rep(0,n))
            }else if(indexJump[iObs]==1){
                return(iIID*iFactor)
            }else{
                return(rowSums(rowMultiply_cpp(t(iIID), iFactor)))
            }
        })
    }

    return(t(do.call(cbind,ls.I)))
}

######################################################################
### ate-iid.R ends here
