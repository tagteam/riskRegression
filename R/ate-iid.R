### ate-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: sep 27 2019 (11:22) 
##           By: Brice Ozenne
##     Update #: 676
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
                   treatment,
                   strata,
                   contrasts,
                   levels,
                   times,
                   cause,
                   n.censor,
                   level.censoring,
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
                   store.jumps,
                   n.obs,
                   n.times,
                   n.jumps,
                   ...){

    n.contrasts <- length(contrasts)
    tol <- 1e-12
    if(inherits(object.event,"CauseSpecificCox")){
        n.cause <- length(object.event$causes)
    }

    iid.treatment <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.treatment) <- contrasts
    iid.censoring <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.censoring) <- contrasts
    iid.survival <- lapply(1:n.contrasts,function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
    names(iid.survival) <- contrasts              

    ## ** Precompute quantities
    if(attr(estimator,"integral")){
        SG <- S.jump*G.jump[,1:n.jumps]
        dM_SG <- dM.jump/SG
        dM_SGG <- dM_SG/G.jump[,1:n.jumps]
        ls.F1tau_F1t <- lapply(1:n.times, function(iT){-colCenter_cpp(F1.jump, center = F1.tau[,iT])})
        ls.F1tau_F1t_dM_SG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]*dM_SG})
        ls.F1tau_F1t_dM_SGG <- lapply(1:n.times, function(iT){ls.F1tau_F1t_dM_SG[[iT]]/G.jump[,1:n.jumps]})
        ls.F1tau_F1t_dM_SSG <- lapply(1:n.times, function(iT){ls.F1tau_F1t_dM_SG[[iT]]/S.jump})
        ls.F1tau_F1t_SG <- lapply(1:n.times, function(iT){ls.F1tau_F1t[[iT]]/SG})
    }
    
    
    ## ** Compute influence function relative to each cox model
    vec.IF.times <- times
    if(attr(estimator,"integral")){
        vec.IF.times <- c(vec.IF.times,time.jumpC)
    }
    vec.IF.times <- sort(unique(vec.IF.times))

    test.Cox <- inherits(object.event,"coxph") || inherits(object.event,"cph") || inherits(object.event,"phreg")
    test.CSC <- inherits(object.event,"CauseSpecificCox")
    test.glm <- inherits(object.event,"glm")

    if(attr(estimator,"Gformula") && (test.Cox || test.CSC) && identical(is.iidCox(object.event),FALSE)){
        object.event <- iidCox(object.event, tau.max = max(times), return.object = TRUE)
    }

    if(attr(estimator,"IPCW")){
        if(identical(is.iidCox(object.censor$iid),FALSE)){
            object.censor$iid <- iidCox(object.censor, tau.max = max(times))
        }
    }
    

    ## ** Compute influence function relative to the prediction
    ## *** outcome model (computation of Prob[T<=t,Delta=1|A,W] = F_1(t|A=a,W))
    if(attr(estimator,"Gformula") && attr(estimator,"integral")){

        ## **** integral outcome model at tau
        ## compute integral
        int.IFF1_tau <- cbind(0,rowCumSum(dM_SG))
        ## extract integral over the right time spand
        for(iTau in 1:n.times){ ## iTau <- 1
            index.col <- prodlim::sindex(jump.times = c(0,time.jumpC), eval.times = pmin(data[[eventVar.time]],times[iTau]))
            factor <- colMultiply_cpp(iW.IPTW, scale = int.IFF1_tau[(1:n.obs) + (index.col-1) * n.obs])
            term.intF1_tau <- attr(predictRiskIID(object.event, newdata = data, times = times[iTau], cause = cause,
                                                  average.iid = TRUE, factor = factor, diag = FALSE),"iid")

            for(iC in 1:n.contrasts){ ## iC <- 1 
                iid.outcome[[iC]][,iTau] <- iid.outcome[[iC]][,iTau] + term.intF1_tau[[iC]]
            }
        }
        
        ## **** integral outcome model at t
        for(iC in 1:n.contrasts){ ## iC <- 1
            factor <- -colMultiply_cpp(dM_SG, scale = iW.IPTW[,iC])
            term.intF1_t <- attr(predictRiskIID(object.event, newdata = data, times = time.jumpC, cause = cause,
                                                average.iid = TRUE, factor = factor, diag = 2), "iid")[[1]]

            iid.outcome[[iC]] <- iid.outcome[[iC]] + calcAugmentation_cpp(term = term.intF1_t,
                                                                          index = store.jumps-1,
                                                                          nObs = n.obs,
                                                                          nTau = n.times)
        }
    }

    ## *** treatment model
    if(attr(estimator,"IPTW")){
        for(iC in 1:n.contrasts){ ## iC <- 1

            if(estimator == "IPTW"){
                factor <- colMultiply_cpp(Y.tau, scale = -iW.IPTW2[,iC])
            }else if(estimator == "IPTW,IPCW"){
                factor <- colMultiply_cpp(Y.tau * iW.IPCW, scale = -iW.IPTW2[,iC])
            }else if(estimator == "AIPTW"){
                factor <- colMultiply_cpp(Y.tau - F1.ctf.tau[[iC]], scale = -iW.IPTW2[,iC])
            }else if(estimator == "AIPTW,AIPCW"){
                factor <- colMultiply_cpp(Y.tau * iW.IPCW - F1.ctf.tau[[iC]] + augTerm, scale = -iW.IPTW2[,iC])
            }
            
            term.treatment <- attr(predictRiskIID(object.treatment,
                                                  newdata = data,
                                                  average.iid = TRUE,
                                                  factor = factor,
                                                  level = contrasts[iC]), "iid")

            iid.treatment[[iC]] <- do.call(cbind,term.treatment)
        }
    }
     
    ## *** survival term
    if(attr(estimator,"integral")){
        for(iTau in 1:n.times){ ## iTau <- 1
            for(iC in 1:n.contrasts){ ## iC <- 1
                factor <- -colMultiply_cpp(ls.F1tau_F1t_dM_SSG[[iTau]], scale = iW.IPTW[,iC])
                term.intS_t <- attr(predictRiskIID(object.event, type = "survival", newdata = data, times = time.jumpC, cause = cause,
                                                   average.iid = TRUE, factor = factor, diag = 2), "iid")[[1]]

                iid.survival[[iC]][,iTau] <- calcAugmentation_cpp(term = term.intS_t,
                                                                  index = store.jumps[,iTau,drop=FALSE]-1,
                                                                  nObs = n.obs,
                                                                  nTau = 1)
            }
        }
    }
    
    ## *** censoring/survival model
    if(attr(estimator,"IPCW")){

        ## **** IPCW
        for(iTau in 1:n.times){ ## iTau <- 1
            factor <- colMultiply_cpp(iW.IPTW, scale =  -iW.IPCW2[,iTau]*Y.tau[,iTau])
            term.censoring <- predictRiskIID(object.censor,
                                             newdata = data,
                                             times = pmin(times[iTau], data[[eventVar.time]] - tol),
                                             diag = TRUE,
                                             average.iid = TRUE,
                                             factor = factor)
            for(iC in 1:n.contrasts){ ## iC <- 1 
                iid.censoring[[iC]][,iTau] <- term.censoring[[iC]]
            }
        }

        ## **** integral term
        if(attr(estimator,"integral")){
            for(iTau in 1:n.times){ ## iTau <- 1
                for(iC in 1:n.contrasts){ ## iC <- 1
                    ## integral censoring denominator
                    factor <- -colMultiply_cpp(ls.F1tau_F1t_dM_SGG[[iTau]], scale = iW.IPTW[,iC])
                    term.intG1 <- attr(predictRiskIID(object.censor, newdata = data, times = time.jumpC, cause = cause,
                                                      average.iid = TRUE, factor = factor, diag = 2), "iid")[[1]]

                    ## integral censoring martingale
                    factor <- TRUE
                    attr(factor,"factor") <- list(-colMultiply_cpp(ls.F1tau_F1t_SG[[iTau]], scale = iW.IPTW[,iC]))
                    term.intG2 <- predictCox(object.censor, newdata = data, times = time.jumpC, type = "hazard",
                                             average.iid = factor, diag = FALSE)$hazard.average.iid[[1]]

                    ## collect
                    iid.censoring[[iC]][,iTau] <- iid.censoring[[iC]][,iTau] + calcAugmentation_cpp(term = term.intG1 + term.intG2,
                                                                                                    index = store.jumps[,iTau,drop=FALSE]-1,
                                                                                                    nObs = n.obs,
                                                                                                    nTau = 1)

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
    
    return(iid.total)            
}

######################################################################
### ate-iid.R ends here
