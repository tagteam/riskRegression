### ate-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: sep  6 2019 (16:56) 
##           By: Brice Ozenne
##     Update #: 517
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
                   index.store,
                   time.jumpC,
                   n.obs,
                   n.times,
                   n.jumps,
                   ...){

    n.contrasts <- length(contrasts)
    tol <- 1e-12

    ## ** Compute influence function relative to the prediction

    ## *** outcome model (computation of Prob[T<=t,Delta=1|A,W] = F_1(t|A=a,W))
    if(attr(estimator,"Gformula")){
        
        for(iC in 1:n.contrasts){ ## iC <- 1 

            if(!is.null(treatment)){
                ## hypothetical world: in which every subject is treated with the same treatment
                data.i <- data.table::copy(data)
                data.i[[treatment]] <- factor(contrasts[iC], levels = levels)
            }else{
                ## hypothetical world: only patients with the same strata variable exist
                index.strata <- which(data[[strata]]==contrasts[iC])
                data.i <- data[index.strata]
                
            }

            if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                factor <- cbind(1-iW.IPTW[,iC])
            }else{
                factor <- matrix(1, nrow = NROW(data.i), ncol = 1)
            }
            iid.ate[[iC]] <- iid.ate[[iC]] + predictRiskIID(object.event, newdata = data.i, times = times,
                                                                      average.iid = TRUE, factor = factor, cause = cause)[[1]]

        }
        

        if(attr(estimator,"integral")){
            F1.all.iid <- predictRiskIID(object.event, newdata = data, times = c(times, time.jumpC), cause = cause,
                                        average.iid = FALSE)
            F1.tau.iid <- F1.all.iid[,1:n.times,,drop=FALSE]
            F1.jump.iid <- F1.all.iid[,n.times + (1:n.jumps),,drop=FALSE]
        }

    }

    ## *** treatment model
    if(attr(estimator,"IPTW")){
        for(iC in 1:n.contrasts){ ## iC <- 1

            if(estimator %in% c("IPTW","IPTW,IPCW")){
                factor <- colMultiply_cpp(Y.tau * iW.IPCW, scale = -iW.IPTW2[,iC])
            }else if(estimator %in% c("AIPTW")){
                factor <- colMultiply_cpp(Y.tau * iW.IPCW - F1.ctf.tau[[iC]], scale = -iW.IPTW2[,iC])
            }else if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                factor <- colMultiply_cpp(Y.tau * iW.IPCW - F1.ctf.tau[[iC]] + augTerm, scale = -iW.IPTW2[,iC])
            }
            
            term.treatment <- predictRiskIID(object.treatment,
                                             newdata = data,
                                             average.iid = TRUE,
                                             factor = factor,
                                             level = contrasts[iC])
            
            iid.ate[[iC]] <- iid.ate[[iC]] + do.call(cbind,term.treatment)
        }
    }

    ## *** survival model
    if(attr(estimator,"integral")){
        if(inherits(object.event,"CauseSpecificCox")){ ## competing risk case
            S.jump.iid <- predict(object.event, type = "survival", newdata = data, times = time.jumpC, product.limit = FALSE, iid = TRUE)$survival.iid
        }else{ ## survival case
            S.jump.iid <- - F1.jump.iid
        }
    }

    ## *** censoring model
    if(attr(estimator,"IPCW")){

        ## at all times of jump of the censoring process
        G.jump.iid <- - predictRiskIID(object.censor, newdata = data, times = c(0,time.jumpC), average.iid = FALSE)

        ## at observation times
        G.T_tau.iid <- array(NA, dim = c(n.obs, n.times, n.obs))
        for(iTau in 1:n.times){  # iTau <- 1
            ## select the jump corresponding to the event time of each observation
            index.col <- prodlim::sindex(jump.times = c(0,time.jumpC),
                                         eval.times = pmin(data[[eventVar.time]] - tol, times[iTau]))
            
            for(iCol in index.col){ ## iCol <- index.col[1]
                index.n <- which(index.col == iCol)
                G.T_tau.iid[index.n,iTau,] <- G.jump.iid[index.n,iCol,]
            }
            ## check
            ## range(G.T_tau.iid[,iTau,,drop=FALSE] - predictCox(object.censor, newdata = data, times = pmin(times[iTau], data[[eventVar.time]] - tol), type = "survival", diag = TRUE, iid = TRUE)$survival.iid)
            for(iC in 1:n.contrasts){ ## iC <- 1
                iid.ate[[iC]][,iTau] <- iid.ate[[iC]][,iTau] - colMeans(colMultiply_cpp(G.T_tau.iid[,iTau,], scale = iW.IPTW[,iC]*iW.IPCW2[,iTau]*Y.tau[,iTau]))
            }
        }
        
        if(attr(estimator,"integral")){

            ## martingale for the censoring process
            dLambda.jump.iid <- predictCox(object.censor, newdata = data, times = time.jumpC, type = "hazard", iid = TRUE)$hazard.iid

            ## check
            ## Lambda.iid_after <- predictCox(object.censor, newdata = data, times = time.jumpC, type = "cumhazard", iid = TRUE)$cumhazard.iid
            ## Lambda.iid_before <- predictCox(object.censor, newdata = data, times = time.jumpC - tol, type = "cumhazard", iid = TRUE)$cumhazard.iid
            ## range(dLambda.jump.iid - (Lambda.iid_after - Lambda.iid_before))
        }
    }

    ## ** iid relative to the augmentation term
    if(attr(estimator,"integral")){
        for(iN in 1:n.obs){ ## iN <- 1
                term.augm <- .calcAugmentationTerm_iid(F1 = F1.jump, F1.tau = F1.tau, S = S.jump, G = G.jump[,1:n.jumps], dM = dM.jump,
                                                       F1.iid = matrix(F1.jump.iid[iN,,,drop=FALSE], nrow = n.obs, ncol = n.jumps, byrow = TRUE),
                                                       F1.tau.iid = matrix(F1.tau.iid[iN,,,drop=FALSE], nrow = n.obs, ncol = n.times, byrow = TRUE),
                                                       S.iid = matrix(S.jump.iid[iN,,,drop=FALSE], nrow = n.obs, ncol = n.jumps, byrow = TRUE),
                                                       G.iid = matrix(G.jump.iid[iN,1:n.jumps,,drop=FALSE], nrow = n.obs, ncol = n.jumps, byrow = TRUE),
                                                       dM.iid = matrix(dLambda.jump.iid[iN,,,drop=FALSE], nrow = n.obs, ncol = n.jumps, byrow = TRUE),
                                                       index.store = index.store, n.obs = n.obs, n.times = n.times, n.jumps = n.jumps)

                for(iC in 1:n.contrasts){ ## iC <- 1
                    iid.ate[[iC]][iN,] <- iid.ate[[iC]][iN,] + colMeans(colMultiply_cpp(term.augm, scale = iW.IPTW[,iC]))
                }
        }
    }

    ## ** export
    return(iid.ate)            
}

## * .calcAugmentationTerm_iid
## computes: int_{min(T,tau)}        (F'(tau) - F'(t)) / (G(t) * S(t)) * dM(t)
##                           - G(t)' (F(tau) - F(t)) / (G(t)^2 * S(t)) * dM(t)         
##                           - S(t)' (F(tau) - F(t)) / (G(t) * S(t)^2) * dM(t)         
##                           -       (F(tau) - F(t)) / (G(t) * S(t)) * dM'(t)         
.calcAugmentationTerm_iid <- function(F1, F1.tau, S, G, dM,
                                      F1.iid, F1.tau.iid, S.iid, G.iid, dM.iid,
                                      index.store, n.obs, n.times, n.jumps){

    K1 <- 1/(S*G)
    K2 <- K1 * dM

    int.tau.iid <- rowCumSum(K2)
    int.tau.t <- - rowCumSum(G.iid/G*K2 + S.iid/S*K2 + K1*dM.iid)
    int.t <- rowCumSum(- F1.iid*K2 + F1*G.iid/G*K2 + F1*S.iid/S*K2 + F1*K1*dM.iid)
    
    iTauStore <- 1:n.times

    out <- matrix(0, nrow = n.obs, ncol = n.times)
    for(iN in 1:n.obs){ ## iN <- 1
        subset <- which(index.store[iN,]>0)
        
        if(length(subset)>0){
            iJump <- index.store[iN,subset] 
            iTau <- iTauStore[subset]

            out[iN,iTau] <- (F1.tau.iid[iN,iTau] * int.tau.iid[iN,iJump] + F1.tau[iN,iTau] * int.tau.t[iN,iJump] + int.t[iN,iJump])
        }
    }

    return(out)
}

######################################################################
### ate-iid.R ends here
