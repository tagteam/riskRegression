## ate-pointEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2019 (10:43) 
## Version: 
## Last-Updated: Oct 17 2024 (11:44) 
##           By: Brice Ozenne
##     Update #: 1041
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * ATE_TI: compute average risk for time independent covariates
ATE_TI <- function(object.event,
                   object.treatment,
                   object.censor,
                   mydata,
                   treatment,
                   strata,
                   contrasts,
                   allContrasts,
                   times,
                   landmark,
                   cause,
                   level.censoring,
                   levels,
                   n.censor,
                   estimator,
                   eventVar.time,
                   eventVar.status,
                   censorVar.time,
                   censorVar.status,
                   return.iid.nuisance,
                   data.index,
                   method.iid,
                   product.limit,
                   store,
                   ...){

    tol <- 1e-12 ## difference in jump time must be above tol
    n.obs <- NROW(mydata)
    n.contrasts <- length(contrasts)
    n.times <- length(times)

    ## ** prepare output
    out <- list(meanRisk = NULL,
                diffRisk = NULL,
                ratioRisk = NULL,
                store = NULL)
    if(attr(estimator,"export.GFORMULA")){
        out$meanRisk <- rbind(out$meanRisk,
                              as.data.table(cbind(estimator = "GFORMULA", expand.grid(time = times, treatment = contrasts), estimate = as.numeric(NA)))
                              )

        out$store$iid.GFORMULA <- lapply(1:n.contrasts, function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
        names(out$store$iid.GFORMULA) <- contrasts
    }
    
    if(attr(estimator,"export.IPTW")){
        out$meanRisk <- rbind(out$meanRisk,
                           as.data.table(cbind(estimator = "IPTW", expand.grid(time = times, treatment = contrasts), estimate = as.numeric(NA)))
                           )

        out$store$iid.IPTW <- lapply(1:n.contrasts, function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
        names(out$store$iid.IPTW) <- contrasts
    }
    
    if(attr(estimator,"export.AIPTW")){
        out$meanRisk <- rbind(out$meanRisk,
                              as.data.table(cbind(estimator = "AIPTW", expand.grid(time = times, treatment = contrasts), estimate = as.numeric(NA)))
                              )
        
        out$store$iid.AIPTW <- lapply(1:n.contrasts, function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
        names(out$store$iid.AIPTW) <- contrasts
    }

    ## ** compute event indicators
    if(attr(estimator,"IPTW")){
        ## *** indicator for the outcome of interest stopped at time tau
        if(inherits(object.event,"glm") || (is.null(object.event) && is.null(object.censor))){
            time.before.tau <- cbind(mydata[[eventVar.status]])
        }else{
            time.before.tau <- sapply(times, function(tau){mydata[[eventVar.time]] <= tau})
        }
        Y.tau <- colMultiply_cpp(time.before.tau,
                                 scale = (mydata[[eventVar.status]] == cause)
                                 )

        ## *** treatment indicator
        M.treatment <- do.call(cbind,lapply(contrasts, "==", mydata[[treatment]]))
    }

    if(attr(estimator,"IPCW")){
        ## *** indicator for no censoring stopped at time tau
        ## C.tau <- colMultiply_cpp(time.before.tau,
        ##                          scale = (mydata[[censorVar.status]] != level.censoring)
        ##                          )
        C.tau <- 1-colMultiply_cpp(time.before.tau,
                                   scale = (mydata[[censorVar.status]] == level.censoring)
                                   )

        ## *** jump time for the censoring process
        time.jumpC <- unique(sort(mydata[[censorVar.time]][(mydata[[censorVar.status]] == level.censoring)]))
        index.obsSINDEXjumpC <- do.call(cbind,lapply(times, function(tau){
            prodlim::sindex(jump.times = time.jumpC, eval.times = pmin(mydata[[censorVar.time]],tau)-tol)
        }))
        index.lastjumpC <- max(index.obsSINDEXjumpC)
        time.jumpC <- time.jumpC[1:index.lastjumpC]
    }

    if(attr(estimator,"integral")){
        ## *** jump time index of the event time stopped at tau
        index.obsSINDEXjumpC.int <- do.call(cbind,lapply(times, function(tau){
            prodlim::sindex(jump.times = time.jumpC, eval.times = pmin(mydata[[eventVar.time]],tau))
        }))
        
        ## *** jump time of the censoring mecanism before event time
        beforeEvent.jumpC <- do.call(cbind,lapply(time.jumpC, function(iJump){iJump <= mydata[[eventVar.time]]}))
        beforeTau.nJumpC <- sapply(times, function(iTau){sum(time.jumpC <= iTau)})
        beforeTau.nJumpC.n0 <- beforeTau.nJumpC[beforeTau.nJumpC!=0]
    }

    ## ** compute predictions
    ## *** treatment model
    if(attr(estimator,"IPTW")){
        iPred <- lapply(contrasts, function(iC){predictRisk(object = object.treatment, newdata = mydata, levels = iC, iid = (method.iid==2)*return.iid.nuisance)})
        pi <- do.call(cbind,iPred)
        if(return.iid.nuisance && (method.iid==2)){
            out$store$iid.nuisance.treatment <- lapply(iPred,attr,"iid")
        }
    
        ## weights relative to the treatment
        iW.IPTW <- M.treatment / pi
    }

    ## *** censoring model
    if(attr(estimator,"IPCW")){
        iPred <- lapply(1:n.times, function(iT){ ## iT <- 1
            1-predictRisk(object.censor, newdata = mydata, times = c(0,time.jumpC)[index.obsSINDEXjumpC[,iT]+1],
                          diag = TRUE, product.limit = product.limit, iid = (method.iid==2)*return.iid.nuisance, store = store)
        })
            
        G.T_tau <- do.call(cbind,iPred)

        if(return.iid.nuisance && (method.iid==2)){
            out$store$iid.nuisance.censoring.diag <- lapply(iPred, function(x){-attr(x,"iid")})
        }

        ## weights relative to the censoring
        ## - because predictRisk outputs the risk instead of the survival
        iW.IPCW <- C.tau / G.T_tau
    }

    ## *** outcome model (computation of Prob[T<=t,Delta=1|A,W] = F_1(t|A=a,W))
    if(attr(estimator,"GFORMULA")){
        n.obs.contrasts <- rep(n.obs, n.contrasts)
        ls.index.strata <- vector(mode = "list", length = n.obs)
        F1.ctf.tau <- lapply(1:n.contrasts, function(x){
            matrix(0, nrow = n.obs, ncol = n.times,
                   dimnames = list(NULL, times))
        })
        names(F1.ctf.tau) <- contrasts
        
        for(iC in 1:n.contrasts){
            
            if(!is.null(treatment)){
                ## hypothetical world: in which every subject is treated with the same treatment
                ls.index.strata[[iC]] <- 1:n.obs
                data.i <- data.table::copy(mydata)
                data.i[[treatment]] <- factor(contrasts[iC], levels = levels)
            }else{
                ## hypothetical world: only patients with the same strata variable exist
                ls.index.strata[[iC]] <- which(mydata[[strata]]==contrasts[iC])
                data.i <- mydata[ls.index.strata[[iC]]]
                n.obs.contrasts[iC] <- length(ls.index.strata[[iC]])
            }
            
            if(return.iid.nuisance){
                factor <- TRUE
                attr(factor,"factor") <- list()

                if(attr(estimator,"export.GFORMULA")){
                    attr(factor,"factor") <- c(attr(factor,"factor"),
                                               list(GFORMULA = matrix(1, nrow =  n.obs.contrasts[iC], ncol = 1))
                                               )
                }
                if(attr(estimator,"export.AIPTW")){
                    attr(factor,"factor") <- c(attr(factor,"factor"),
                                               list(AIPTW = cbind(1-iW.IPTW[ls.index.strata[[iC]],iC]))
                                               )
                }
            }else{
                factor <- FALSE
            }
            outRisk <- predictRisk(object.event, newdata = data.i, times = times,
                                   average.iid = factor, cause = cause,
                                   product.limit = product.limit, store = store)

            F1.ctf.tau[[iC]][ls.index.strata[[iC]],] <- outRisk
            if(return.iid.nuisance){
                if(attr(estimator,"export.GFORMULA")){
                    out$store$iid.GFORMULA[[iC]] <- attr(outRisk,"average.iid")[["GFORMULA"]]
                }
                if(attr(estimator,"export.AIPTW")){
                    out$store$iid.AIPTW[[iC]] <- attr(outRisk,"average.iid")[["AIPTW"]]
                }
            }
        }
    }

    ## ** Compute augmentation term 
    if(attr(estimator,"integral")){
        ## integral is sum over individuals from 0 to min(T_i,\tau) of dM_i^C. It is therefore non-0 only if:
        ## - some of the requested times are after the first censoring time 
        ## - some of the censoring jump times are before the event times
        ## this would correspond to index.lastjumpC>0 && any(beforeEvent.jumpC)
        ## the first should be automatically enforced in ate_initArgs when creating attr(estimator,"integral")
        
        augTerm <- matrix(0, nrow = n.obs, ncol = n.times)

        ## absolute risk at event times and jump times of the censoring process
        predTempo <- predictRisk(object.event, newdata = mydata, times = c(times, time.jumpC), cause = cause, product.limit = product.limit,
                                 iid = (method.iid==2)*return.iid.nuisance, store = store)
        F1.tau <- predTempo[,1:n.times,drop=FALSE]
        F1.jump <- predTempo[,n.times + (1:index.lastjumpC),drop=FALSE]
        if((method.iid==2)*return.iid.nuisance){
            out$store$iid.nuisance.outcome <- attr(predTempo,"iid")
        }
        
        ## survival at all jump of the censoring process
        S.jump <- predictRisk(object.event, type = "survival", newdata = mydata, times = time.jumpC-tol, product.limit = product.limit,
                              iid = (method.iid==2)*return.iid.nuisance, store = store)
        if((method.iid==2)*return.iid.nuisance){
            out$store$iid.nuisance.survival <- attr(S.jump,"iid")
            attr(S.jump,"iid") <- NULL
        }

        ## martingale for the censoring process
        ## at all times of jump of the censoring process
        G.jump <- 1-predictRisk(object.censor, newdata = mydata, times = if(index.lastjumpC>1){c(0,time.jumpC[1:(index.lastjumpC-1)])}else{0},
                                product.limit = product.limit, iid = (method.iid==2)*return.iid.nuisance, store = store)
        
        if(return.iid.nuisance && (method.iid==2)){
            out$store$iid.nuisance.censoring <- -attr(G.jump,"iid")
            attr(G.jump,"iid") <- NULL
        }
        dN.jump <- do.call(cbind,lapply(time.jumpC, function(iJump){(mydata[[eventVar.time]] == iJump)*(mydata[[eventVar.status]] == level.censoring)}))
        dLambda.jump <- predictCox(object.censor, newdata = mydata, times = time.jumpC, type = "hazard", iid = (method.iid==2)*return.iid.nuisance, store = store)
        if((method.iid==2)*return.iid.nuisance){
            out$store$iid.nuisance.martingale <- dLambda.jump$hazard.iid
        }
        dM.jump <- dN.jump - dLambda.jump$hazard
        ## integral = \int_0^min(T_i,\tau) (F1(\tau|A_i,W_i)-F1(t|A_i,W_i)) / S(t|A_i,W_i) * dM_i^C(t)/Gc(t|A_i,W_i)
        ##          = F1(\tau|A_i,W_i) \int_0^min(T_i,\tauf) dM_i^C(t) / (S(t|A_i,W_i) * Gc(t|A_i,W_i)) - \int_0^min(T_i,\tauf) F1(t|A_i,W_i) dM_i^C(t) / (S(t|A_i,W_i) * Gc(t|A_i,W_i))
        integrand <- matrix(0, nrow = n.obs, ncol = index.lastjumpC)
        integrand2 <- matrix(0, nrow = n.obs, ncol = index.lastjumpC)
        index.beforeEvent.jumpC <- which(beforeEvent.jumpC) ## identify which increment to put in the integral, i.e. when the individual was still at risk of being censored
        integrand[index.beforeEvent.jumpC] <- dM.jump[index.beforeEvent.jumpC] / (G.jump[index.beforeEvent.jumpC] * S.jump[index.beforeEvent.jumpC])
        integrand2[index.beforeEvent.jumpC] <- F1.jump[index.beforeEvent.jumpC] * integrand[index.beforeEvent.jumpC]
        integral <- rowCumSum(integrand)
        integral2 <- rowCumSum(integrand2)
        augTerm[,beforeTau.nJumpC!=0] <- F1.tau[,beforeTau.nJumpC!=0,drop=FALSE] * integral[,beforeTau.nJumpC.n0,drop=FALSE] - integral2[,beforeTau.nJumpC.n0,drop=FALSE]
    }

    ## ** Compute individual contribution to the ATE + influence function for the G-formula
    for(iC in 1:n.contrasts){ ## iC <- 1
        if(attr(estimator,"export.GFORMULA")){
            if(!is.null(treatment)){
                iIID.ate <- F1.ctf.tau[[iC]]
                iATE <- colSums(iIID.ate)/n.obs
                out$meanRisk[list("GFORMULA",contrasts[iC]), c("estimate") := iATE, on = c("estimator","treatment")] 
                out$store$iid.GFORMULA[[iC]][data.index,] <- out$store$iid.GFORMULA[[iC]][data.index,] + rowCenter_cpp(iIID.ate, center = iATE)/n.obs
            }else{
                iIID.ate <- F1.ctf.tau[[iC]][ls.index.strata[[iC]],,drop=FALSE]
                iATE <- colSums(iIID.ate)/n.obs.contrasts[iC]                
                out$meanRisk[list("GFORMULA",contrasts[iC]), c("estimate") := iATE, on = c("estimator","treatment")]
                out$store$iid.GFORMULA[[iC]][data.index[ls.index.strata[[iC]]],] <- out$store$iid.GFORMULA[[iC]][data.index[ls.index.strata[[iC]]],] + rowCenter_cpp(iIID.ate, center = iATE)/n.obs.contrasts[iC]
            }
        }
        if(attr(estimator,"export.IPTW")){

            if(attr(estimator,"IPCW")){
                iIID.ate <- colMultiply_cpp(iW.IPCW * Y.tau, scale = iW.IPTW[,iC])
            }else{
                iIID.ate <- colMultiply_cpp(Y.tau, scale = iW.IPTW[,iC])
            }
            iATE <- colSums(iIID.ate)/n.obs
            if(attr(estimator,"monotone")){ ## ensure monotonicity over time (not accounted for in the standard error)
                out$meanRisk[list("IPTW",contrasts[iC]), c("estimate") := .monotonize(iATE), on = c("estimator","treatment")]
            }else{
                out$meanRisk[list("IPTW",contrasts[iC]), c("estimate") := iATE, on = c("estimator","treatment")]
            }
            out$store$iid.IPTW[[iC]][data.index,] <- out$store$iid.IPTW[[iC]][data.index,]  + rowCenter_cpp(iIID.ate, center = iATE)/n.obs
        }
        if(attr(estimator,"export.AIPTW")){
            if(attr(estimator,"integral")){ ## full augmentation
                iIID.ate <- F1.ctf.tau[[iC]] + colMultiply_cpp(iW.IPCW * Y.tau - F1.ctf.tau[[iC]] + augTerm, scale = iW.IPTW[,iC])
            }else if(inherits(object.event,"wglm") && attr(estimator,"IPCW")){ ## not full augmentation
                iIID.ate <- F1.ctf.tau[[iC]] + colMultiply_cpp(iW.IPCW * Y.tau - F1.ctf.tau[[iC]], scale = iW.IPTW[,iC])
            }else{
                iIID.ate <- F1.ctf.tau[[iC]] + colMultiply_cpp(Y.tau - F1.ctf.tau[[iC]], scale = iW.IPTW[,iC])
            }
            iATE <- colSums(iIID.ate)/n.obs
            if(attr(estimator,"monotone")){ ## ensure monotonicity over time (not accounted for in the standard error)
                out$meanRisk[list("AIPTW",contrasts[iC]), c("estimate") := .monotonize(iATE), on = c("estimator","treatment")]
            }else{
                out$meanRisk[list("AIPTW",contrasts[iC]), c("estimate") := iATE, on = c("estimator","treatment")]
            }
            out$store$iid.AIPTW[[iC]][data.index,] <- out$store$iid.AIPTW[[iC]][data.index,] + rowCenter_cpp(iIID.ate, center = iATE)/n.obs
        }
    }
        
    ## ** save quantities useful for the calculation of iid.nuisance
    if(return.iid.nuisance){
        out$store$n.obs <- n.obs
        out$store$n.times <- n.times
        
        if(attr(estimator,"GFORMULA")){
            out$store$F1.ctf.tau <- F1.ctf.tau            
        }
        
        if(attr(estimator,"IPTW")){
            out$store$iW.IPTW <- iW.IPTW
            out$store$iW.IPTW2 <- iW.IPTW / pi
            out$store$Y.tau <- Y.tau
        }

        if(attr(estimator,"IPCW")){
            out$store$iW.IPCW <- iW.IPCW
            out$store$iW.IPCW2 <- iW.IPCW / G.T_tau

            out$store$time.jumpC <- time.jumpC
            out$store$n.jumps <- index.lastjumpC
            out$store$index.obsSINDEXjumpC <- index.obsSINDEXjumpC
        }
        
        if(attr(estimator,"integral")){
            out$store$augTerm <- augTerm
            out$store$F1.tau <- F1.tau
            out$store$F1.jump <- F1.jump
            out$store$S.jump <- S.jump
            out$store$G.jump <- G.jump
            out$store$dM.jump <- dM.jump

            out$store$beforeEvent.jumpC <- beforeEvent.jumpC
            out$store$beforeTau.nJumpC <- beforeTau.nJumpC
            out$store$index.obsSINDEXjumpC.int <- index.obsSINDEXjumpC.int
        }
    }

    ## ** Compute risk comparisons
    out[c("diffRisk","ratioRisk")] <- ATE_COMPARISONS(out$meanRisk, TD = FALSE, allContrasts = allContrasts)
    
    ## ** Export
    return(out)            
}

## * ATE_TD: compute average risk for time dependent covariates (using G-formula)
ATE_TD <- function(object.event,
                   mydata,
                   formula,
                   treatment,
                   contrasts,
                   allContrasts,
                   times,
                   landmark,
                   cause,
                   n.contrasts,
                   levels,
                   ...){

    n.contrasts <- length(contrasts)

    response <- eval(formula[[2]],envir=mydata)
    time <- response[,"time"]
    entry <- response[,"entry"]
    if(inherits(x=object.event,what="coxph")){
        riskhandler <- "predictRisk.coxphTD"
    }else{
        riskhandler <- "predictRisk.CSCTD"
    }
    ## prediction for the hypothetical worlds in which every subject is treated with the same treatment
    out <- list(meanRisk = NULL,
                diffRisk = NULL,
                ratioRisk = NULL,
                store = NULL)
    out$meanRisk <- data.table::rbindlist(lapply(1:n.contrasts,function(i){
        data.i <- mydata
        data.i[[treatment]] <- factor(contrasts[i], levels = levels)
        data.table::rbindlist(lapply(landmark,function(lm){
            atrisk <- (entry <= lm & time >= lm)
            risk.i <- colMeans(do.call(riskhandler,
                                       args = list(object.event,
                                                   newdata = data.i[atrisk,],
                                                   times = times,
                                                   cause = cause,
                                                   landmark=lm,
                                                   ...)))
            data.table::data.table(estimator = "GFORMULA",
                                   treatment=contrasts[[i]],
                                   time=times,
                                   landmark=lm,
                                   estimate=risk.i)
        }))
    }))

    ## ** Compute risk comparisons
    out[c("diffRisk","ratioRisk")] <- ATE_COMPARISONS(out$meanRisk, TD = TRUE, allContrasts = allContrasts)

    ## ** Export
    return(out)
}

## * ATE_COMPARISONS: compute average risk for time dependent covariates (using G-formula)
ATE_COMPARISONS <- function(data, TD, allContrasts){
    ## duplicate
    dataA <- copy(data)
    setnames(dataA, old = c("treatment","estimate"), new = c("A","estimate.A"))
    dataB <- copy(data)
    setnames(dataB, old = c("treatment","estimate"), new = c("B","estimate.B"))

    ## perform all pairwise combinations
    mdata <- do.call(rbind,apply(allContrasts, 2, function(iC){ merge(dataA[iC[1],.SD, on = "A"],dataB[iC[2],.SD, on = "B"],
                                                                      by = c("estimator","time",if(TD){"landmark"}))
    })) ## iC <- c("T0","T1")

    ## re-order by estimator
    mdata <- mdata[order(factor(mdata$estimator, levels = unique(data$estimator)))]
    
    ## compute stats
    out <- list(diffRisk = data.table::copy(mdata),
                ratioRisk = data.table::copy(mdata))
    out$diffRisk[, c("estimate") := .SD$estimate.B - .SD$estimate.A]
    out$ratioRisk[, c("estimate") := .SD$estimate.B / .SD$estimate.A]

    ## export
    return(out)
}
## * .monotonize: enforce monotone constrain over time
## .monotonize(c(0.3708385, 0.4035446, 0.4346159, 0.4591131, 0.4844758, 0.4894976 ))
## .monotonize(c(0.3708385, 0.4035446, 0.4346159, 0.4591131, 0.4844758, 0.47 ))
## .monotonize(c(0.3708385, 0.37, 0.4346159, 0.4591131, 0.4844758, 0.4894976 ))
## .monotonize(c(0.3708385, 0.4035446, 0.40, 0.4591131, 0.4844758, 0.4894976 ))
## .monotonize(c(0.4, 0.35, 0.38, 0.37, 0.4844758, 0.4894976 ))
.monotonize <- function(x){
    if(length(x)==1){return(x)}

    p <- length(x)
    A <- matrix(0,p,p)
    A[lower.tri(A, diag = TRUE)] <- 1

    fit <- nnls::nnls(A = A, b = x)$fitted
    return(fit)
}
######################################################################
### ate-pointEstimate.R ends here
