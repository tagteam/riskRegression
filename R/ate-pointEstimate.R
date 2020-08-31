## ate-pointEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2019 (10:43) 
## Version: 
## Last-Updated: aug 31 2020 (15:25) 
##           By: Brice Ozenne
##     Update #: 751
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

                                        # {{{ Gformula: time dependent covariates
## * ATE_TD
ATE_TD <- function(object.event,
                   mydata,
                   formula,
                   treatment,
                   contrasts,
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
    if(class(object.event)[[1]]=="coxph"){
        riskhandler <- "predictRisk.coxphTD"
    }else{
        riskhandler <- "predictRisk.CSCTD"
    }
    ## prediction for the hypothetical worlds in which every subject is treated with the same treatment
    dt.meanRisk <- data.table::rbindlist(lapply(1:n.contrasts,function(i){
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
            data.table::data.table(treatment=contrasts[[i]],
                                   time=times,
                                   landmark=lm,
                                   meanRisk.Gformula=risk.i)
        }))
    }))
    riskComparison <- data.table::rbindlist(lapply(1:(n.contrasts-1),function(i){
        data.table::rbindlist(lapply(((i+1):n.contrasts),function(j){
            ## compute differences between all pairs of treatments
            iDT <- dt.meanRisk[dt.meanRisk$treatment==contrasts[[i]]]
            setnames(iDT,"treatment","treatment.A")
            baseRisk <- iDT$meanRisk.Gformula
            iDT[,c("meanRisk.Gformula"):=NULL]

            newRisk <- dt.meanRisk[treatment==contrasts[[j]], .SD$meanRisk.Gformula]

            iDT[,c("treatment.B","diff.Gformula","ratio.Gformula") := list(contrasts[[j]],newRisk-baseRisk,newRisk/baseRisk)]
            return(iDT[])
        }))}))
    setcolorder(riskComparison, neworder = c("treatment.A","treatment.B", setdiff(names(riskComparison),c("treatment.A","treatment.B"))))
    out <- list(meanRisk = dt.meanRisk,
                riskComparison = riskComparison,
                treatment = treatment,
                strata = strata)
    return(out)
}

# }}}

                                        # {{{ Gformula: time independent covariates
## * ATE_TI
ATE_TI <- function(object.event,
                   object.treatment,
                   object.censor,
                   mydata,
                   treatment,
                   strata,
                   contrasts,
                   times,
                   landmark,
                   cause,
                   level.censoring,
                   n.contrasts,
                   levels,
                   n.censor,
                   estimator,
                   eventVar.time,
                   eventVar.status,
                   censorVar.time,
                   censorVar.status,
                   type.multistate,
                   return.iid.nuisance,
                   data.index,
                   method.iid,
                   product.limit,
                   ...){

    tol <- 1e-12 ## difference in jump time must be above tol
    n.obs <- NROW(mydata)
    n.contrasts <- length(contrasts)
    n.times <- length(times)

    ## ** prepare output
    out <- list()
    meanRisk <- list()
    
    if(attr(estimator,"export.Gformula")){
        meanRisk <- c(meanRisk,
                      Gformula = list(matrix(0, nrow = n.contrasts, ncol = n.times,
                                             dimnames = list(contrasts, times))))

        attr(out,"iid.Gformula") <- lapply(1:n.contrasts, function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
        names(attr(out,"iid.Gformula")) <- contrasts
    }
    
    if(attr(estimator,"export.IPTW")){
        meanRisk <- c(meanRisk,
                      IPTW = list(matrix(0, nrow = n.contrasts, ncol = n.times,
                                         dimnames = list(contrasts, times))))

        attr(out,"iid.IPTW") <- lapply(1:n.contrasts, function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
        names(attr(out,"iid.IPTW")) <- contrasts
    }
    
    if(attr(estimator,"export.AIPTW")){
        meanRisk <- c(meanRisk,
                      AIPTW = list(matrix(0, nrow = n.contrasts, ncol = n.times,
                                          dimnames = list(contrasts, times))))

        attr(out,"iid.AIPTW") <- lapply(1:n.contrasts, function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
        names(attr(out,"iid.AIPTW")) <- contrasts
    }

    ## ** compute event indicators
    if(attr(estimator,"IPTW")){
        ## *** indicator for the outcome of interest stopped at time tau
        if(inherits(object.event,"glm")){
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
        C.tau <- colMultiply_cpp(time.before.tau,
                                 scale = (mydata[[censorVar.status]] != level.censoring)
                                 )

        ## *** jump time for the censoring process
        time.jumpC <- sort(mydata[[censorVar.time]][(mydata[[censorVar.status]] == level.censoring)])
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
            attr(out,"iid.nuisance.treatment") <- lapply(iPred,attr,"iid")
        }
    
        ## weights relative to the treatment
        iW.IPTW <- M.treatment / pi
    }

    ## *** censoring model
    if(attr(estimator,"IPCW")){
        iPred <- lapply(1:n.times, function(iT){
            1-predictRisk(object.censor, newdata = mydata, times = c(0,time.jumpC)[index.obsSINDEXjumpC[,iT]+1],
                        diag = TRUE, product.limit = product.limit, iid = (method.iid==2)*return.iid.nuisance)
        })
            
        G.T_tau <- do.call(cbind,iPred)

        if(return.iid.nuisance && (method.iid==2)){
            attr(out,"iid.nuisance.censoring.diag") <- lapply(iPred, function(x){-attr(x,"iid")})
        }

        ## weights relative to the censoring
        ## - because predictRisk outputs the risk instead of the survival
        iW.IPCW <- C.tau / G.T_tau
    }

    ## *** outcome model (computation of Prob[T<=t,Delta=1|A,W] = F_1(t|A=a,W))
    if(attr(estimator,"Gformula")){
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
                
                if(attr(estimator,"export.Gformula")){
                    attr(factor,"factor") <- c(attr(factor,"factor"),
                                               list(Gformula = matrix(1, nrow =  n.obs.contrasts[iC], ncol = 1))
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
                                   product.limit = product.limit)

            F1.ctf.tau[[iC]][ls.index.strata[[iC]],] <- outRisk
            if(return.iid.nuisance){
                if(attr(estimator,"export.Gformula")){
                    attr(out,"iid.Gformula")[[iC]] <- attr(outRisk,"average.iid")[["Gformula"]]
                }
                if(attr(estimator,"export.AIPTW")){
                    attr(out,"iid.AIPTW")[[iC]] <- attr(outRisk,"average.iid")[["AIPTW"]]
                }
            }
        }
    }
    
    ## ** Compute augmentation term
    if(attr(estimator,"integral")){
        if(inherits(object.event,"glm")){ ## WARNING: we need a proper estimator of F1
            browser()
            predTempo <- predictRisk(object.event, newdata = mydata, iid = (method.iid==2)*return.iid.nuisance)
            augTerm <- colMultiply_cpp( (1-iW.IPCW) * F1.ctf.tau[[iC]], scale = iW.IPTW[,iC])
        }else{
            ## absolute risk at event times
            predTempo <- predictRisk(object.event, newdata = mydata, times = c(times, time.jumpC), cause = cause, product.limit = product.limit,
                                     iid = (method.iid==2)*return.iid.nuisance)
            F1.tau <- predTempo[,1:n.times,drop=FALSE]
            F1.jump <- predTempo[,n.times + (1:index.lastjumpC),drop=FALSE]
            if((method.iid==2)*return.iid.nuisance){
                attr(out,"iid.nuisance.outcome") <- attr(predTempo,"iid")
            }
        
            ## survival at all jump of the censoring process
            S.jump <- predictRisk(object.event, type = "survival", newdata = mydata, times = time.jumpC-tol, product.limit = product.limit,
                                  iid = (method.iid==2)*return.iid.nuisance)
            if((method.iid==2)*return.iid.nuisance){
                attr(out,"iid.nuisance.survival") <- attr(S.jump,"iid")
                attr(S.jump,"iid") <- NULL
            }

            ## martingale for the censoring process
            ## at all times of jump of the censoring process
            G.jump <- 1-predictRisk(object.censor, newdata = mydata, times = if(index.lastjumpC>0){c(0,time.jumpC[1:(index.lastjumpC-1)])}else{0},
                                    product.limit = product.limit, iid = (method.iid==2)*return.iid.nuisance)
        
            if(return.iid.nuisance && (method.iid==2)){
                attr(out,"iid.nuisance.censoring") <- -attr(G.jump,"iid")
                attr(G.jump,"iid") <- NULL
            }
        
            dN.jump <- do.call(rbind,lapply(1:n.obs, function(iObs){(mydata[[eventVar.time]][iObs] == time.jumpC)*(mydata[[eventVar.status]][iObs] == level.censoring)}))
            dLambda.jump <- predictCox(object.censor, newdata = mydata, times = time.jumpC, type = "hazard", iid = (method.iid==2)*return.iid.nuisance)
            if((method.iid==2)*return.iid.nuisance){
                attr(out,"iid.nuisance.martingale") <- dLambda.jump$hazard.iid
            }

            dM.jump <- dN.jump - dLambda.jump$hazard

            ## integral
            integrand <- dM.jump * beforeEvent.jumpC / (G.jump * S.jump)
            integrand2 <- F1.jump * integrand
            integral <- rowCumSum(integrand)
            integral2 <- rowCumSum(integrand2)
            augTerm <- matrix(0, nrow = n.obs, ncol = n.times)
            augTerm[,beforeTau.nJumpC!=0] <- F1.tau[,beforeTau.nJumpC!=0,drop=FALSE] * integral[,beforeTau.nJumpC.n0,drop=FALSE] - integral2[,beforeTau.nJumpC.n0,drop=FALSE]
        }
    }
       
    ## ** Compute individual contribution to the ATE + influence function for the Gformula
    for(iC in 1:n.contrasts){ ## iC <- 1
        if(attr(estimator,"export.Gformula")){
            if(!is.null(treatment)){
                iIID.ate <- F1.ctf.tau[[iC]]
                meanRisk$Gformula[iC,] <- colSums(iIID.ate)/n.obs
                attr(out,"iid.Gformula")[[iC]][data.index,] <- attr(out,"iid.Gformula")[[iC]][data.index,] + rowCenter_cpp(iIID.ate, center = meanRisk$Gformula[iC,])/n.obs
            }else{
                iIID.ate <- F1.ctf.tau[[iC]][ls.index.strata[[iC]],,drop=FALSE]
                meanRisk$Gformula[iC,] <- colSums(iIID.ate)/n.obs.contrasts[iC]
                attr(out,"iid.Gformula")[[iC]][data.index[ls.index.strata[[iC]]],] <- attr(out,"iid.Gformula")[[iC]][data.index[ls.index.strata[[iC]]],] + rowCenter_cpp(iIID.ate, center = meanRisk$Gformula[iC,])/n.obs.contrasts[iC]
            }
        }
        if(attr(estimator,"export.IPTW")){
            if(attr(estimator,"IPCW")){
                iIID.ate <- colMultiply_cpp(iW.IPCW * Y.tau, scale = iW.IPTW[,iC])
            }else{
                iIID.ate <- colMultiply_cpp(Y.tau, scale = iW.IPTW[,iC])
            }
            
            meanRisk$IPTW[iC,] <- colSums(iIID.ate)/n.obs
            attr(out,"iid.IPTW")[[iC]][data.index,] <- attr(out,"iid.IPTW")[[iC]][data.index,]  + rowCenter_cpp(iIID.ate, center = meanRisk$IPTW[iC,])/n.obs
        }
        if(attr(estimator,"export.AIPTW")){
            if(attr(estimator,"IPCW")){
                iIID.ate <- F1.ctf.tau[[iC]] + colMultiply_cpp(iW.IPCW * Y.tau - F1.ctf.tau[[iC]] + augTerm, scale = iW.IPTW[,iC])
            }else{
                iIID.ate <- F1.ctf.tau[[iC]] + colMultiply_cpp(Y.tau - F1.ctf.tau[[iC]], scale = iW.IPTW[,iC])
            }

            meanRisk$AIPTW[iC,] <- colSums(iIID.ate)/n.obs
            attr(out,"iid.AIPTW")[[iC]][data.index,] <- attr(out,"iid.AIPTW")[[iC]][data.index,] + rowCenter_cpp(iIID.ate, center = meanRisk$AIPTW[iC,])/n.obs
        }
    }

    ## ** save quantities useful for the calculation of iid.nuisance
    if(return.iid.nuisance){
        attr(out,"n.obs") <- n.obs
        attr(out,"n.times") <- n.times
        
        if(attr(estimator,"Gformula")){
            attr(out,"F1.ctf.tau") <- F1.ctf.tau            
        }
        
        if(attr(estimator,"IPTW")){
            attr(out,"iW.IPTW") <- iW.IPTW
            attr(out,"iW.IPTW2") <- iW.IPTW / pi
            attr(out,"Y.tau") <- Y.tau
        }

        if(attr(estimator,"IPCW")){
            attr(out,"iW.IPCW") <- iW.IPCW
            attr(out,"iW.IPCW2") <- iW.IPCW / G.T_tau

            attr(out,"time.jumpC") <- time.jumpC
            attr(out,"n.jumps") <- index.lastjumpC
            attr(out,"index.obsSINDEXjumpC") <- index.obsSINDEXjumpC
        }
        
        if(attr(estimator,"integral")){
            attr(out,"augTerm") <- augTerm
            attr(out,"F1.tau") <- F1.tau
            attr(out,"F1.jump") <- F1.jump
            attr(out,"S.jump") <- S.jump
            attr(out,"G.jump") <- G.jump
            attr(out,"dLambda.jump") <- dLambda.jump$hazard
            attr(out,"dM.jump") <- dM.jump

            attr(out,"beforeEvent.jumpC") <- beforeEvent.jumpC
            attr(out,"beforeTau.nJumpC") <- beforeTau.nJumpC
            attr(out,"index.obsSINDEXjumpC.int") <- index.obsSINDEXjumpC.int
        }
    }

    ## ** reshape results before exporting
    meanRiskL <- lapply(names(meanRisk), function(iE){ ## iE <- "Gformula"
        iDT <- melt(data.table(treatment = rownames(meanRisk[[iE]]),meanRisk[[iE]]),
                    id.vars = "treatment",
                    value.name = paste0("meanRisk.",iE),
                    variable.name = "time")
        if(iE==names(meanRisk)[[1]]){
            return(iDT)
        }else{
            return(iDT[,.SD,.SDcols = paste0("meanRisk.",iE)])
        }
    })
    out$meanRisk <- do.call(cbind,meanRiskL)
    if(all(is.na(times))){
        out$meanRisk[, time := as.numeric(NA)]
    }else{
        out$meanRisk[, time := as.numeric(levels(time))[time]] ## recommanded way to convert from factor to numeric (https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information)
        ## range(as.numeric(as.character(out$meanRisk$timeChar))-out$meanRisk$time)
    }
    out$riskComparison <- data.table::rbindlist(lapply(1:(n.contrasts-1),function(i){ ## i <- 1
        data.table::rbindlist(lapply(((i+1):n.contrasts),function(j){ ## j <- 2

            ## compute differences between all pairs of treatments
            iDT <- data.table(treatment.A=contrasts[i],
                              treatment.B=contrasts[j],
                              time=times)
            if(attr(estimator,"export.Gformula")){
                iDT[, c("diff.Gformula","ratio.Gformula") := list(meanRisk$Gformula[j,]-meanRisk$Gformula[i,],
                                                                  meanRisk$Gformula[j,]/meanRisk$Gformula[i,])]
            }
            if(attr(estimator,"export.IPTW")){
                iDT[, c("diff.IPTW","ratio.IPTW") := list(meanRisk$IPTW[j,]-meanRisk$IPTW[i,],
                                                          meanRisk$IPTW[j,]/meanRisk$IPTW[i,])]
            }
            if(attr(estimator,"export.AIPTW")){
                iDT[, c("diff.AIPTW","ratio.AIPTW") := list(meanRisk$AIPTW[j,]-meanRisk$AIPTW[i,],
                                                            meanRisk$AIPTW[j,]/meanRisk$AIPTW[i,])]
            }
            return(iDT)
        }))}))

    return(out)            
}
                                        # }}}


######################################################################
### ate-pointEstimate.R ends here
