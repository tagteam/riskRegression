## ate-pointEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2019 (10:43) 
## Version: 
## Last-Updated: okt  4 2019 (15:59) 
##           By: Brice Ozenne
##     Update #: 533
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
                   data,
                   formula,
                   treatment,
                   contrasts,
                   times,
                   landmark,
                   cause,
                   n.contrasts,
                   levels,
                   ...){

    Treatment <- Treatment.B <- meanRisk <- ratio <- NULL ## [:forCRANcheck:]
    n.contrasts <- length(contrasts)

    response <- eval(formula[[2]],envir=data)
    time <- response[,"time"]
    entry <- response[,"entry"]
    if(class(object.event)[[1]]=="coxph"){
        riskhandler <- "predictRisk.coxphTD"
    }else{
        riskhandler <- "predictRisk.CSCTD"
    }
    ## prediction for the hypothetical worlds in which every subject is treated with the same treatment
    dt.meanRisk <- data.table::rbindlist(lapply(1:n.contrasts,function(i){
        data.i <- data
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
            data.table::data.table(Treatment=contrasts[[i]],time=times,landmark=lm,meanRisk=risk.i)
        }))
    }))
    riskComparison <- data.table::rbindlist(lapply(1:(n.contrasts-1),function(i){
        data.table::rbindlist(lapply(((i+1):n.contrasts),function(j){
            ## compute differences between all pairs of treatments
            RC <- dt.meanRisk[Treatment==contrasts[[i]]]
            setnames(RC,"Treatment","Treatment.A")
            RC[,Treatment.B:=contrasts[[j]]]
            RC[,diff:=dt.meanRisk[Treatment==contrasts[[j]],meanRisk]-meanRisk]
            RC[,ratio:=dt.meanRisk[Treatment==contrasts[[j]],meanRisk]/meanRisk]
            RC[,meanRisk:=NULL]
            RC[]
        }))}))
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
                   data,
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
                   return.iid,
                   return.iid.nuisance,
                   method.iid,
                   product.limit,
                   ...){

    tol <- 1e-12 ## difference in jump time must be above tol
    n.obs <- NROW(data)
    n.contrasts <- length(contrasts)
    n.times <- length(times)

    ## ** prepare output
    out <- list()
    if(return.iid){
        attr(out,"iid.ate") <- vector(mode = "list", length = n.contrasts)
        names(attr(out,"iid.ate")) <- contrasts
    }
    if(return.iid.nuisance){
        attr(out,"iid.outcome") <- lapply(1:n.contrasts, function(iC){matrix(0, nrow = n.obs, ncol = n.times)})
        names(attr(out,"iid.outcome")) <- contrasts
    }
    ## point estimate
    meanRisk <- matrix(NA, nrow = n.contrasts, ncol = n.times,
                       dimnames = list(contrasts, times))
    
    ## ** compute event indicators
    if(attr(estimator,"IPTW")){
        ## *** indicator for the outcome of interest stopped at time tau
        if(inherits(object.event,"glm")){
            time.before.tau <- cbind(data[[eventVar.status]])
        }else{
            time.before.tau <- sapply(times, function(tau){data[[eventVar.time]] <= tau})
        }
        
        Y.tau <- colMultiply_cpp(time.before.tau,
                                 scale = (data[[eventVar.status]] == cause)
                                 )

        ## *** treatment indicator
        M.treatment <- do.call(cbind,lapply(contrasts, "==", data[[treatment]]))
    }

    if(attr(estimator,"IPCW")){
        ## *** indicator for no censoring stopped at time tau
        C.tau <- colMultiply_cpp(time.before.tau,
                                 scale = (data[[eventVar.status]] != level.censoring)
                                 )

        ## *** jump time for the censoring process
        time.jumpC <- sort(data[[eventVar.time]][(data[[eventVar.status]] == level.censoring)])

        index.obsSINDEXjumpC <- do.call(cbind,lapply(times, function(tau){
            prodlim::sindex(jump.times = time.jumpC, eval.times = pmin(data[[eventVar.time]],tau))
        }))
        index.lastjumpC <- max(index.obsSINDEXjumpC)
        time.jumpC <- time.jumpC[1:index.lastjumpC]

    }
    if(attr(estimator,"integral")){
        ## *** jump time of the censoring mecanism before event time
        beforeEvent.jumpC <- do.call(cbind,lapply(time.jumpC, function(iJump){iJump <= data[[eventVar.time]]}))
        beforeTau.nJumpC <- sapply(times, function(iTau){sum(time.jumpC <= iTau)})
        beforeTau.nJumpC.n0 <- beforeTau.nJumpC[beforeTau.nJumpC!=0]
    }

    ## ** compute predictions
    ## *** treatment model
    if(attr(estimator,"IPTW")){
        iPred <- lapply(contrasts, function(iC){predictRisk(object = object.treatment, newdata = data, levels = iC, iid = (method.iid==2)*return.iid.nuisance)})
        pi <- do.call(cbind,iPred)
        if(return.iid.nuisance && (method.iid==2)){
            attr(out,"iid.nuisance.treatment") <- lapply(iPred,attr,"iid")
        }
    
        ## weights relative to the treatment
        iW.IPTW <- M.treatment / pi
    }

    ## *** censoring model
    if(attr(estimator,"IPCW")){

        ## at all times of jump of the censoring process
        if(product.limit){
            G.jump <- predictCoxPL(object.censor, newdata = data, times = time.jumpC, iid = (method.iid==2)*return.iid.nuisance)
        }else{
            G.jump <- predictCox(object.censor, newdata = data, times = time.jumpC, iid = (method.iid==2)*return.iid.nuisance)
        }
        if(return.iid.nuisance && (method.iid==2)){
            attr(out,"iid.nuisance.censoring") <- G.jump$survival.iid
        }
    
        
        ## select the jump corresponding to the event time of each observation
        G.T_tau <- apply(index.obsSINDEXjumpC, 2, function(iCol){ ## iCol <- 2
            iOut <- rep(1,n.obs)
            iN0 <- which(iCol!=0)
            iOut[iN0] <- G.jump$survival[iN0 + (iCol[iN0]-1) * n.obs]
            return(iOut)
        })

        ## weights relative to the censoring
        iW.IPCW <- C.tau / G.T_tau

        ## set G at t-
        G.jump$survival <- cbind(1,G.jump$survival[,1:(index.lastjumpC-1),drop=FALSE])
    }

    ## *** outcome model (computation of Prob[T<=t,Delta=1|A,W] = F_1(t|A=a,W))
    n.obs.contrasts <- rep(n.obs, n.contrasts)
    if(attr(estimator,"Gformula")){
        F1.ctf.tau <- lapply(1:n.contrasts, function(x){
            matrix(0, nrow = n.obs, ncol = n.times,
                   dimnames = list(NULL, times))
        })
        names(F1.ctf.tau) <- contrasts
        
        for(iC in 1:n.contrasts){
            
            if(!is.null(treatment)){
                ## hypothetical world: in which every subject is treated with the same treatment
                index.strata <- 1:n.obs
                data.i <- data.table::copy(data)
                data.i[[treatment]] <- factor(contrasts[iC], levels = levels)
            }else{
                ## hypothetical world: only patients with the same strata variable exist
                index.strata <- which(data[[strata]]==contrasts[iC])
                data.i <- data[index.strata]
                n.obs.contrasts[iC] <- length(index.strata)
            }
            
            if(return.iid.nuisance){
                average.iid <- TRUE
                if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                    attr(average.iid,"factor") <- list(cbind(1-iW.IPTW[index.strata,iC]))
                }else{
                    attr(average.iid,"factor") <- list(matrix(1, nrow =  NROW(data.i), ncol = 1))
                }
            }else{
                average.iid <- FALSE
            }            
            outRisk <- predictRisk(object.event, newdata = data.i, times = times,
                                   average.iid = average.iid, cause = cause,
                                   product.limit = product.limit)
            F1.ctf.tau[[iC]][index.strata,] <- outRisk
            if(return.iid.nuisance){
                attr(out,"iid.outcome")[[iC]] <- attr(outRisk,"average.iid")[[1]]
            }

        }
    }
    
    ## ** Compute augmentation term    
    if(attr(estimator,"integral")){
        ## absolute risk at event times
        predTempo <- predictRisk(object.event, newdata = data, times = c(times, time.jumpC), cause = cause, product.limit = product.limit,
                                 iid = (method.iid==2)*return.iid.nuisance)
        F1.tau <- predTempo[,1:n.times,drop=FALSE]
        F1.jump <- predTempo[,n.times + (1:index.lastjumpC),drop=FALSE]
        if((method.iid==2)*return.iid.nuisance){
            attr(out,"iid.nuisance.outcome") <- attr(predTempo,"iid")
        }
        
        ## survival
        if(inherits(object.event,"CauseSpecificCox")){ ## competing risk case
            S.jump <- predict(object.event, type = "survival", newdata = data, times = time.jumpC-tol, product.limit = product.limit,
                              iid = (method.iid==2)*return.iid.nuisance)
        }else if(product.limit){ ## survival case
            S.jump <- predictCoxPL(object.event, type = "survival", newdata = data, times = time.jumpC-tol,
                                   iid = (method.iid==2)*return.iid.nuisance)
        }else{
            S.jump <- predictCox(object.event, type = "survival", newdata = data, times = time.jumpC-tol,
                                 iid = (method.iid==2)*return.iid.nuisance)
        }
        if((method.iid==2)*return.iid.nuisance){
            attr(out,"iid.nuisance.survival") <- S.jump$survival.iid
        }

        ## martingale for the censoring process
        dN.jump <- do.call(rbind,lapply(1:n.obs, function(iObs){(data[[eventVar.time]][iObs] == time.jumpC)*(data[[eventVar.status]][iObs] == level.censoring)}))
        dLambda.jump <- predictCox(object.censor, newdata = data, times = time.jumpC, type = "hazard", iid = (method.iid==2)*return.iid.nuisance)
        if((method.iid==2)*return.iid.nuisance){
            attr(out,"iid.nuisance.martingale") <- dLambda.jump$hazard.iid
        }

        dM.jump <- dN.jump - dLambda.jump$hazard

        ## integral
        integrand <- dM.jump * beforeEvent.jumpC / (G.jump$survival * S.jump$survival)
        integrand2 <- F1.jump * integrand
        integral <- rowCumSum(integrand)
        integral2 <- rowCumSum(integrand2)

        augTerm <- matrix(0, nrow = n.obs, ncol = n.times)
        augTerm[,beforeTau.nJumpC!=0] <- F1.tau[,beforeTau.nJumpC!=0,drop=FALSE] * integral[,beforeTau.nJumpC.n0,drop=FALSE] - integral2[,beforeTau.nJumpC.n0,drop=FALSE]
    }
       
    ## ** Compute individual contribution to the ATE
    for(iC in 1:n.contrasts){ ## iC <- 1
        ## compute influence function
        if(estimator == "Gformula"){
            iid.ate <- F1.ctf.tau[[iC]]
        }else if(estimator == "IPTW"){
            iid.ate <- colMultiply_cpp(Y.tau, scale = iW.IPTW[,iC])
        }else if(estimator == "AIPTW"){
            iid.ate <- F1.ctf.tau[[iC]] + colMultiply_cpp(Y.tau - F1.ctf.tau[[iC]], scale = iW.IPTW[,iC])
        }else if(estimator == "IPTW,IPCW"){
            iid.ate <- colMultiply_cpp(iW.IPCW * Y.tau, scale = iW.IPTW[,iC])
        }else if(estimator == "AIPTW,AIPCW"){
            iid.ate <- F1.ctf.tau[[iC]] + colMultiply_cpp(iW.IPCW * Y.tau - F1.ctf.tau[[iC]] + augTerm, scale = iW.IPTW[,iC])
        }
        
        ## estimate ate
        meanRisk[iC,] <- colSums(iid.ate)/n.obs.contrasts[iC]

        ## first term of the iid decomposition
        if(return.iid){
            ## center and scale iid decomposition for the functional delta method
            attr(out,"iid.ate")[[iC]] <- rowCenter_cpp(iid.ate * (n.obs/n.obs.contrasts[iC]), center = meanRisk[iC,])/n.obs
            dimnames(attr(out,"iid.ate")[[iC]]) <- list(NULL, times)
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
            attr(out,"S.jump") <- S.jump$survival
            attr(out,"G.jump") <- G.jump$survival
            attr(out,"dM.jump") <- dM.jump

            attr(out,"beforeEvent.jumpC") <- beforeEvent.jumpC
            attr(out,"beforeTau.nJumpC") <- beforeTau.nJumpC
        }
    }

    ## ** reshape results before exporting
    out$meanRisk <- melt(data.table(Treatment = rownames(meanRisk),meanRisk),
                         id.vars = "Treatment",
                         value.name = "meanRisk",
                         variable.name = "timeChar")
    out$meanRisk[,c("time") := times,by="Treatment"]
    out$meanRisk[,c("timeChar") := NULL]
    ## range(as.numeric(as.character(out$meanRisk$timeChar))-out$meanRisk$time)
    
    out$riskComparison <- data.table::rbindlist(lapply(1:(n.contrasts-1),function(i){ ## i <- 1
        data.table::rbindlist(lapply(((i+1):n.contrasts),function(j){ ## j <- 2
            ## compute differences between all pairs of treatments
            data.table(Treatment.A=contrasts[i],
                       Treatment.B=contrasts[j],
                       time=times,
                       diff=meanRisk[j,]-meanRisk[i,],
                       ratio=meanRisk[j,]/meanRisk[i,])
        }))}))

    return(out)            
}
                                        # }}}


######################################################################
### ate-pointEstimate.R ends here
