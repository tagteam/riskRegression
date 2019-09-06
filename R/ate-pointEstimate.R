### ate-pointEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2019 (10:43) 
## Version: 
## Last-Updated: sep  6 2019 (17:20) 
##           By: Brice Ozenne
##     Update #: 380
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
                   product.limit,
                   ...){

    tol <- 1e-12 ## difference in jump time must be above tol
    n.obs <- NROW(data)
    n.contrasts <- length(contrasts)
    n.times <- length(times)
    
    ## ** prepare output
    out <- list()
    ## point estimate
    meanRisk <- matrix(NA, nrow = n.contrasts, ncol = n.times,
                       dimnames = list(contrasts, times))

    ## ** compute event indicators
    if(attr(estimator,"IPTW")){
        ## *** indicator for the outcome of interest stopped at time tau
        time.before.tau <- sapply(times, function(tau){data[[eventVar.time]] <= tau})
        
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
        ## check
        ## X.censor <- coxModelFrame(object.censor)
        ## new.time <- X.censor$stop
        ## new.status <- X.censor$status
        ## jump.time <- sort(X.censor$stop[X.censor$status == 1]) ##  only select jumps

        index.maxjumpC <- do.call(cbind,lapply(times, function(tau){
            prodlim::sindex(jump.times = time.jumpC, eval.times = pmin(data[[eventVar.time]],tau))
        }))
        index.lastjumpC <- max(index.maxjumpC)
        time.jumpC <- time.jumpC[1:index.lastjumpC]

        timesInJump <- prodlim::sindex(jump.times = time.jumpC, eval.times = times)
    }

    ## ** compute predictions
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
                data.i <- data.table::copy(data)
                data.i[[treatment]] <- factor(contrasts[iC], levels = levels)
                F1.ctf.tau[[iC]] <- predictRisk(object.event, newdata = data.i, times = times, cause = cause, product.limit = product.limit)
            }else{
                ## hypothetical world: only patients with the same strata variable exist
                index.strata <- which(data[[strata]]==contrasts[iC])
                data.i <- data[index.strata]
                n.obs.contrasts[iC] <- length(index.strata)
                F1.ctf.tau[[iC]] <- matrix(0, nrow = n.obs, ncol = n.times)
                F1.ctf.tau[[iC]][index.strata,] <- predictRisk(object.event, newdata = data.i, times = times, cause = cause, product.limit = product.limit)
            }

        }

        if(attr(estimator,"integral")){
            predTempo <- predictRisk(object.event, newdata = data, times = c(times, time.jumpC), cause = cause, product.limit = product.limit)
            F1.tau <- predTempo[,1:n.times,drop=FALSE]
            F1.jump <- predTempo[,n.times + (1:index.lastjumpC),drop=FALSE]

            if(inherits(object.event,"CauseSpecificCox")){ ## competing risk case
                S.jump <- predict(object.event, type = "survival", newdata = data, times = time.jumpC, product.limit = product.limit, iid = FALSE)$survival
            }else{ ## survival case
                S.jump <- 1 - F1.jump
            }
        }

    }
    ## *** treatment model
    if(attr(estimator,"IPTW")){
        pi <- do.call(cbind,lapply(contrasts, function(iC){ ## iC <- 1
            predictRisk(object.treatment, newdata = data, level = iC)
        }))
        ## weights relative to the treatment
        iW.IPTW <- M.treatment / pi
    }

    ## *** censoring model
    if(attr(estimator,"IPCW")){

        ## at all times of jump of the censoring process
        G.jump <- 1 - predictRisk(object.censor, newdata = data, times = c(0,time.jumpC), product.limit = product.limit)
        
        ## at observation times
        G.T_tau <- matrix(NA, nrow = n.obs, ncol = n.times)
        for(iTau in 1:n.times){  # iTau <- 1
            ## select the jump corresponding to the event time of each observation
            index.col <- prodlim::sindex(jump.times = c(0,time.jumpC),
                                         eval.times = pmin(data[[eventVar.time]] - tol, times[iTau]))

            G.T_tau[,iTau] <- G.jump[(1:n.obs) + (index.col-1) * n.obs]
            ## G.T_tau[,iTau] - predictCox(object.censor, newdata = data, times = pmin(data[[eventVar.time]] - tol, times[iTau]), type = "survival", diag = TRUE)$survival[, 1]
            ## G.T_tau[,iTau] - predictCoxPL(object.censor, newdata = data, times = pmin(data[[eventVar.time]] - tol, times[iTau]), type = "survival", diag = TRUE)$survival[, 1]
        }
        
        ## weights relative to the censoring
        iW.IPCW <- C.tau / G.T_tau
        
#### for the augmentation term relative to the censoring process
        if(attr(estimator,"integral")){

            ## martingale for the censoring process
            dN.jump <- matrix(0, nrow = n.obs, ncol = index.lastjumpC,
                              dimnames = list(NULL, time.jumpC))
            row.index <- which( (data[[eventVar.status]] == level.censoring)*(data[[eventVar.time]] <= max(time.jumpC)) == 1 )
            for(iN in row.index){ ## iN <- row.index[1]
                dN.jump[iN,which((time.jumpC == data[[eventVar.time]][iN]))] <- 1
            }

            dLambda.jump <- predictCox(object.censor, newdata = data, times = time.jumpC, type = "hazard", ...)$hazard

            dM.jump <- dN.jump - dLambda.jump

        }
    }

    ## ** Compute augmentation term    
    if(attr(estimator,"integral")){
        augTerm <- .calcAugmentationTerm(G = G.jump[,1:index.lastjumpC], ## at t-
                                         F1 = F1.jump,
                                         F1.tau = F1.tau,
                                         S = S.jump, 
                                         dM = dM.jump,
                                         index.store = index.maxjumpC,
                                         n.obs = n.obs,
                                         n.times = n.times,
                                         n.jumps = max(index.lastjumpC))
    }
       
    ## cumsum( (F1.tau[1,1] - F1.jump[1,])/(S.jump[1,] * G.jump[1,1:index.lastjumpC]) * dM.jump[1,] )
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
        ## colMeans(colMultiply_cpp(iW.IPCW * Y.tau, scale = iW.IPTW[,iC]))
        
        ## estimate ate
        meanRisk[iC,] <- colSums(iid.ate)/n.obs.contrasts[iC]

        ## first term of the iid decomposition
        if(return.iid){
            ## center and scale iid decomposition for the functional delta method
            attr(out,"iid.ate")[[iC]] <- rowCenter_cpp(iid.ate * (n.obs/n.obs.contrasts), center = meanRisk[iC,])/n.obs
            dimnames(attr(out,"iid.ate")[[iC]]) <- list(NULL, times)
        }        
    }

    ## ** save quantities useful for the calculation of iid.nuisance
    if(return.iid){
        names(attr(out,"iid.ate")) <- contrasts
        
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
            attr(out,"n.obs") <- n.obs
            attr(out,"n.times") <- n.times
        }
        
        if(attr(estimator,"integral")){
            attr(out,"augTerm") <- augTerm
            attr(out,"F1.tau") <- F1.tau
            attr(out,"F1.jump") <- F1.jump
            attr(out,"S.jump") <- S.jump
            attr(out,"G.jump") <- G.jump
            attr(out,"dM.jump") <- dM.jump

            attr(out,"index.store") <- index.maxjumpC
            attr(out,"n.jumps") <- max(index.lastjumpC)
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

## * .calcAugmentationTerm
## computes: int_{min(T,tau)} (F(tau) - F(t)) / (G(t) * S(t)) * dM(t) = F(tau) sum_{t \in min(T,tau)} K1(t) - sum_{t \in min(T,tau)} K2(t) 
.calcAugmentationTerm <- function(F1, F1.tau, S, G, dM, 
                                  index.store, n.obs, n.times, n.jumps){

    K1 <- rowCumSum(dM / (G*S))
    K2 <- rowCumSum(F1 * dM / (G*S))
    
    iTauStore <- 1:n.times

    out <- matrix(0, nrow = n.obs, ncol = n.times)
    for(iN in 1:n.obs){ ## iN <- 1
        subset <- which(index.store[iN,]>0)
        
        if(length(subset)>0){
            iJump <- index.store[iN,subset] 
            iTau <- iTauStore[subset]

            out[iN,iTau] <- (F1.tau[iN,iTau] * K1[iN,iJump] - K2[iN,iJump])
        }
    }

    return(out)
}

######################################################################
### ate-pointEstimate.R ends here
