### ate-pointEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2019 (10:43) 
## Version: 
## Last-Updated: jun 28 2019 (15:53) 
##           By: Brice Ozenne
##     Update #: 86
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
                   type.multistate,
                   return.iid,
                   ...){

    n.obs <- NROW(data)
    n.contrasts <- length(contrasts)
    n.times <- length(times)
    out <- list()
    
    ## ** compute augmentation term (AIPCW)
    if(estimator == "AIPTW,AIPCW"){
        augTerm <- matrix(0, nrow = n.obs, ncol = n.times)
        prob.event <- predictRisk(object.event, newdata = data, times = times, cause = cause,...)
        for(iTau in 1:n.times){ ## iTau <- 2
            if(n.censor[iTau]>0){
                data$prob.event <- prob.event[,iTau]
                augTerm[,iTau] <- .calcLterm(data = data, n.obs = n.obs, times = times[iTau],
                                             model.censor = object.censor,
                                             model.event = object.event,
                                             type = type.multistate,
                                             predictor.cox = "predictCox",
                                             product.limit = switch(type.multistate,
                                                                    "survival" = FALSE,
                                                                    "competing.risks" = TRUE),
                                             cause = cause)
            }
        }
    }

    ## ** estimator of the average risk
    meanRisk <- matrix(NA, nrow = n.contrasts, ncol = n.times,
                       dimnames = list(contrasts, times))
    if(return.iid){
        attr(out,"iid") <- lapply(1:n.contrasts, function(x){matrix(NA, nrow = n.obs, ncol = n.times)})
        if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
            attr(out,"prob.event") <- lapply(1:n.contrasts, function(x){matrix(NA, nrow = n.obs, ncol = n.times)})
        }
        if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
            attr(out,"prob.treatment") <- lapply(1:n.contrasts, function(x){matrix(NA, nrow = n.obs, ncol = n.times)})
        }
        if(estimator %in% c("IPTW,IPCW","AIPTW,AIPCW")){
            attr(out,"prob.censor") <- lapply(1:n.contrasts, function(x){matrix(NA, nrow = n.obs, ncol = n.times)})
        }
        if(estimator == "AIPTW,AIPCW"){
            attr(out,"augTerm") <- augTerm
        }
    }
    
    for(iC in 1:n.contrasts){ ## iC <- 1
        iIF <- matrix(0, nrow = n.obs, ncol = n.times)
        
        ## ** hypothetical world
        if(!is.null(treatment)){
            ## in which every subject is treated with the same treatment
            data.i <- data
            data.i[[treatment]] <- factor(contrasts[iC], levels = levels)
        }else{
            ## only patients with the same strata variable exist
            data.i <- data[data[[strata]]==contrasts[iC]]
        }
        
        ## ** G-formula
        if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
            iProb.event <- predictRisk(object.event, newdata = data.i, times = times, cause = cause,...)
            if(!is.matrix(iProb.event)){iProb.event <- cbind(iProb.event)}
            if(return.iid){attr(out,"prob.event")[[iC]] <- iProb.event}
            iIF <- iIF + iProb.event
        }
        
        ## ** Inverse probability weighting        
        if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
            for(iTau in 1:n.times){ ## iTau <- 1
                ## random variables stopped at times
                status.tau <- (data.i[[eventVar.time]] <= times[iTau]) * (data.i[[eventVar.time]] == cause)

                ## ** IPTW,IPCW
                ## compute IPTW
                iProb.treatment <- predictRisk(object.treatment, newdata = data)
                if(return.iid){attr(out,"prob.treatment")[[iC]][,iTau] <- iProb.treatment}
                iW.IPTW <- (data[[treatment]] == contrasts[iC]) / iProb.treatment 
                
                ## compute IPCW
                if(n.censor[iTau]==0){
                    iW.IPCW <- rep(1,NROW(data))
                }else{
                    Ncensoring.tau <- (data.i[[eventVar.time]] <= times[iTau]) * (data.i[[eventVar.time]] != level.censoring)
                    time.tau <- pmin(data.i[[eventVar.time]], times[iTau])

                    iProb.censor <- predictCox(object.censor, newdata = data, times = time.tau-(1e-10), type = "survival", diag = TRUE)$survival[,1]
                    if(return.iid){attr(out,"prob.censor")[[iC]][,iTau] <- iProb.censor}
                    iW.IPCW <- Ncensoring.tau / iProb.censor
                }

                ## assemble
                iIF[,iTau] <- iIF[,iTau] + status.tau * iW.IPCW * iW.IPTW ## IPCW,IPTW

                ## ** AIPTW,AIPCW
                ## assemble
                if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                    iIF[,iTau] <- iIF[,iTau] + iProb.event[,iTau] * iW.IPTW ## AIPTW
                }
                if(estimator %in% c("AIPTW,AIPCW")){
                    iIF[,iTau] <- iIF[,iTau] + augTerm[,iTau] * iW.IPTW ## AIPCW
                }
            }
        }

        meanRisk[iC,] <- colMeans(iIF)
        if(return.iid){attr(out,"iid")[[iC]] <- iIF}
    }
    
    riskComparison <- data.table::rbindlist(lapply(1:(n.contrasts-1),function(i){ ## i <- 1
        data.table::rbindlist(lapply(((i+1):n.contrasts),function(j){ ## j <- 2
            ## compute differences between all pairs of treatments
            data.table(Treatment.A=contrasts[i],
                       Treatment.B=contrasts[j],
                       time = times,
                       diff=meanRisk[j,]-meanRisk[i,],
                       ratio=meanRisk[j,]/meanRisk[i,])
        }))}))

    ## reshape for export
    name.strata <- unlist(lapply(1:n.contrasts, function(c){rep(contrasts[c],length(meanRisk[[c]]))}))

    out$meanRisk <- data.table(Treatment=name.strata,
                               time = times,
                               meanRisk=as.double(meanRisk))
    out$riskComparison <- riskComparison
    return(out)            
}
                                        # }}}



######################################################################
### ate-pointEstimate.R ends here
