### ate-pointEstimate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2019 (10:43) 
## Version: 
## Last-Updated: jul  4 2019 (11:03) 
##           By: Brice Ozenne
##     Update #: 187
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
                   type.multistate,
                   return.iid,
                   ...){

    n.obs <- NROW(data)
    n.contrasts <- length(contrasts)
    n.times <- length(times)

    ## ** prepare output
    out <- list()
    ## point estimate
    meanRisk <- matrix(NA, nrow = n.contrasts, ncol = n.times,
                       dimnames = list(contrasts, times))
    if(return.iid){ ## iid decomposition + useful quantities
        attr(out,"iid") <- vector(mode = "list", length = n.contrasts)
        if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
            attr(out,"prob.event") <- lapply(1:n.contrasts, function(x){matrix(NA, nrow = n.obs, ncol = n.times)})
        }
        if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
            attr(out,"prob.treatment") <- matrix(NA, nrow = n.obs, ncol = n.contrasts)
        }
    }

    ## ** compute indicators
    if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
        ## indicator for the outcome of interest stopped at time tau
        time.before.tau <- sapply(times, function(tau){data[[eventVar.time]] <= tau})
        
        status.tau <- colMultiply_cpp(time.before.tau,
                                      scale = (data[[eventVar.status]] == cause)
                                      )
        if(return.iid){
            attr(out,"status.tau") <- status.tau
        }

        if(estimator %in% c("IPTW,IPCW","AIPTW,AIPCW")){
            ## indicator for no censoring stopped at time tau
            Ncensoring.tau <- colMultiply_cpp(time.before.tau,
                                              scale = (data[[eventVar.status]] != level.censoring)
                                              )
            if(return.iid){
                attr(out,"Ncensoring.tau") <- Ncensoring.tau
            }
        }
    }

    ## ** IPCW                           
    if(estimator %in% c("IPTW,IPCW","AIPTW,AIPCW")){
        iProb.censor <- sapply(times, function(tau){ ## tau <- 1
            predictCox(object.censor,
                       newdata = data,
                       times = pmin(tau,data[[eventVar.time]])-(1e-10),
                       type = "survival",
                       diag = TRUE)$survival[,1]
        })
        if(return.iid){
            attr(out,"prob.censor") <- iProb.censor
        }
        iW.IPCW <- Ncensoring.tau / iProb.censor
        
    }else if(estimator %in% c("IPTW","AIPTW")){
        iW.IPCW <- matrix(1, nrow = n.obs, ncol = n.times)
    }

    ## ** compute augmentation term    
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
        if(return.iid){
            attr(out,"augTerm") <- augTerm
        }
    }else if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW")){
        augTerm <- matrix(0, nrow = n.obs, ncol = n.times)
    }

    ## ** estimator of the average risk
    for(iC in 1:n.contrasts){ ## iC <- 1
        iIF <- matrix(0, nrow = n.obs, ncol = n.times)
        
        ## ** hypothetical world
        if(!is.null(treatment)){
            ## in which every subject is treated with the same treatment
            data.i <- data
            data.i[[treatment]] <- factor(contrasts[iC], levels = levels)
            n.strata <- n.obs
        }else{
            ## only patients with the same strata variable exist
            index.strata <- which(data[[strata]]==contrasts[iC])
            n.strata <- length(index.strata)
            data.i <- data[index.strata]
        }
        
        ## ** G-formula
        if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
            iProb.event <- predictRisk(object.event, newdata = data.i, times = times, cause = cause,...)
            if(!is.matrix(iProb.event)){iProb.event <- cbind(iProb.event)}

            if(!is.null(treatment)){
                iIF <- iIF + iProb.event
            }else{
                iIF[index.strata,] <- iIF[index.strata,] + iProb.event
            }
        }else if(estimator %in% c("IPTW","IPTW,IPCW")){
            iProb.event <- matrix(0, nrow = n.obs, ncol = n.times)
        }
        
        ## ** Inverse probability weighting        
        if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
            ## IPTW
            iProb.treatment <- predictRisk(object.treatment, newdata = data, level = contrasts[iC])
            iW.IPTW <- (data[[treatment]] == contrasts[iC]) / iProb.treatment 

            ## assemble (also with augmentation terms)
            iIF <- iIF + colMultiply_cpp(status.tau * iW.IPCW - iProb.event + augTerm,
                                         scale= iW.IPTW)

        }

        ## ** estimate ATE
        meanRisk[iC,] <- colSums(iIF)/n.strata

        ## **  first term of the iid decomposition
        if(return.iid){
            ## center and scale iid decomposxition for the functional delta method
            attr(out,"iid")[[iC]] <- rowCenter_cpp(iIF * (n.obs/n.strata), center = meanRisk[iC,])/n.obs
            names(attr(out,"iid")) <- contrasts

            ## save other useful quantities
            if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
                if(!is.null(treatment)){
                    attr(out,"prob.event")[[iC]] <- iProb.event
                }else{
                    attr(out,"prob.event")[[iC]][index.strata,] <- iProb.event
                }
            }
            if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
                attr(out,"prob.treatment")[,iC] <- iProb.treatment
            }
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
