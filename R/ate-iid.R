### ate-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: jul  2 2019 (15:58) 
##           By: Brice Ozenne
##     Update #: 385
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * iidATE
iidATE <- function(meanRisk,
                   riskComparison,
                   object.event,
                   object.treatment,
                   object.censor,
                   data,
                   treatment,
                   strata,
                   contrasts,
                   times,
                   cause,
                   level.censoring,
                   n.contrasts,
                   levels,
                   n.censor,
                   estimator,
                   eventVar.time,
                   eventVar.status,
                   prob.event,
                   prob.treatment,
                   prob.censor,
                   augTerm,
                   iid,
                   store.iid,
                   export,
                   ...){

    names(list(...))
    n.obs <- NROW(data)
    n.contrasts <- length(contrasts)
    n.times <- length(times)
    iidTotal <- iid
    
    ## ** Compute influence function for each modality
    
    for(iC in 1:n.contrasts){

        if(!is.null(treatment)){
            data.i <- data
            data.i[[treatment]] <- factor(contrasts[iC], levels = levels)
            n.strata <- n.obs
        }else{
            index.strata <- which(data[[strata]]==contrasts[iC])
            n.strata <- length(index.strata)
            data.i <- data[index.strata]
        }
        
        if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
            iW.IPTW <- (data[[treatment]] == contrasts[iC])/prob.treatment[,iC]
        }
        
        
        ## *** outcome model
        if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
            
            if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                factor <- cbind(1-iW.IPTW)
            }else{
                factor <- matrix(1, nrow = n.obs, ncol = 1)
            }

            iid.outcome <- predictRiskIID(model.event,
                                          newdata = data.i,
                                          times = times,
                                          average.iid = TRUE,
                                          factor = factor,
                                          cause = cause)

            iidTotal[[contrasts[iC]]] <- iidTotal[[contrasts[iC]]] + iid.outcome[[1]]
            ## dim(iid.outcome[[1]])
            ## dim(iidTotal)
        }

        
        ## *** Inverse probability weighting of treatment       
        if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
            browser()
            
            predictRiskIID(object, newdata, average.iid, factor)
            factor <- cbind("IPW0" = data[, .SD$weights * .SD$treatment.bin0 * (.SD$status.tau==1) / .SD$prob.treatment0^2],
                            "IPW1" = data[, .SD$weights * .SD$treatment.bin1 * (.SD$status.tau==1) / .SD$prob.treatment1^2],
                            "AIPW0" = data[, .SD$treatment.bin0 * .SD$prob.event0 / .SD$prob.treatment0^2],
                            "AIPW1" = data[, .SD$treatment.bin1 * .SD$prob.event1 / .SD$prob.treatment1^2])
        
        
            prediction.treatment.iid <- predictRiskIID(object.treatment, newdata = data.i, average.iid = average.iid)
            iidIPW.treatment0 <-  prediction.treatment.iid[,1]
            iidIPW.treatment1 <- -prediction.treatment.iid[,2]
            iidAIPW.treatment0 <- -prediction.treatment.iid[,3]
            iidAIPW.treatment1 <- prediction.treatment.iid[,4]
        
            for(iTau in 1:n.times){ ## iTau <- 1
                ## random variables stopped at times
                status.tau <- (data.i[[eventVar.time]] <= times[iTau]) * (data.i[[eventVar.time]] == cause)

                ## ** IPTW,IPCW
                ## compute IPTW
                iProb.treatment <- predictRisk(object.treatment, newdata = data)
                iW.IPTW <- (data[[treatment]] == contrasts[iC]) / iProb.treatment 
                if(return.iid){
                    ls.prob.treatment[[iC]][,iTau] <- iProb.treatment
                }
            
                ## compute IPCW
                if(n.censor[iTau]==0){
                    iW.IPCW <- rep(1,NROW(data))
                }else{
                    Ncensoring.tau <- (data.i[[eventVar.time]] <= times[iTau]) * (data.i[[eventVar.time]] != level.censoring)
                    time.tau <- pmin(data.i[[eventVar.time]], times[iTau])

                    iProb.censor <- predictCox(object.censor, newdata = data, times = time.tau-(1e-10), type = "survival", diag = TRUE)$survival[,1]
                    iW.IPCW <- Ncensoring.tau / iProb.censor
                    if(return.iid){
                        ls.prob.censor[[iC]][,iTau] <- iProb.event
                    }
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
    }

    return(iidTotal)            
}
     
######################################################################
### ate-iid.R ends here
