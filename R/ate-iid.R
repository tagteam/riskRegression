### ate-iid.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: jul  9 2019 (09:40) 
##           By: Brice Ozenne
##     Update #: 430
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
                   status.tau,
                   Ncensoring.tau,
                   ...){

    n.obs <- NROW(data)
    n.contrasts <- length(contrasts)
    n.times <- length(times)
    iidTotal <- iid

    ## ** Compute influence function relative to each model
    
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
        if(estimator %in% c("IPTW,IPCW","AIPTW,AIPCW")){
            iW.IPCW <- Ncensoring.tau / prob.censor
        }else if(estimator %in% c("IPTW","AIPTW")){
            iW.IPCW <- matrix(1, nrow = n.obs, ncol = n.times)
        }
        
        ## *** outcome model
        ## terms relative to the augmentation term are missing (neglected)
        if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
            
            if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                factor <- cbind(1-iW.IPTW)
            }else if(estimator %in% "Gformula"){
                factor <- matrix(1, nrow = n.strata, ncol = 1)
            }
            iid.outcome <- predictRiskIID(object.event,
                                          newdata = data.i,
                                          times = times,
                                          average.iid = TRUE,
                                          factor = factor,
                                          cause = cause)

            iidTotal[[contrasts[iC]]] <- iidTotal[[contrasts[iC]]] + iid.outcome[[1]]
            ## dim(iid.outcome[[1]])
            ## dim(iidTotal)
        }

        
        ## *** Treatment model
        if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){

            if(estimator %in% c("IPTW","IPTW,IPCW")){
                factor <- colMultiply_cpp(status.tau * iW.IPCW, scale = -iW.IPTW/prob.treatment[,iC])
            }else if(estimator %in% c("AIPTW")){
                factor <- colMultiply_cpp(status.tau * iW.IPCW - prob.event[[iC]], scale = -iW.IPTW/prob.treatment[,iC])
            }else if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                factor <- colMultiply_cpp(status.tau * iW.IPCW - prob.event[[iC]] + augTerm, scale = -iW.IPTW/prob.treatment[,iC])
            }
            
            iid.treatment <- predictRiskIID(object.treatment,
                                            newdata = data,
                                            average.iid = TRUE,
                                            factor = factor,
                                            level = contrasts[iC])
            
            iidTotal[[contrasts[iC]]] <- iidTotal[[contrasts[iC]]] + do.call(cbind,iid.treatment)
        }

        ## *** censoring model
        ## neglected

        colnames(iidTotal[[iC]]) <- times
    }
    return(iidTotal)            
}
     
######################################################################
### ate-iid.R ends here
