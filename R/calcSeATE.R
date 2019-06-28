### calcSeATE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: jun 28 2019 (16:02) 
##           By: Brice Ozenne
##     Update #: 325
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * calcSeATE_TI
calcSeATE_TI <- function(object.event,
                         object.treatment,
                         object.censor,
                         data,
                         treatment,
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

    ## ** Compute influence function for each modality
    for(iC in 1:n.contrasts){
        data.i <- data
        data.i[[treatment]] <- factor(contrasts[iC], levels = levels)

        if(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW")){
            iW.IPTW <- sweep(1/prob.treatment[[iC]], FUN = "*", MARGIN = 1, STATS = (data[[treatment]] == contrasts[iC]))
        }
        
        ## *** outcome model
        if(estimator %in% c("Gformula","AIPTW","AIPTW,AIPCW")){
            if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                factor <- cbind(1, 1-iW.IPTW)
            }else{
                factor <- NULL
            }

            iid.outcome <- predictRiskIID(model.event,
                                          newdata = data.i,
                                          times = times,
                                          iid = FALSE,
                                          average.iid = TRUE,
                                          factor = factor,
                                          cause = cause)

            if(estimator %in% c("AIPTW","AIPTW,AIPCW")){
                iidG.event <- iid.outcome[[1]]
                ls.iidAIPW.event <- iid.outcome[-1]
            }else{
                iidG.event <- iid.outcome$absRisk.average.iid
            }
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

        }
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

            meanRisk[iC,] <- colMeans(iIF)
            if(return.iid){
                ls.iid[[iC]] <- iIF
            }
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

    out <- list(meanRisk = data.table(Treatment=name.strata,
                                      time = times,
                                      meanRisk=as.double(meanRisk)),
                riskComparison = riskComparison)
    if(estimator %in% c("AIPTW","AIPTW,AIPCW")){    
        attr(out,"iid") <- ls.iid
        attr(out,"augTerm") <- ls.augTerm
        attr(out,"augTerm") <- ls.augTerm
        attr(out,"augTerm") <- ls.augTerm
    }
    return(out)            
}


## * calcSeATE
calcSeATE <- function(object.event, object.treatment,
                      iid, data, times, cause,
                      treatment, contrasts, strata, levels, n.obs,
                      pointEstimate, export, store.iid){

    lowerBand <- upperBand <- diffBand.lower <- diffBand.upper <- ratioBand.lower <- ratioBand.upper <- NULL ## [:forCRANcheck:]
    Treatment <- meanRisk <- . <- Treatment.A <- Treatment.B <- .GRP <- NULL ## [:forCRANcheck:]
    lower <- upper <- diff.se <- diff.lower <- diff.upper <- diff.p.value <- ratio.se <- ratio.lower <- ratio.upper <- ratio.p.value <- NULL ## [:forCRANcheck:]

    out <- list()
    n.contrasts <- length(contrasts)
    n.times <- length(times)
    browser()
    
                                        # {{{ 1- influence function for the individual predictions
    ## in hypothetical worlds in which every subject is treated with the same treatment
    if(!is.null(treatment)){        
        n.iid <- n.contrasts
        average.iid <- TRUE
        iid.argPredict <- FALSE
    }else{
        n.iid <- 1
        average.iid <- FALSE
        iid.argPredict <- TRUE
    }
    IFrisk <- lapply(1:n.iid,function(i){
        data.i <- data
        if(!is.null(treatment)){
            data.i[[treatment]] <- factor(contrasts[i], levels = levels)
        }
        ## influence function for the absolute risk
        if ("CauseSpecificCox" %in% class(object)){
            pred.i <- predict(object,
                              newdata = data.i,
                              times = times,
                              cause = cause,
                              se = FALSE,
                              iid = iid.argPredict,
                              keep.times = FALSE,
                              store.iid = store.iid,
                              average.iid = average.iid)
            risk.i <- pred.i$absRisk
            if(!is.null(treatment)){
                attr(risk.i,"iid") <- pred.i[["absRisk.average.iid"]]
            }else{
                attr(risk.i,"iid") <- pred.i[["absRisk.iid"]]
            }
        } else if(any(c("coxph","cph") %in% class(object))){
            pred.i <- predictCox(object,
                                 newdata = data.i,
                                 times = times,
                                 se = FALSE,
                                 iid = iid.argPredict,
                                 keep.times = FALSE,
                                 type = "survival",
                                 store.iid = store.iid,
                                 average.iid = average.iid)
            risk.i <- 1-pred.i$survival
            if(!is.null(treatment)){
                attr(risk.i,"iid") <- -pred.i[["survival.average.iid"]]
            }else{
                attr(risk.i,"iid") <- -pred.i[["survival.iid"]]
            }
        }else if("glm" %in% class(object)){
            risk.i <- predictGLM(object, newdata = data.i)
            if(average.iid){
                attr(risk.i,"iid") <- colMeans(attr(risk.i,"iid"))
            }else{
                attr(risk.i,"iid") <- array(attr(risk.i,"iid"), dim = c(NROW(attr(risk.i,"iid")),1,NCOL(attr(risk.i,"iid"))))
            }
            ## se.pred sqrt(colSums(iid.pred^2))
        }
        return(risk.i)
    })

    if(is.null(treatment)){
        IFrisk <- lapply(1:n.contrasts, function(iC){
            ## iid [pred,time,train]
            indexC <- which(data[[strata]]==contrasts[iC])
            iOut <- IFrisk[[1]][indexC,,drop=FALSE]
            ## for each time and initial sample average over the levels of the covariates
            attr(iOut,"iid") <- apply(attr(IFrisk[[1]],"iid")[indexC,,,drop=FALSE], MARGIN = 2:3, FUN = mean)
            return(iOut)
        })        
    }
    
                                        # }}}
    
                                        # {{{ 2- influence function for the average treatment effect
    ## IF had dimension n.predictions (row), n.times (columns), n.dataTrain (length)
    name.treatmentTime <- paste0(pointEstimate$meanRisk[[1]],".",pointEstimate$meanRisk[[2]])
    n.treatmentTime <- length(name.treatmentTime)

    out$meanRisk.iid <- matrix(NA, nrow = n.treatmentTime, ncol = n.obs,
                               dimnames = list(name.treatmentTime, NULL))
    if("se" %in% export){
        out$meanRisk.se <- matrix(NA, nrow = n.treatmentTime, ncol = 1,
                                  dimnames = list(name.treatmentTime, NULL))
    }
    
    for(iTreat in 1:n.contrasts){ # iTreat <- 1
        term1 <- t(attr(IFrisk[[iTreat]],"iid"))
        ## note: here IFrisk[[iTreat]] is the risk (the influence function is in the attribute "iid")
        if(is.null(treatment)){
            indexC <- which(data[[strata]]==contrasts[iTreat])
            term2full <- matrix(-pointEstimate$meanRisk[Treatment==contrasts[iTreat],meanRisk], ncol = n.times, nrow = n.obs,
                                byrow = TRUE)
            term2full[indexC,] <- term2full[indexC,] + IFrisk[[iTreat]] * n.obs/length(indexC)
            iid.tempo <- term1 + term2full/n.obs
        }else{
            term2 <- rowCenter_cpp(IFrisk[[iTreat]], center = pointEstimate$meanRisk[Treatment==contrasts[iTreat],meanRisk])
            ## we get n * IF instead of IF for the absolute risk. This is why the second term need to be rescaled
            iid.tempo <- term1 + t(term2)/n.obs
        }
        out$meanRisk.iid[(iTreat-1)*n.times + 1:n.times,] <- iid.tempo        
    }
    if("se" %in% export){
        out$meanRisk.se[] <- sqrt(rowSums(out$meanRisk.iid^2))
    }
                                        # }}}

                                        # {{{ 3- influence function for the difference/ratio in average treatment effect
    name.T2Time <- paste0(pointEstimate$riskComparison[[1]],".",pointEstimate$riskComparison[[2]],".",pointEstimate$riskComparison[[3]])
    n.T2Time <- length(name.T2Time)

    out$diffRisk.iid <- matrix(NA, nrow = n.T2Time, ncol = n.obs,
                               dimnames = list(name.T2Time, NULL))
    out$ratioRisk.iid <- matrix(NA, nrow = n.T2Time, ncol = n.obs,
                                dimnames = list(name.T2Time, NULL))
    if("se" %in% export){
        out$diffRisk.se <- matrix(NA, nrow = n.T2Time, ncol = 1,
                                  dimnames = list(name.T2Time, NULL))
        out$ratioRisk.se <- matrix(NA, nrow = n.T2Time, ncol = 1,
                                   dimnames = list(name.T2Time, NULL))
    }
    
    for(iContrast in 1:n.T2Time){ ## iContrast <- 1

        iTreatmentA <- pointEstimate$riskComparison[[1]][iContrast]
        iTreatmentB <- pointEstimate$riskComparison[[2]][iContrast]
        iTime <- pointEstimate$riskComparison[[3]][iContrast]

        ## use identical instead of which to be compatible with NA
        iIndexTime <- which(sapply(pointEstimate$meanRisk[[2]],identical,iTime))
        iIndexA <- intersect(which(pointEstimate$meanRisk[[1]]==iTreatmentA), iIndexTime)
        iIndexB <- intersect(which(pointEstimate$meanRisk[[1]]==iTreatmentB), iIndexTime)

        ## Compute the iid function of the average treatment effect (difference)
        out$diffRisk.iid[iContrast,] <- out$meanRisk.iid[iIndexB,] - out$meanRisk.iid[iIndexA,]
          
        ## Compute the iid function of the average treatment effect (ratio)
        ## IF(A/B) = IF(A)/B-IF(B)A/B^2
        term1 <- out$meanRisk.iid[iIndexB,] / pointEstimate$meanRisk[["meanRisk"]][iIndexA]
        term2 <- out$meanRisk.iid[iIndexA,] * pointEstimate$meanRisk[["meanRisk"]][iIndexB] / pointEstimate$meanRisk[["meanRisk"]][iIndexA]^2
        out$ratioRisk.iid[iContrast,] <- term1 - term2
        
    }
    if("se" %in% export){
        out$diffRisk.se[] <- sqrt(rowSums(out$diffRisk.iid^2))
        out$ratioRisk.se[] <- sqrt(rowSums(out$ratioRisk.iid^2))
    }
                                        # }}}
    ## export
    if("iid" %in% export == FALSE){
        out$meanRisk.iid <- NULL
        out$diffRisk.iid <- NULL
        out$ratiAte.iid <- NULL
    }
    return(out)
    
}
     
######################################################################
### calcSeATE.R ends here
