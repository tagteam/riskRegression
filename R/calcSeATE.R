### calcSeATE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: jul  6 2018 (09:23) 
##           By: Brice Ozenne
##     Update #: 277
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

calcSeATE <- function(object, data, times, cause,
                      treatment, contrasts, strata, n.contrasts, levels, n.times, n.obs,
                      pointEstimate, export, store.iid){

    lowerBand <- upperBand <- diffBand.lower <- diffBand.upper <- ratioBand.lower <- ratioBand.upper <- NULL ## [:forCRANcheck:]
    Treatment <- meanRisk <- . <- Treatment.A <- Treatment.B <- .GRP <- NULL ## [:forCRANcheck:]
    lower <- upper <- diff.se <- diff.lower <- diff.upper <- diff.p.value <- ratio.se <- ratio.lower <- ratio.upper <- ratio.p.value <- NULL ## [:forCRANcheck:]

    out <- list()
    
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
            risk.i <- .predictGLM(object, newdata = data.i)
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

        iIndexA <- intersect(which(pointEstimate$meanRisk[[1]]==iTreatmentA),
                             which(pointEstimate$meanRisk[[2]]==iTime))
        iIndexB <- intersect(which(pointEstimate$meanRisk[[1]]==iTreatmentB),
                             which(pointEstimate$meanRisk[[2]]==iTime))

        ## Compute the iid function of the average treatment effect (difference)
        out$diffRisk.iid[iContrast,] <- out$meanRisk.iid[iIndexA,] - out$meanRisk.iid[iIndexB,]
          
        ## Compute the iid function of the average treatment effect (ratio)
        ## IF(A/B) = IF(A)/B-IF(B)A/B^2
        term1 <- out$meanRisk.iid[iIndexA,] / pointEstimate$meanRisk[["meanRisk"]][iIndexB]
        term2 <- out$meanRisk.iid[iIndexB,] * pointEstimate$meanRisk[["meanRisk"]][iIndexA] / pointEstimate$meanRisk[["meanRisk"]][iIndexB]^2
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
