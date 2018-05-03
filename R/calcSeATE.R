### calcSeATE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (17:01) 
## Version: 
## Last-Updated: maj  3 2018 (18:12) 
##           By: Brice Ozenne
##     Update #: 190
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
                      pointEstimate,
                      alpha, conf.level,
                      se, band, nsim.band, store.iid){

    lowerBand <- upperBand <- diffBand.lower <- diffBand.upper <- ratioBand.lower <- ratioBand.upper <- NULL ## [:forCRANcheck:]
    Treatment <- meanRisk <- . <- Treatment.A <- Treatment.B <- .GRP <- NULL ## [:forCRANcheck:]
    lower <- upper <- diff.se <- diff.lower <- diff.upper <- diff.p.value <- ratio.se <- ratio.lower <- ratio.upper <- ratio.p.value <- NULL ## [:forCRANcheck:]
    
                                        # {{{ 1- influence function for the individual predictions
    ## in hypothetical worlds in which every subject is treated with the same treatment
    if(!is.null(treatment)){
        n.iid <- n.contrasts
        average.iid <- TRUE
        iid <- FALSE
    }else{
        n.iid <- 1
        average.iid <- FALSE
        iid <- TRUE
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
                              cause=cause,
                              se = FALSE,
                              iid = iid,
                              keep.times = FALSE,
                              log.transform = FALSE,
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
                                 iid = iid,
                                 keep.times = FALSE,
                                 log.transform = FALSE,
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
            risk.i <- cbind(predict(object, type = "response", newdata = data.i, se=FALSE))
                
            ## compute influence function of the coefficients using lava
            iid.beta <- lava::iid(object)
            newX <- model.matrix(object$formula, data.i)
            if(object$family$link=="logit"){
                ## 1/(1+exp(-Xbeta)) - risk.i
                ## newX %*% coef(object) - Xbeta
                Xbeta <- predict(object, type = "link", newdata = data.i, se=FALSE)
                iid.pred <- sapply(1:n.obs, function(iObs){ ## iObs <- 1
                    iid.beta %*% cbind(newX[iObs,]) * exp(-Xbeta[iObs])/(1+exp(-Xbeta[iObs]))^2
                })                    
            }else if(object$family$link=="identity"){
                iid.pred <- apply(newX, 1, function(iRow){ ## iRow <- newX[1,]
                    iid.beta %*% cbind(iRow)
                })
            }else {
                stop("Cannot handle ",object$family$link," \n",
                     "Only handle the following link function: identity, logit \n")
            }
            if(average.iid){
                attr(risk.i,"iid") <- rowMeans(iid.pred)
            }else{
                attr(risk.i,"iid") <- array(iid.pred, dim = c(NROW(iid.pred),1,NCOL(iid.pred)))
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
    iid.treatment <- array(NA, dim = c(n.contrasts, n.times, n.obs))
    sdIF.treatment <- matrix(NA, nrow = n.contrasts, ncol = n.times)
    for(iTreat in 1:n.contrasts){ # iTreat <- 1
        term1 <- t(attr(IFrisk[[iTreat]],"iid"))
        ## note: here IFrisk[[iTreat]] is the risk (the influence function is in the attribute "iid")
        term2 <- rowCenter_cpp(IFrisk[[iTreat]], center = pointEstimate$meanRisk[Treatment==contrasts[iTreat],meanRisk])
        
        ## we get n * IF instead of IF for the absolute risk. This is why the second term need to be rescaled
        if(is.null(treatment)){
            indexC <- which(data[[strata]]==contrasts[iTreat])
            term2full <- matrix(0, ncol = n.times, nrow = n.obs) 
            term2full[indexC,] <- term2/length(indexC) 
            iid.treatment[iTreat,,] <- term1 + term2full
        }else{
            iid.treatment[iTreat,,] <- term1 + t(term2)/n.obs
        }
        sdIF.treatment[iTreat,] <- apply(iid.treatment[iTreat,,,drop=FALSE],2, ## MARGIN=2 and drop=FALSE to deal with the case of one timepoint
                                         function(x){sqrt(sum(x^2))}
                                         )
    }
                                        # }}}

                                        # {{{ 3- influence function for the difference/ratio in average treatment effect
    nall.contrasts <- n.contrasts*(n.contrasts-1)/2
    iid_diff.contrasts <- array(NA, dim = c(nall.contrasts, n.times, n.obs))
    sdIF_diff.contrasts <- matrix(NA, nrow = nall.contrasts, ncol = n.times)
    iid_ratio.contrasts <- array(NA, dim = c(nall.contrasts, n.times, n.obs))
    sdIF_ratio.contrasts <- matrix(NA, nrow = nall.contrasts, ncol = n.times)
    iiCon <- 0
    sdIF.fct <- pointEstimate$riskComparison[,.(Treatment.A,Treatment.B,time)]
    if(se){
        sdIF.fct[,c("diff.se","ratio.se") := as.double(NA)]
    }
    if(band){
        sdIF.fct[,c("diffBand.quantile","ratioBand.quantile") := as.double(NA)]
    }

    for(iCon in 1:((n.contrasts-1))){ # iCon <- 1
        for(iCon2 in (iCon+1):n.contrasts){ # iCon2 <- 2
            iiCon <- iiCon + 1 ## index of which comparison is performed - used to store the results

            ## Compute the iid function of the average treatment effect (difference)
            iid_diff.contrasts[iiCon,,] <- iid.treatment[iCon,,] - iid.treatment[iCon2,,]
            sdIF_diff.contrasts[iiCon,] <- apply(iid_diff.contrasts[iiCon,,,drop=FALSE],2,
                                                 function(x){sqrt(sum(x^2))}
                                                 )
          
            ## IF(A/B) = IF(A)/B-IF(B)A/B^2
            ate.iCon <- pointEstimate$meanRisk[Treatment == contrasts[iCon],meanRisk]
            ate.iCon2 <- pointEstimate$meanRisk[Treatment == contrasts[iCon2],meanRisk]

            term1 <- sweep(iid.treatment[iCon,,,drop=FALSE], MARGIN = 2:3, STATS = ate.iCon2, FUN = "/")
            term2 <- sweep(iid.treatment[iCon2,,,drop=FALSE], MARGIN = 2:3, STATS = ate.iCon / ate.iCon2^2, FUN = "*")
                    
            iid_ratio.contrasts[iiCon,,] <- term1 - term2
            sdIF_ratio.contrasts[iiCon,] <- apply(iid_ratio.contrasts[iiCon,,,drop=FALSE],2,
                                                  function(x){sqrt(sum(x^2))}
                                                  )
                
            ## store the result
            if(se){
                sdIF.fct[Treatment.A==contrasts[[iCon]] & Treatment.B==contrasts[[iCon2]],
                         c("diff.se","ratio.se") := .(sdIF_diff.contrasts[iiCon,],sdIF_ratio.contrasts[iiCon,])]
            }                   
        }
    }
                                        # }}}

                                        # {{{ 4- compute quantiles for the confidence bands
    if(band){ # nsim.band <- 500
            quantileIF <- confBandCox(iid = abind::abind(iid.treatment, iid_diff.contrasts, iid_ratio.contrasts, along = 1),
                                      se = rbind(sdIF.treatment, sdIF_diff.contrasts, sdIF_ratio.contrasts),
                                      n.sim = nsim.band,
                                      conf.level = conf.level)
        
            qIF.treatment <- quantileIF[1:n.contrasts]                
            sdIF.fct[,c("diffBand.quantile","ratioBand.quantile") := .(quantileIF[n.contrasts+.GRP],
                                                                       quantileIF[n.contrasts+nall.contrasts+.GRP]),
                     by = c("Treatment.A","Treatment.B")]
        
        
    }
                                        # }}}
    
                                        # {{{ 5- compute confidence intervals and confidence bands
    crisks <- sdIF.fct[,.(Treatment.A,Treatment.B,time)]
    mrisks <- data.table::data.table(Treatment = pointEstimate$meanRisk$Treatment,
                                     time = times)
      
    mrisks[, meanRisk := pointEstimate$meanRisk$meanRisk]
    if(se){
        mrisks[, se := sdIF.treatment[.GRP,], by = "Treatment"]            
        mrisks[, lower := meanRisk + qnorm(alpha/2) * sdIF.treatment[.GRP,], by = "Treatment"]
        mrisks[, upper := meanRisk + qnorm(1-alpha/2) * sdIF.treatment[.GRP,], by = "Treatment"]

        crisks[, diff.se := sdIF.fct$diff.se]
        crisks[, diff.lower := pointEstimate$riskComparison$diff - qnorm(1-alpha/2) * sdIF.fct$diff.se]
        crisks[, diff.upper := pointEstimate$riskComparison$diff + qnorm(1-alpha/2) * sdIF.fct$diff.se]
        crisks[, diff.p.value := 2*(1-pnorm(abs(pointEstimate$riskComparison$diff), sd = sdIF.fct$diff.se))]

        crisks[, ratio.se := sdIF.fct$ratio.se]
        crisks[, ratio.lower := pointEstimate$riskComparison$ratio - qnorm(1-alpha/2) * sdIF.fct$ratio.se]
        crisks[, ratio.upper := pointEstimate$riskComparison$ratio + qnorm(1-alpha/2) * sdIF.fct$ratio.se]
        crisks[, ratio.p.value := 2*(1-pnorm(abs(pointEstimate$riskComparison$ratio-1), sd = sdIF.fct$ratio.se))]                    
      }
        if(band){
            mrisks[, lowerBand := meanRisk - qIF.treatment[.GRP] * sdIF.treatment[.GRP,], by = "Treatment"]
            mrisks[, upperBand := meanRisk + qIF.treatment[.GRP] * sdIF.treatment[.GRP,], by = "Treatment"]
        
            crisks[, diffBand.lower := pointEstimate$riskComparison$diff - sdIF.fct$diffBand.quantile * sdIF.fct$diff.se]
            crisks[, diffBand.upper := pointEstimate$riskComparison$diff + sdIF.fct$diffBand.quantile * sdIF.fct$diff.se]
        
            crisks[, ratioBand.lower := pointEstimate$riskComparison$ratio - sdIF.fct$ratioBand.quantile * sdIF.fct$ratio.se]
            crisks[, ratioBand.upper := pointEstimate$riskComparison$ratio + sdIF.fct$ratioBand.quantile * sdIF.fct$ratio.se]
        }
      
    mrisks[, meanRisk := NULL]
                                        # }}}

    ## export
    return(list(mrisks = mrisks,
                crisks = crisks))
    
}
     
######################################################################
### calcSeATE.R ends here
