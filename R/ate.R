### ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 23 2016 (08:53) 
## Version: 
## last-updated: Jun 29 2017 (18:01) 
##           By: Thomas Alexander Gerds
##     Update #: 245
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Compute the average treatment effect using CSC.
#'
#' @description Use the g-formula to estimate the average treatment effect
#' @param object outcome model which describes how event risk depends on treatment and covariates.
#' The object carry its own call and have a \code{predictRisk} method. See examples.
#' @param data data set in which to evaluate risk predictions based on the outcome model
#' @param treatment name of the treatment variable
#' @param contrasts the levels of treatment variable to be compared
#' @param times time points at which to evaluate risks
#' @param cause the cause of interest
#' @param conf.level Numeric. Confidence level of the confidence intervals.
#' @param se Logical. If \code{TRUE} add the standard errors and confidence intervals
#' to the output.
#' @param B the number of bootstrap replications used to compute the
#' confidence intervals. If it equals 0, then Wald-type confidence intervals are computed.
#' They rely on the standard error estimated using the influence function of the estimator.
#' @param band Logical. If \code{TRUE} add the confidence bands to the output.
#' @param nSim.band the number of simulations used to compute the quantiles for the confidence bands.
#' @param seed An integer used to generate seeds for bootstrap and to
#'     achieve reproducibility of the bootstrap confidence intervals.
#' @param handler parallel handler for bootstrap. Either "mclapply" or
#'     "foreach". If "foreach" use \code{doParallel} to create a cluster.
#' @param mc.cores Passed to \code{parallel::mclapply} or
#'     \code{doParallel::registerDoParallel}. The number of cores to use, i.e. at
#'     most how many child processes will be run simultaneously.  The
#'     option is initialized from environment variable MC_CORES if
#'     set.
#' @param verbose Logical. If \code{TRUE} inform about estimated run
#'     time.
#' @param logTransform Should the confidence interval for the ratio be computed using a log-tranformation. Only active if Wald-type confidence intervals are computed.
#' @param store.iid Implementation used to estimate the standard error. Can be \code{"full"} or \code{"minimal"}.
#' \code{"minimal"} requires less memory but can only estimate the standard for the difference between treatment effects (and not for the ratio).
#' @param ... passed to predictRisk
#' @author Brice Ozenne \email{broz@@sund.ku.dk} and Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @return A list with: point estimates, bootstrap quantile confidence
#'     intervals model: the CSC model (optional)
#' 
#' @details WARNING: the p.value and the confidence intervals for the ratio using Wald-type approximations are still experimental.
#' 
#' @examples 
#' library(survival)
#' library(rms)
#' 
#' set.seed(10)
#' n <- 100
#' 
#' ## Cox model
#' dtS <- sampleData(n,outcome="survival")
#' dtS$time <- round(dtS$time,1)
#' dtS$X1 <- factor(rbinom(n, prob = c(0.3,0.4) , size = 2), labels = paste0("T",0:2))
#'
#' fit=cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
#'
#' \dontrun{
#' ateFit1 <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
#'         times = 5:8, B = 1e3, y = TRUE,  mc.cores=1)
#' 
#' ateFit1 <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
#'         times = 5:8, B = 1e1, y = TRUE,  mc.cores=1)
#' ateFit2 <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
#'         times = 5:8, B = 0, y = TRUE, band = TRUE, mc.cores=1)
#'
#' ateFit3 <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
#'            times = 5:8, B = 0, y = TRUE, band = TRUE, mc.cores=1,
#'            store.iid = "minimal")
#' }
#' 
#' ## Competing risks: Cause specific Cox regression
#' \dontrun{
#' set.seed(17)
#' n=100
#' dt <- sampleData(n,outcome="competing.risks")
#' dt$time <- round(dt$time,1)
#' dt$X1 <- factor(rbinom(n, prob = c(0.2,0.3,0.2) , size = 3), labels = paste0("T",0:3))
#' fitCR= CSC(Hist(time,event)~ X1+X8,data=dt,cause=1)
#' ate(fitCR, data = dt, treatment = "X1", contrasts = NULL,
#'         times = 7, cause = 1, B = 2, mc.cores=1)
#'
#' atefit=ate(fitCR, data = dt, treatment = "X1", contrasts = NULL,
#'         times = 1:7, cause = 1, mc.cores=1, se = FALSE, band = FALSE)
#'
#' 
#'  ate(fitCR, data = dt, treatment = "X1", contrasts = NULL,
#'         times = 5:7, cause = 1, B = 0, se = TRUE, band = TRUE, mc.cores=1)
#' }
#' 
#' @export
ate <- function(object,
                data,
                treatment,
                contrasts = NULL,
                times,
                cause,
                conf.level = 0.95,
                se = TRUE,
                band = FALSE,
                B = 0,
                nSim.band = 1e3,
                seed,
                handler=c("mclapply","foreach"),
                mc.cores = 1,
                verbose=TRUE,
                logTransform=FALSE,
                store.iid="full",
                ...){
    
    meanRisk=Treatment=ratio=Treatment.A=Treatment.B=b <- NULL
    .=.I <- NULL
    diff.se=ratio.se=.GRP=lower=upper=diff.lower=diff.upper=diff.p.value=ratio.lower=ratio.upper=ratio.p.value <- NULL
    lowerBand=upperBand=diffBand.lower=diffBand.upper=ratioBand.lower=ratioBand.upper <- NULL
    
    
    # {{{ Prepare
    dots <- list(...)
    
    if(band && B>0){
        stop("the confidence bands cannot be computed when using the bootstrap approach \n",
             "set argument \'band\' to FALSE to not compute the confidence bands \n",
             "or set argument \'B\' to 0 to use the influence function instead of the bootstrap\n")
    }
    if(treatment %in% names(data) == FALSE){
        stop("The data set does not seem to have a variable ",treatment," (argument: treatment). \n")
    }
    test.CR <- !missing(cause) # test whether the argument cause has been specified, i.e. it is a competing risk model
    if(test.CR==FALSE){cause <- NA}

    if(B==0 && (se || band)){
        validClass <- c("CauseSpecificCox","coxph","cph","phreg")
        if(all(validClass %in% class(object) == FALSE)){
            stop("Standard error based on the influence function only implemented for \n",
                 paste(validClass, collapse = " ")," objects \n",
                 "set argument \'B\' to a positive integer to use a boostrap instead \n")
        }
    }
    data[[treatment]] <- factor(data[[treatment]])
  
    if(is.null(contrasts)){
        levels <- levels(data[[treatment]])
        contrasts <- levels(data[[treatment]])
        ## if (length(contrasts)>50) warning("Treatment variable has more than 50 levels.\nIf this is not a mistake,
        ## you should use the argument `contrasts'.")
    }else{levels <- contrasts}
    n.contrasts <- length(contrasts)
    n.times <- length(times)
    n.obs <- NROW(data)

    # }}}
    
    # {{{ Checking the model

    # for predictRisk S3-method
    allmethods <- utils::methods(predictRisk)
    candidateMethods <- paste("predictRisk",class(object),sep=".")
    if (all(match(candidateMethods,allmethods,nomatch=0)==0))
        stop(paste("Could not find predictRisk S3-method for ",class(object),collapse=" ,"),sep="")
    # for compatibility with resampling
    if(is.null(object$call))
        stop(paste("The object does not contain its own call, which is needed to refit the model in the bootstrap loop."))

  # }}}
  
    # {{{ calc G formula
  Gformula <- function(object, data, treatment, contrasts, times, cause, ...){
      meanRisk <- lapply(1:n.contrasts,function(i){
          ## prediction for the hypothetical worlds in which every subject is treated with the same treatment
          data.i <- data
          data.i[[treatment]] <- factor(contrasts[i], levels = levels)
          risk.i <- colMeans(do.call("predictRisk",args = list(object, newdata = data.i, times = times, cause = cause, ...)))
      })
      riskComparison <- data.table::rbindlist(lapply(1:(n.contrasts-1),function(i){
          data.table::rbindlist(lapply(((i+1):n.contrasts),function(j){
              ## compute differences between all pairs of treatments
              data.table(Treatment.A=contrasts[[i]],
                         Treatment.B=contrasts[[j]],
                         time = times,
                         diff=meanRisk[[i]]-meanRisk[[j]],
                         ratio=meanRisk[[i]]/meanRisk[[j]])
          }))}))
      name.Treatment <- unlist(lapply(1:n.contrasts, function(c){rep(contrasts[c],length(meanRisk[[c]]))}))
      out <- list(meanRisk = data.table(Treatment=name.Treatment, time = times, meanRisk=unlist(meanRisk)),
                  riskComparison = riskComparison)
      out
  }
    # }}}
    
    # {{{ point estimate
    estimateTime <- system.time(
        pointEstimate <- Gformula(object=object,
                                  data=data,
                                  treatment=treatment,
                                  contrasts=contrasts,
                                  times=times,
                                  cause=cause,
                                  dots))
    # }}}
    
    # {{{ Confidence interval

    ##### Confidence interval
    if(se || band){
    
    alpha <- 1-conf.level
    
        if(B>0){
            # {{{ Bootstrap
            if (verbose==TRUE)
                message(paste0("Approximated bootstrap netto run time (without time for copying data to cores):\n",
                               round(estimateTime["user.self"],2),
                               " seconds times ",
                               B,
                               " bootstraps / ",
                               mc.cores,
                               " cores = ",
                               round(estimateTime["user.self"]*B/mc.cores,2)," seconds.\n"))
            x.cores <- parallel::detectCores()
            if(mc.cores > x.cores){
                warning("Not enough available cores \n",
                        "available: ",parallel::detectCores()," | requested: ",mc.cores,"\n")
                mc.cores=x.cores
            }
            if (!missing(seed)) set.seed(seed)
            bootseeds <- sample(1:1000000,size=B,replace=FALSE)
            if (handler[[1]]=="foreach"){
        
                cl <- parallel::makeCluster(mc.cores)
                doParallel::registerDoParallel(cl)
                pp <- find(as.character(object$call[[1]]))
                addPackage <- if(grep("package:",pp)){gsub("package:","",pp[grep("package:",pp)])}else{NULL}
                boots <- foreach::`%dopar%`(foreach::foreach(b=1:B,.packages=c("riskRegression",addPackage),.export=NULL), {
                    set.seed(bootseeds[[b]])
                    dataBoot <- data[sample(1:n.obs, size = n.obs, replace = TRUE),]
                    object$call$data <- dataBoot
                    objectBoot <- try(eval(object$call),silent=TRUE)
                    if ("try-error" %in% class(objectBoot)){
                        stop(paste0("Failed to fit model ",class(object),ifelse(try(b>0,silent=TRUE),paste0(" in bootstrap step ",b,"."))))
                    }
                    tryCatch(Gformula(object=objectBoot,
                                      data=dataBoot,
                                      treatment=treatment,
                                      contrasts=contrasts,
                                      times=times,
                                      cause=cause,
                                      dots),
                             error = function(x){return(NULL)})
                })
        
                parallel::stopCluster(cl)
        
            } else {
                if(Sys.info()["sysname"] == "Windows" && mc.cores>1){
                    message("mclapply cannot perform parallel computations on Windows \n",
                            "consider setting argument handler to \"foreach\" \n")
                    mc.cores <- 1
                }
        
                boots <- parallel::mclapply(1:B, function(b){
                    set.seed(bootseeds[[b]])
                    dataBoot <- data[sample(1:n.obs, size = n.obs, replace = TRUE),]
                    object$call$data <- dataBoot
                    objectBoot <- try(eval(object$call),silent=TRUE)
                    if ("try-error" %in% class(objectBoot)){
                        stop(paste0("Failed to fit model",ifelse(try(b>0,silent=TRUE),paste0(" in bootstrap step ",b,"."))))
                    }
                    tryCatch(Gformula(object=objectBoot,
                                      data=dataBoot,
                                      treatment=treatment,
                                      contrasts=contrasts,
                                      times=times,
                                      cause=cause,
                                      dots),
                             error = function(x){return(NULL)})
                }, mc.cores = mc.cores)
            }
      
            ## gc()
            meanRisksBoot <- data.table::rbindlist(lapply(boots,function(x)x$meanRisk))
            riskComparisonsBoot <- data.table::rbindlist(lapply(boots,function(x)x$riskComparison))
      
            if(NROW(meanRisksBoot)==0){
                stop("no successful bootstrap \n")
            }

            mrisks <- meanRisksBoot[,data.table::data.table(meanRiskBoot=mean(meanRisk, na.rm = TRUE),
                                                            lower=quantile(meanRisk,alpha/2, na.rm = TRUE),
                                                            upper=quantile(meanRisk,1-(alpha/2), na.rm = TRUE),
                                                            n.boot=sum(!is.na(meanRisk))),
                                    keyby=c("Treatment","time")]
      
            crisks <- riskComparisonsBoot[,data.table::data.table(diffMeanBoot=mean(diff, na.rm = TRUE),
                                                                  diff.lower=quantile(diff,alpha/2, na.rm = TRUE),
                                                                  diff.upper=quantile(diff,1-(alpha/2), na.rm = TRUE),
                                                                  diff.p.value=findP1(diff, alternative = "two.sided"),
                                                                  ratioMeanBoot=mean(ratio, na.rm = TRUE),
                                                                  ratio.lower=quantile(ratio,alpha/2, na.rm = TRUE),
                                                                  ratio.upper=quantile(ratio,1-(alpha/2), na.rm = TRUE),
                                                                  ratio.p.value=findP1(ratio-1, alternative = "two.sided"),
                                                                  n.boot=sum(!is.na(diff))),
                                          keyby=c("Treatment.A","Treatment.B","time")]
            # }}}
        } else {
            
            # {{{ Influence function and variance
            average.iid <- store.iid=="minimal"
            iid <- store.iid!="minimal"
            
            IFrisk <- lapply(1:n.contrasts,function(i){
                #### influence function for the hypothetical worlds in which every subject is treated with the same treatment
                data.i <- data
                data.i[[treatment]] <- factor(contrasts[i], levels = levels)
                ## influence function for the absolute risk
                if ("CauseSpecificCox" %in% class(object)){
                    pred.i <- do.call("predict",args = list(object,
                                                            newdata = data.i,
                                                            times = times,
                                                            cause=cause,
                                                            se=FALSE,
                                                            iid=iid,
                                                            keep.times=FALSE,
                                                            logTransform=FALSE,
                                                            store.iid=store.iid,
                                                            average.iid=average.iid))
                    risk.i <- pred.i$absRisk
                    attr(risk.i,"iid") <- pred.i$absRisk.iid
                } else{
                    pred.i <- do.call("predictCox",args = list(object,
                                                               newdata = data.i,
                                                               times = times,
                                                               se=FALSE,
                                                               iid=iid,
                                                               keep.times=FALSE,
                                                               logTransform=FALSE,
                                                               type="survival",
                                                               store.iid=store.iid,
                                                               average.iid=average.iid))
                    risk.i <- 1-pred.i$survival
                    attr(risk.i,"iid") <- pred.i$survival.iid
                }
                return(risk.i)
            })

            ## influence function for the average treatment effect
            iid.treatment <- array(NA, dim = c(n.contrasts, n.times, n.obs))
            sdIF.treatment <- matrix(NA, nrow = n.contrasts, ncol = n.times)
            for(iTreat in 1:n.contrasts){ # iTreat <- 1
                if(average.iid){
                     term1 <- t(attr(IFrisk[[iTreat]],"iid"))
                }else{
                    term1 <- apply(attr(IFrisk[[iTreat]],"iid"),2:3,mean)
                }

                term2 <- rowCenter_cpp(IFrisk[[iTreat]], center = pointEstimate$meanRisk[Treatment==contrasts[iTreat],meanRisk])

                # we get n * IF instead of IF for the absolute risk. This is why the second term need to be rescaled
                iid.treatment[iTreat,,] <- term1 + t(term2)/sqrt(n.obs)
                sdIF.treatment[iTreat,] <- apply(iid.treatment[iTreat,,,drop=FALSE],2, ## MARGIN=2 and drop=FALSE to deal with the case of one timepoint
                                                 function(x){sqrt(sum(x^2))}
                                                 )
            }

            ## influence function for the difference/ratio in average treatment effect
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
                    ## compute differences between all pairs of treatments
                    # IF had dimension n.predictions (row), n.times (columns), n.dataTrain (length)
                    # apply 2 do for each time: extract the value of IF with n.predictions (row) and n.dataTrain (length)
                    # colSums: compute the IF for the G formula for each observation in the training data set
                    # sqrt(sum 2)) compute the variance of the estimator over the observations in the training data set
                    iiCon <- iiCon + 1
                    if(average.iid){
                        term1 <- t(attr(IFrisk[[iCon]],"iid")-attr(IFrisk[[iCon2]],"iid"))
                    }else{
                        term1 <- apply(attr(IFrisk[[iCon]],"iid")-attr(IFrisk[[iCon2]],"iid"), 2:3, mean)
                    }
                    term2 <- rowCenter_cpp(IFrisk[[iCon]]-IFrisk[[iCon2]],
                                           center = pointEstimate$riskComparison[Treatment.A==contrasts[iCon] & Treatment.B==contrasts[iCon2],diff])
                    iid_diff.contrasts[iiCon,,] <- term1 + t(term2)/sqrt(n.obs)
                    sdIF_diff.contrasts[iiCon,] <- apply(iid_diff.contrasts[iiCon,,,drop=FALSE],2,
                                                         function(x){sqrt(sum(x^2))}
                                                         )

                    # IF(A/B) = IF(A)/B-IF(B)A/B^2
                    if(average.iid){
                        sdIF_ratio.contrasts[iiCon,] <- as.numeric(NA)
                    }else{
                        iidTempo1 <- aperm(attr(IFrisk[[iCon]],"iid"), c(3,2,1))
                        term1 <- aperm(sliceScale_cpp(iidTempo1, IFrisk[[iCon2]]), c(3,2,1))

                    iidTempo2 <- aperm(attr(IFrisk[[iCon2]],"iid"), c(3,2,1))
                    term2 <- aperm(sliceMultiply_cpp(iidTempo2, IFrisk[[iCon]]/IFrisk[[iCon2]]^2), c(3,2,1))

                    term3 <- rowCenter_cpp(IFrisk[[iCon]]/IFrisk[[iCon2]],
                                           center = pointEstimate$riskComparison[Treatment.A==contrasts[iCon] & Treatment.B==contrasts[iCon2],ratio])
                    
                    iid_ratio.contrasts[iiCon,,] <- apply(term1 - term2, 2:3, mean) + t(term3)/sqrt(n.obs)
                    sdIF_ratio.contrasts[iiCon,] <- apply(iid_ratio.contrasts[iiCon,,,drop=FALSE],2,
                                                          function(x){sqrt(sum(x^2))}
                                                          )
                    }
                    ## store the result
                    index.contrasts <- sdIF.fct[, .I[Treatment.A==contrasts[[iCon]] & Treatment.B==contrasts[[iCon2]]]]
                    if(se){
                        sdIF.fct[index.contrasts,diff.se := sdIF_diff.contrasts[iiCon,]]  
                        sdIF.fct[index.contrasts,ratio.se := sdIF_ratio.contrasts[iiCon,]]  
                    }                   
                }
            }            
            # }}}

            # {{{ confidence bands
            if(band){ # nSim.band <- 500

                quantileIF <- confBandCox(iid = abind::abind(iid.treatment, iid_diff.contrasts, iid_ratio.contrasts, along = 1),
                                          se = rbind(sdIF.treatment, sdIF_diff.contrasts, sdIF_ratio.contrasts),
                                          times = times,
                                          n.sim = nSim.band,
                                          conf.level = conf.level)
                
                qIF.treatment <- quantileIF[1:n.contrasts]                
                sdIF.fct[,c("diffBand.quantile","ratioBand.quantile") := .(quantileIF[n.contrasts+.GRP],
                                                                           quantileIF[n.contrasts+nall.contrasts+.GRP]),
                         by = c("Treatment.A","Treatment.B")]

                
            }
            # }}}
            
            # {{{ compute confidence intervals and bands

            crisks <- sdIF.fct[,.(Treatment.A,Treatment.B,time)]
            mrisks <- data.table::data.table(Treatment = pointEstimate$meanRisk$Treatment,
                                             time = times)
            
            mrisks[, meanRisk := pointEstimate$meanRisk$meanRisk]
            if(logTransform){
                stop("not implemented yet \n")
                ##     ratio.lower <- exp(log(pointEstimate$riskComparison$ratio) + qnorm(alpha/2) * sdIF.fct$ratio)
                ##     ratio.upper <- exp(log(pointEstimate$riskComparison$ratio) + qnorm(1-alpha/2) * sdIF.fct$ratio)
                ##     ratio.p.value <- 2*(1-pnorm(abs(log(pointEstimate$riskComparison$ratio)), sd = sdIF.fct$ratio))
            }else{
                if(se){                    
                    mrisks[, lower := meanRisk + qnorm(alpha/2) * sdIF.treatment[.GRP,], by = "Treatment"]
                    mrisks[, upper := meanRisk + qnorm(1-alpha/2) * sdIF.treatment[.GRP,], by = "Treatment"]

                    crisks[, diff.lower := pointEstimate$riskComparison$diff - qnorm(1-alpha/2) * sdIF.fct$diff.se]
                    crisks[, diff.upper := pointEstimate$riskComparison$diff + qnorm(1-alpha/2) * sdIF.fct$diff.se]
                    crisks[, diff.p.value := 2*(1-pnorm(abs(pointEstimate$riskComparison$diff), sd = sdIF.fct$diff.se))]

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
            }            
            mrisks[, meanRisk := NULL]

            # }}}
            
            bootseeds <- NULL
        }} else{
      
             mrisks <- data.table::data.table(Treatment = pointEstimate$meanRisk$Treatment,
                                              time = times)
             crisks <- data.table::data.table(Treatment.A = pointEstimate$riskComparison$Treatment.A,
                                              Treatment.B = pointEstimate$riskComparison$Treatment.B,
                                              time = times)
             bootseeds <- NULL

             # }}}
         }

    ## merge with pointEstimate
    mrisks <- merge(pointEstimate$meanRisk,mrisks,by=c("Treatment","time"))
    crisks <- merge(pointEstimate$riskComparison,crisks,by=c("Treatment.A","Treatment.B","time"))
    out <- list(meanRisk=mrisks,
                riskComparison=crisks,
                treatment=treatment,
                contrasts=contrasts,
                times=times,
                se = se,
                n.bootstrap=B,
                band = band,
                nSim.band = nSim.band,
                seeds=bootseeds,
                conf.level=conf.level,
                logTransform = logTransform)
  
    class(out) <- c("ate",class(object))
    out
}


#' @title Compute the p.value from the distribution under H1
#' @description Compute the p.value from the distribution under H1
#' 
#' @param x the sample
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less"
#' 
#' @examples 
#' set.seed(10)
#' 
#' # no effect
#' x <- rnorm(1e3) 
#' riskRegression:::findP1(x, alternative = "two.sided")
#' riskRegression:::findP1(x, alternative = "greater")
#' riskRegression:::findP1(x, alternative = "less")
#' 
#' # effect
#' x <- rnorm(1e3, mean = 1) 
#' riskRegression:::findP1(x, alternative = "two.sided")
#' riskRegression:::findP1(x, alternative = "greater") # pnorm(q = 0, mean = 1)
#' riskRegression:::findP1(x, alternative = "less")
#' 
#' x <- rnorm(1e3, mean = -1) 
#' riskRegression:::findP1(x, alternative = "two.sided") 
#' riskRegression:::findP1(x, alternative = "greater")
#' riskRegression:::findP1(x, alternative = "less") # pnorm(q = 0, mean = -1)
#' 
findP1 <- function(x, alternative = "two.sided"){ 
  
  x <- na.omit(x)
  if(length(x)==0 || length(x) < 10){
    return(as.numeric(NA))
  }else if(all(x>0)){
    p.value <- switch(alternative,
                      "two.sided" = 0,
                      "less" = 1,
                      "greater" = 1)
  } else if(all(x<0)){
    p.value <- switch(alternative,
                      "two.sided" = 0,
                      "less" = 0,
                      "greater" = 1)
  }else if(all(x==0)){
    p.value <- switch(alternative,
                      "two.sided" = 0,
                      "less" = 0,
                      "greater" = 0)
  }else{
    fn <- function(p){abs(quantile(x, probs = p))}
    
    optimum <- optim(fn = fn, par = 0.5, lower = 0, upper = 1, method = "L-BFGS-B")
    
    
    p.value <- switch(alternative,
                      "two.sided" = if(optimum$par < 0.5){2*optimum$par}else{2*(1-optimum$par)},
                      "less" = 1-optimum$par,
                      "greater" = optimum$par)
    
  }
  
  return(p.value)
}


#----------------------------------------------------------------------
### ate.R ends here
