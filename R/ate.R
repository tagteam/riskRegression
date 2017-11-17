### ate.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 23 2016 (08:53) 
## Version: 
## last-updated: nov 14 2017 (19:34) 
##           By: Brice Ozenne
##     Update #: 379
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @title Compute the average treatment effects via the g-formula
#'
#' @description Use the g-formula to estimate the average treatment
#'     effect based on Cox regression with or without competing risks
#' @param object outcome model which describes how event risk depends
#'     on treatment and covariates.  The object carry its own call and
#'     have a \code{predictRisk} method. See examples.
#' @param data data set in which to evaluate risk predictions based on
#'     the outcome model
#' @param formula For analyses with time-dependent covariates, the response formula. See examples.
#' @param treatment name of the treatment variable
#' @param contrasts the levels of the treatment variable to be
#'     compared
#' @param times time points at which to evaluate risks
#' @param cause the cause of interest
#' @param landmark for models with time-dependent covariates the landmark time(s) of evaluation.
#'        In this case, argument \code{time} may only be one value and for the prediction of risks
#'        it is assumed that that the covariates do not change between landmark and landmark+time.
#' @param conf.level Numeric value between 0 and 1 (default is 0.05). Confidence level of the confidence intervals.
#' @param se Logical. If \code{TRUE} compute standard errors and confidence intervals
#' @param band Logical. If \code{TRUE} compute the confidence bands
#' @param B the number of bootstrap replications used to compute the
#'     confidence intervals. If it equals 0, then Wald-type confidence
#'     intervals are computed.  They rely on the standard error
#'     estimated using the influence function of the estimator.
#' @param nsim.band the number of simulations used to compute the
#'     quantiles for the confidence bands.
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
#' @param store.iid Implementation used to estimate the standard error. Can be \code{"full"} or \code{"minimal"}.
#' \code{"minimal"} requires less memory but can only estimate the standard for the difference between treatment effects (and not for the ratio).
#' @param ... passed to predictRisk
#' @author Brice Ozenne \email{broz@@sund.ku.dk} and Thomas Alexander
#'     Gerds \email{tag@@biostat.ku.dk}
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
#' fit <- cph(formula = Surv(time,event)~ X1+X2,data=dtS,y=TRUE,x=TRUE)
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
#' ## time-dependent covariates
#' \dontrun{
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ celltype+karno + age + trt, veteran)
#' vet2 <- survSplit(Surv(time, status) ~., veteran,
#'                        cut=c(60, 120), episode ="timegroup")
#' fitTD <- coxph(Surv(tstart, time, status) ~ celltype+karno + age + trt,
#'                data= vet2,x=1)
#' set.seed(16)
#' resVet <- ate(fitTD,formula=Hist(entry=tstart,time=time,event=status)~1,
#'           data = vet2, treatment = "celltype", contrasts = NULL,
#'         times=5,verbose=1,
#'         landmark = c(0,30,60,90), cause = 1, B = 20, se = 1,
#'         band = FALSE, mc.cores=1)
#' resVet
#' }
#' \dontrun{
#' set.seed(137)
#' d=sampleDataTD(127)
#' library(survival)
#' d[,status:=1*(event==1)]
#' ## ignore competing risks
#' cox1TD <- coxph(Surv(start,time, status,type="counting") ~ X3+X5+X6+X8, data=d)
#' resTD1 <- ate(cox1TD,formula=Hist(entry=start,time=time,event=status)~1,
#'         data = d, treatment = "X3", contrasts = NULL,
#'         times=.5,verbose=1,
#'         landmark = c(0,0.5,1), B = 20, se = 1,
#'         band = FALSE, mc.cores=1)
#' resTD1
#' ## adjust for competing risks
#' cscTD <- CSC(Hist(time=time, event=event,entry=start) ~ X3+X5+X6+X8, data=d)
#' set.seed(16)
#' resTD <- ate(cscTD,formula=Hist(entry=start,time=time,event=event)~1,
#'         data = d, treatment = "X3", contrasts = NULL,
#'         times=.5,verbose=1,
#'         landmark = c(0,0.5,1), cause = 1, B = 20, se = 1,
#'         band = FALSE, mc.cores=1)
#' resTD
#' }
#' @export
ate <- function(object,
                data,
                formula,
                treatment,
                contrasts = NULL,
                times,
                cause,
                landmark,
                conf.level = 0.95,
                se = TRUE,
                band = FALSE,
                B = 0,
                nsim.band = ifelse(band,1e3,0),
                seed,
                handler="foreach",
                mc.cores = 1,
                verbose=TRUE,
                store.iid="full",
                ...){
  
  meanRisk=Treatment=ratio=Treatment.A=Treatment.B=b <- NULL
  .=.I <- NULL
  diff.se=ratio.se=.GRP=lower=upper=diff.lower=diff.upper=diff.p.value=ratio.lower=ratio.upper=ratio.p.value <- NULL
  lowerBand=upperBand=diffBand.lower=diffBand.upper=ratioBand.lower=ratioBand.upper <- NULL
  
  handler <- match.arg(handler, c("foreach","mclapply"))
  # {{{ checking for time-dependent covariates (left-truncation)
  TD <- switch(class(object)[[1]],"coxph"=(attr(object$y,"type")=="counting"),
               "CauseSpecificCox"=(attr(object$models[[1]]$y,"type")=="counting"),FALSE)
  if (TD){
    if (missing(formula))
      stop("Need formula to do landmark analysis.")
    if (missing(landmark))
      stop("Need landmark time(s) to do landmark analysis.")
    if(length(times)!=1){
      stop("In settings with time-dependent covariates argument 'time' must be a single value, argument 'landmark' may be a vector of time points.")
    }
  }else{
    landmark=NULL
  }
  # }}}
  
    # {{{ Prepare
    dots <- list(...)

    if(se==0 && B>0){
        warning("argument 'se=0' means 'no standard errors' so number of bootstrap repetitions is forced to B=0.")
    }
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
  if (TD){
    Gformula <- function(object,
                         data,
                         treatment,
                         contrasts,
                         times,
                         landmark,
                         cause,
                         ...){
      response <- eval(formula[[2]],envir=data)
      time <- response[,"time"]
      entry <- response[,"entry"]
      if(class(object)[[1]]=="coxph"){
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
                                     args = list(object,
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
          RC[,diff:=meanRisk-dt.meanRisk[Treatment==contrasts[[j]],meanRisk]]
          RC[,ratio:=meanRisk/dt.meanRisk[Treatment==contrasts[[j]],meanRisk]]
          RC[,meanRisk:=NULL]
          RC[]
        }))}))
      out <- list(meanRisk = dt.meanRisk, riskComparison = riskComparison)
      out
    }
  }else{
    Gformula <- function(object,
                         data,
                         treatment,
                         contrasts,
                         times,
                         landmark,
                         cause,
                         ...){
      meanRisk <- lapply(1:n.contrasts,function(i){
        ## prediction for the hypothetical worlds in which every subject is treated with the same treatment
        data.i <- data
        data.i[[treatment]] <- factor(contrasts[i], levels = levels)
        risk.i <- colMeans(do.call("predictRisk",args = list(object,newdata = data.i,times = times,cause = cause,...)))
        risk.i
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
                              landmark=landmark,
                              dots))
  # }}}
  
    # {{{ Confidence interval    
    if(se || band){
        if (TD){
            key1 <- c("Treatment","landmark")
            key2 <- c("Treatment.A","Treatment.B","landmark")
        }
        else{
            key1 <- c("Treatment","time")
            key2 <- c("Treatment.A","Treatment.B","time")
        }
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
        
            if(verbose){
                cl <- parallel::makeCluster(mc.cores, outfile = "")
                pb <- txtProgressBar(max = B, style = 3)          
            }else{
                cl <- parallel::makeCluster(mc.cores)
            }
            doParallel::registerDoParallel(cl)
            pp <- find(as.character(object$call[[1]]))
            addPackage <- if(grep("package:",pp)){gsub("package:","",pp[grep("package:",pp)])}else{NULL}
        
        boots <- foreach::`%dopar%`(foreach::foreach(b=1:B,.packages=unique(c("riskRegression","survival",addPackage)),
                                                     .export=NULL), {
          set.seed(bootseeds[[b]])
          if(verbose){setTxtProgressBar(pb, b)}
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
                            landmark=landmark,
                            dots),
                   error = function(x){return(NULL)})
        })
        if(verbose){close(pb)}
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
                            landmark=landmark,
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
                              keyby=key1]
      crisks <- riskComparisonsBoot[,data.table::data.table(diffMeanBoot=mean(diff, na.rm = TRUE),
                                                            diff.lower=quantile(diff,alpha/2, na.rm = TRUE),
                                                            diff.upper=quantile(diff,1-(alpha/2), na.rm = TRUE),
                                                            diff.p.value=findP1(diff, alternative = "two.sided"),
                                                            ratioMeanBoot=mean(ratio, na.rm = TRUE),
                                                            ratio.lower=quantile(ratio,alpha/2, na.rm = TRUE),
                                                            ratio.upper=quantile(ratio,1-(alpha/2), na.rm = TRUE),
                                                            ratio.p.value=findP1(ratio-1, alternative = "two.sided"),
                                                            n.boot=sum(!is.na(diff))),
                                    keyby=key2]
      ## merge with pointEstimate
      mrisks <- merge(pointEstimate$meanRisk,mrisks,by=key1)
      crisks <- merge(pointEstimate$riskComparison,crisks,by=key2)
      # }}}
    } else {
      
      # {{{ Influence function and variance
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
                                                  iid=FALSE,
                                                  keep.times=FALSE,
                                                  log.transform=FALSE,
                                                  store.iid=store.iid,
                                                  average.iid=TRUE))
          risk.i <- pred.i$absRisk
          attr(risk.i,"iid") <- pred.i$absRisk.average.iid
        } else{
          pred.i <- do.call("predictCox",args = list(object,
                                                     newdata = data.i,
                                                     times = times,
                                                     se=FALSE,
                                                     iid=FALSE,
                                                     keep.times=FALSE,
                                                     log.transform=FALSE,
                                                     type="survival",
                                                     store.iid=store.iid,
                                                     average.iid=TRUE))
          risk.i <- 1-pred.i$survival
          attr(risk.i,"iid") <- -pred.i$survival.average.iid
        }
        return(risk.i)
      })
      
        ## influence function for the average treatment effect
        ## IF had dimension n.predictions (row), n.times (columns), n.dataTrain (length)
        iid.treatment <- array(NA, dim = c(n.contrasts, n.times, n.obs))
        sdIF.treatment <- matrix(NA, nrow = n.contrasts, ncol = n.times)
        for(iTreat in 1:n.contrasts){ # iTreat <- 1
            term1 <- t(attr(IFrisk[[iTreat]],"iid"))
            ## note: here IFrisk[[iTreat]] is the risk (the influence function is in the attribute "iid")
            term2 <- rowCenter_cpp(IFrisk[[iTreat]], center = pointEstimate$meanRisk[Treatment==contrasts[iTreat],meanRisk])
        
        # we get n * IF instead of IF for the absolute risk. This is why the second term need to be rescaled
        iid.treatment[iTreat,,] <- term1 + t(term2)/n.obs
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
                index.contrasts <- sdIF.fct[, .I[Treatment.A==contrasts[[iCon]] & Treatment.B==contrasts[[iCon2]]]]
                if(se){
                    sdIF.fct[index.contrasts,diff.se := sdIF_diff.contrasts[iiCon,]]  
                    sdIF.fct[index.contrasts,ratio.se := sdIF_ratio.contrasts[iiCon,]]  
                }                   
            }
        }            
                                        # }}}
                                        # {{{ confidence bands
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
      
      # {{{ compute confidence intervals and bands
      crisks <- sdIF.fct[,.(Treatment.A,Treatment.B,time)]
      mrisks <- data.table::data.table(Treatment = pointEstimate$meanRisk$Treatment,
                                       time = times)
      
      mrisks[, meanRisk := pointEstimate$meanRisk$meanRisk]
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
      
      mrisks[, meanRisk := NULL]
      
      # }}}
      bootseeds <- NULL
      ## merge with pointEstimate
      mrisks <- merge(pointEstimate$meanRisk,mrisks,by=key1)
      crisks <- merge(pointEstimate$riskComparison,crisks,by=key2)            
    }
  } else{
    mrisks <- pointEstimate$meanRisk
    crisks <- pointEstimate$riskComparison
    bootseeds <- NULL
    # }}}
  }
  out <- list(meanRisk=mrisks,
              riskComparison=crisks,
              treatment=treatment,
              contrasts=contrasts,
              times=times,
              se = se,
              n.bootstrap=B,
              band = band,
              nsim.band = nsim.band,
              seeds=bootseeds,
              conf.level=conf.level)
  
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
