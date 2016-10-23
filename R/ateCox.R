#' Compute the average treatment effect using CSC
#'
#' Use the g-formula to estimate the average treatment effect
#' @param formula formula argument for the model
#' @param model the function to be called to estimate the outcome
#' @param predictor a method to be applied to the model to obtain the predicted outcome for each time and patient.
#' @param data data argument for CSC
#' @param treatment name of the column containing the treatment
#'     variable
#' @param contrasts the levels of treatment variable to be compared
#' @param times time points at which to evaluate risks
#' @param cause the cause of interest
#' @param B the number of bootstrap replications used to
#'     compute the confidence interval.
#' @param seed An integer used to generate seeds for bootstrap and to achieve
#'        reproducibility of the bootstrap confidence intervals. 
#' @param handler parallel handler for bootstrap. Either "mclapply" or
#'     "sfClusterApplyLB"
#' @param mc.cores Passed on to \code{mclapply} or
#'     \code{sfClusterApplyLB}. The number of cores to use, i.e. at
#'     most how many child processes will be run simultaneously.  The
#'     option is initialized from environment variable MC_CORES if
#'     set.
#' @param conf.level the level for bootstrap confidence intervals
#' @param return.model Logical. if \code{TRUE} the fitted CSC is
#'     returned as part of the output.
#' @param verbose Logical. If \code{TRUE} inform about estimated run
#'     time.
#' @param ... passed to CSC
#' @return A list with: point estimates, bootstrap quantile confidence
#'     intervals model: the CSC model (optional)
#' 
#' @details 
#' The predictor must have (at least) object, newdata, times, and cause as arguments.
#' 
#' @examples 
#' library(survival)
#' library(rms)
#' 
#' set.seed(10)
#' n <- 1e3
#' 
#' ## Cox model
#' dtS <- sampleData(n,outcome="survival")
#' dtS$time <- round(dtS$time,1)
#' seqtimes <- sample(x = unique(sort(dtS$time)), size = 100) 
#' dtS$X1 <- factor(rbinom(n, prob = c(0.2,0.3,0.2) , size = 3), labels = paste0("T",0:3))
#' 
#' ateCox(formula = Surv(time,event)~ X1+X2, model = "cph", data = dtS, treatment = "X1", contrasts = NULL,
#'         times = 5:7, B = 3, y = TRUE, mc.cores=1)
#' 
#' ## Cause specific cox model
#' dt <- sampleData(1e3,outcome="competing.risks")
#' dt$time <- round(dt$time,1)
#' seqtimes <- sample(x = unique(sort(dt$time)), size = 100) 
#' dt$X1 <- factor(rbinom(1e3, prob = c(0.2,0.3,0.2) , size = 3), labels = paste0("T",0:3))
#' ateCox(formula = Hist(time,event)~ X1+X2, model = "CSC", data = dt, treatment = "X1", contrasts = NULL,
#'         times = 5:7, cause = 1, B = 3, mc.cores=1)
#'
#' @export
ateCox <- function(formula,
                   model, 
                   predictor,
                   data,
                   treatment,
                   contrasts = NULL,
                   times,
                   cause,
                   B = 0,
                   seed,
                   handler=c("mclapply","foreach"),
                   mc.cores = 1,
                   conf.level = .95,
                   return.model=TRUE,
                   verbose=TRUE,
                   ...){
  meanRisk=Treatment=ratio=Treatment.A=Treatment.B=b=NULL
  #### Prepare
  if(treatment %in% names(data) == FALSE){
    stop("The data set does not seem to have a variable ",treatment," (argument: treatment). \n")
  }
  
  test.CR <- !missing(cause) # test whether the argument cause has been specified, i.e. it is a competing risk model
  if(test.CR==FALSE){cause <- NA}
  
  if(missing(predictor)){
    if(model %in% c("cph","coxph")){
      if(model == "coxph") loadNamespace("survival")
      if(model == "cph") loadNamespace("rms")
      predictor <- function(object, newdata, times, cause){
        riskRegression::predictCox(object, newdata = newdata, times= times, type = "survival", keep.times = FALSE, keep.strata = FALSE)$survival
      }
    }else if(model == "CSC"){
      predictor <- function(object, newdata, times, cause){
        stats::predict(object, newdata = newdata, times= times, cause = cause)
      }
    }else if(model == "LRR"){
      predictor <- pec::predictEventProb
    }else if(model == "ranger"){
      predictor <- function(object, newdata, times, cause){
        pred <- stats::predict(object, data = newdata)
        indexTime <- prodlim::sindex(jump.times = pred$unique.death.times, eval.times = times)
        return(pred$survival[,indexTime,drop=FALSE])
      }
    }else if(model == "rfsrc"){
      predictor <- function(object, newdata, times, cause){
        pred <- stats::predict(object, newdata = newdata, times= times)
        indexTime <- prodlim::sindex(jump.times = pred$time.interest, eval.times = times)
        return(pred$survival[,indexTime,drop=FALSE])
      }
    }else{
      predictor <- "predict"
    }
  }
  data[[treatment]] <- factor(data[[treatment]])
  
  if(is.null(contrasts)){
    levels <- levels(data[[treatment]])
    contrasts <- levels(data[[treatment]])
    if (length(contrasts)>5) stop("Treatment variable has more than 5 levels.\nIf this is not a mistake, use argument `contrasts'.")
  }
  n.contrasts <- length(contrasts)
  n.times <- length(times)
  n.obs <- NROW(data)
  
  #### calc G formula
  Gformula <- function(data, treatment, contrasts, formula, times, cause, model, predictor, return.model, ...){
    
    ## fit
    if(model == "LRR"){
      #   model.fit <- lapply(times, function(x){if(x==0){NA}else{LRR(formula = formula, time = c(x,, args = list(formula = formula, data = data, ...) ) 
      model.fit <- do.call(model, args = list(formula = formula, data = data, times = times, ...) ) 
    }else{
      model.fit <- do.call(model, args = list(formula = formula, data = data, ...) ) 
    }
    
    # }
    
    meanRisk <- lapply(1:n.contrasts,function(i){
      ## prediction for the hypothetical worlds in which every subject is treated with the same treatment
      data.i <- data
      data.i[[treatment]] <- factor(contrasts[i], levels = levels)
      if (length(times)>1){
        colMeans(do.call(predictor, args = list(model.fit, newdata = data.i, times = times, cause = cause)))
      } else {
        mean(do.call(predictor, args = list(model.fit, newdata = data.i, times = times, cause = cause)))
      }
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
    if(return.model){out <- c(out,list(model=model.fit))}
    out
  }
  #### point estimate
  estimateTime <- system.time(pointEstimate <- Gformula(data=data,
                                                        treatment=treatment,
                                                        contrasts=contrasts,
                                                        formula=formula,
                                                        times=times,
                                                        cause=cause,
                                                        model=model,
                                                        predictor=predictor,
                                                        return.model=return.model,
                                                        ...))
  #### Bootstrap
  if(B>0){
    if (verbose==TRUE)
      message(paste0("Approximated bootstrap netto run time (without time for copying data to cores):\n",
                     round(estimateTime["user.self"],2),
                     " seconds times ",
                     B,
                     " bootstraps / ",
                     mc.cores,
                     " cores = ",
                     round(estimateTime["user.self"]*B/mc.cores,2)," seconds.\n",
                     "To reduce computation time you may consider a coarser time grid,\n e.g., round event times to weeks, months or years ...\n"))
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
      addPackage <- if(grep("package:",find(model))){gsub("package:","",find(model))}else{NULL}
      
      boots <- foreach::`%dopar%`(foreach::foreach(b=1:B,.packages=c("riskRegression",addPackage),.export=NULL), {
        set.seed(bootseeds[[b]])
        dataBoot <- data[sample(1:n.obs, size = n.obs, replace = TRUE),]
        tryCatch(Gformula(data=dataBoot, treatment=treatment, contrasts=contrasts, formula=formula, times=times, cause=cause, model=model, predictor=predictor, return.model = FALSE, ...),
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
        tryCatch(Gformula(data=dataBoot, treatment=treatment, contrasts=contrasts, formula=formula, times=times, cause=cause, model=model, predictor=predictor, return.model = FALSE, ...),
                 error = function(x){return(NULL)})
      }, mc.cores = mc.cores)
    }
    
    ## gc()
    meanRisksBoot <- data.table::rbindlist(lapply(boots,function(x)x$meanRisk))
    riskComparisonsBoot <- data.table::rbindlist(lapply(boots,function(x)x$riskComparison))
    alpha <- 1-conf.level
    mrisks <- meanRisksBoot[,data.table::data.table(meanRiskBoot=mean(meanRisk, na.rm = TRUE),
                                                    lower=quantile(meanRisk,alpha/2, na.rm = TRUE),
                                                    upper=quantile(meanRisk,1-(alpha/2), na.rm = TRUE),
                                                    n.boot=sum(!is.na(meanRisk))),
                            keyby=c("Treatment","time")]
    
    crisks <- riskComparisonsBoot[,data.table::data.table(diffMeanBoot=mean(diff, na.rm = TRUE),
                                                          diff.lower=quantile(diff,alpha/2, na.rm = TRUE),
                                                          diff.upper=quantile(diff,1-(alpha/2), na.rm = TRUE),
                                                          ratioMeanBoot=mean(ratio, na.rm = TRUE),
                                                          ratio.lower=quantile(ratio,alpha/2, na.rm = TRUE),
                                                          ratio.upper=quantile(ratio,1-(alpha/2), na.rm = TRUE),
                                                          n.boot=sum(!is.na(diff))),
                                  keyby=c("Treatment.A","Treatment.B","time")]
  }else{
    
    mrisks <- data.table::data.table(Treatment = pointEstimate$meanRisk$Treatment, time = times, meanRiskBoot = NA, lower = NA, upper = NA)
    crisks <- data.table::data.table(Treatment.A = pointEstimate$riskComparison$Treatment.A, Treatment.B = pointEstimate$riskComparison$Treatment.B,  time = times,
                                     diffMeanBoot = NA, diff.lower = NA, diff.upper = NA, ratioMeanBoot = NA, ratio.lower = NA, ratio.upper = NA)
    bootseeds <- NULL
  }
  
  ## merge bootstrap with pointEstimate
  mrisks <- merge(pointEstimate$meanRisk,mrisks,by=c("Treatment","time"))
  crisks <- merge(pointEstimate$riskComparison,crisks,by=c("Treatment.A","Treatment.B","time"))
  out <- list(meanRisk=mrisks,
              riskComparison=crisks,
              treatment=treatment,
              contrasts=contrasts,
              times=times,
              n.bootstrap=B,
              seeds=bootseeds)
  
  class(out) <- "ateCox"
  out
}
