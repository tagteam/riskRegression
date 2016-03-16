#' @export
predictHazardRR <- function (x, ...) {
  UseMethod("predictHazardRR", x)
}

#' @export
predictSurvProbRR <- function (x, ...) {
  UseMethod("predictSurvProbRR", x)
}

#' @title Predicting hazard or cumulative hazard
#' 
#' @aliases predictSurvProbRR predictSurvProbRR.coxph predictSurvProbRR.cph
#
#' @param object The fitted coxph model
#' @param newdata A data frame containing the values of the variables in the right hand side of 'coxph' for each patient.
#' @param times Vector of times at which to return the estimated hazards/survival
#' @param type should the hazard or the cumulative hazard be returned
#' @param method.baseHaz The implementation to be used for computing the baseline hazard: "dt" or "cpp"
#' @param col.strata Should an additional column containing the strata level be returned
#' 
#' @details 
#' Not suited for time varying cox models
#' 
#' @return A data table (or a list of) containing the predictions for each patient (in rows) and each time (in columns) and the strata (if requested). 
#' 
#' @examples 
#' 
#' library(survival)
#' 
#' set.seed(10)
#' d <- SimSurv(1e3)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d)
#' # table(duplicated(d$time))
#' seq_times <- sample(x = unique(sort(d$time)), size = 100) 
#'
#' res1 <- predictSurvProbRR(fit, newdata = d, times = 10)
#' res2 <- predictSurvProbRR(fit, newdata = d, times = seq_times)
#' res3 <- predictSurvProbRR(fit, newdata = d, times = 10, type = "cumHazard")
#' res4 <- predictSurvProbRR(fit, newdata = d, times = seq_times, type = "cumHazard")
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- predictSurvProbRR(fitS, newdata = d, times = 10)
#' res2S <- predictSurvProbRR(fitS, newdata = d, times = seq_times)
#' res3S <- predictSurvProbRR(fitS, newdata = d, times = 10, type = "cumHazard")
#' res4S <- predictSurvProbRR(fitS, newdata = d, times = seq_times, type = "cumHazard")
#' @export
predictSurvProbRR.coxph <- function (object, newdata, times, type = "Survival",
                                     method.baseHaz = "cpp", col.strata = FALSE) {
  
  times <- sort(times)
  if(any(duplicated(times))){
    stop("predictHazardRR.coxph: argument \"times\" must not contain duplicates \n",
         "sum(duplicated(times)): ",sum(duplicated(times)),"\n")
  }
  
  strataspecials <- attr(object$terms,"specials")$strata 
  test.Nstrata <- is.null(strataspecials)
  
  ## linear predictor
  BaseVar <- attr(object$terms, "term.labels")
  
  if(test.Nstrata){ 
    
    stratalabels <- NULL
    strataVar <- NULL
    
  } else {
    
    stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
    BaseVar <- BaseVar[ BaseVar %in% stratalabels == FALSE]
    
    names.newdata <- names(newdata)
    strataVar <- names.newdata[which(paste0("strata(", names.newdata,")") %in% stratalabels)]
    
  }
  
  ## main
  return(predictSurvProbRR_internal(Lambda0 = baseHazRR(object, method = method.baseHaz, centered = TRUE, lasttime =  times[length(times)], addFirst = TRUE, addLast = TRUE),
                                    newdata = newdata, times = times, type = type, 
                                    Xb = if(length(BaseVar) == 0){ rep(0, nrow(newdata))  
                                      }else if(test.Nstrata){ survival:::predict.coxph(object, newdata, type = "lp") 
                                      }else { rowSums(survival:::predict.coxph(object, newdata = newdata, type = "terms")) },
                                    test.Nstrata = test.Nstrata, stratalabels = stratalabels, strataVar = strataVar, 
                                    col.strata = col.strata)
  )
  
  
}

predictSurvProbRR.cph <- function (object, newdata, times, type = "Survival",
                                   method.baseHaz = "cpp", col.strata = FALSE) {
  
  times <- sort(times)
  if(any(duplicated(times))){
    stop("predictHazardRR.coxph: argument \"times\" must not contain duplicates \n",
         "sum(duplicated(times)): ",sum(duplicated(times)),"\n")
  }

  strataspecials <- attr(object$terms,"specials")$strat
  test.Nstrata <- is.null(strataspecials)
 
  ## linear predictor
  BaseVar <- attr(object$terms, "term.labels")
  
  if(test.Nstrata){ 
    
    stratalabels <- NULL
    strataVar <- NULL
      
  } else {
    
    stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
    BaseVar <- BaseVar[ BaseVar %in% stratalabels == FALSE]
    names.newdata <- names(newdata)
    strataVar <- names.newdata[which(paste0("strat(", names.newdata,")") %in% stratalabels)]
    
  }
  
  ## main
  return(predictSurvProbRR_internal(Lambda0 = baseHazRR(object, method = method.baseHaz, centered = TRUE, lasttime =  times[length(times)], addFirst = TRUE, addLast = TRUE),
                                    newdata = newdata, times = times, type = type, 
                                    Xb = if(length(BaseVar) > 0){predict(object, newdata, type = "lp")}else{ rep(0, nrow(newdata)) },
                                    test.Nstrata = test.Nstrata, stratalabels = stratalabels, strataVar = strataVar, 
                                    col.strata = col.strata)
  )
  
  
}

## Warning times must be in ascending order
predictSurvProbRR_internal <- function (Lambda0, newdata, times, type, Xb,
                                    test.Nstrata, stratalabels, strataVar, col.strata) {
  
  ## preparation
  match.arg(type, choices = c("Survival", "cumHazard", "hazard"), several.ok = TRUE)
  n.type <- length(type)

  resPred <- list()
  
  if(test.Nstrata){ ## no strata
    
    for(iterType in 1:n.type){
      
      if(type[iterType] == "Survival"){ # step function interpolation
        time.index <- prodlim::sindex(jump.times = Lambda0$time, eval.times = times)
        resPred[[iterType]] <- data.table::data.table( Lambda0[time.index, exp(- exp(Xb) %o% .SD[[1]]), .SDcols = "cumHazard"] )
      }else if(type[iterType] == "cumHazard"){ # step function interpolation
        time.index <- prodlim::sindex(jump.times = Lambda0$time, eval.times = times)
        resPred[[iterType]] <- data.table::data.table( Lambda0[time.index, exp(Xb) %o% .SD[[1]], .SDcols = "cumHazard"] )
      }else if(type[iterType] == "hazard"){ # dirac function interpolation
        resPred[[iterType]] <- data.table::data.table(matrix(0, nrow = nrow(newdata), ncol = length(times)))
        time.index <- which(Lambda0$time %in% times)
        if(length(time.index)>0){
          resPred[[iterType]][, which(times %in% Lambda0$time) := data.frame(Lambda0[time.index, exp(Xb) %o% .SD[[1]], .SDcols = "hazard"]), with = FALSE ]
        }
      }
      
      data.table::setnames(resPred[[iterType]], paste0("t",times))
      
      if(col.strata == TRUE){
        resPred[[iterType]][, strata := "1"]
      }
    }
    
  }else{ ## strata
    

      
    nStrata <- length(strataVar)
    #     XXstrata <- interaction(lapply(1:length(strataVar),function(x){paste0(strataVar[x], "=", newdata[[strataVar[x]]])}), drop = TRUE, sep = ", ")
      
    if("XXstrata" %in% names(newdata)){
          stop("predictSurvProbRR.coxph: \'newdata\' must not contains a column named \"XXstrata\" \n")
    }
    
    ##### define the strata # [COULD BE OPTIMIZED in C]
    if(!data.table::is.data.table(newdata)){
      newdata <- data.table::as.data.table(newdata)
    }
    sapply(1:nStrata, function(x){
      newdata[, strataVar[x] := paste0(strataVar[x],"=",.SD[[1]]),.SDcols = strataVar[x], with = FALSE]
    })
    newdata[, XXstrata := interaction(.SD, drop = TRUE, sep = ", ") ,.SDcols = strataVar]
    
    ## compute hazard
    levelsStrata <- levels(newdata$XXstrata)
#      levelsStrata <- levels(XXstrata)
     
    for(iterType in 1:n.type){
      
      resPred[[iterType]] <- data.table::data.table()
      
      for(iterS in levelsStrata){
        
#          indexNew <- which(XXstrata == iterS)
         indexNew <- newdata[,.I[XXstrata == iterS]]
        indexHaz <- Lambda0[,.I[strata == iterS]]
        if(length(indexHaz) == 0){
          stop("predictSurvProbRR_internal: unknown strata ",iterS," \n",
               "existing strata: \"",paste(levels(Lambda0$strata), collapse = "\" \""),"\"\n")
        }
        
        if(type[iterType] == "Survival"){ # step function interpolation
          time.index <- prodlim::sindex(jump.times = Lambda0[indexHaz,time], eval.times=times)
          res.tempo <- Lambda0[indexHaz[time.index], exp(- exp(Xb[indexNew]) %o% .SD[[1]]), .SDcols = "cumHazard"]
        }else if(type[iterType] == "cumHazard"){ # step function interpolation
          time.index <- prodlim::sindex(jump.times = Lambda0[indexHaz,time], eval.times=times)
          res.tempo <- Lambda0[indexHaz[time.index], exp(Xb[indexNew]) %o% .SD[[1]], .SDcols = "cumHazard"]
        }else if(type[iterType] == "hazard"){ # diract function interpolation
          res.tempo <- matrix(0, nrow = length(indexNew), ncol = length(times))
          time.index <- which(Lambda0[indexHaz,time] %in% times)
          if(length(time.index)>0){
            res.tempo[,which(times %in% Lambda0[indexHaz,time])] <- Lambda0[indexHaz[time.index], exp(Xb[indexNew]) %o% .SD[[1]], .SDcols ="hazard"]
          }
        }
        
        resPred[[iterType]] <-  data.table::rbindlist(list(resPred[[iterType]], 
                                                           data.table::data.table(res.tempo, strata = iterS, index = indexNew))
        )
        
      }
      
      data.table::setkey(resPred[[iterType]],index)
      resPred[[iterType]][,index := NULL]
      data.table::setnames(resPred[[iterType]], c(paste0("t",times),"strata"))
      if(col.strata == FALSE){
        resPred[[iterType]][,strata := NULL]
      }
    }
    
  }
  
  ## export
  names(resPred) <- type
  switch(as.character(n.type),
         "1" = return(resPred[[1]]),
         resPred
  )
  
}
  
  


