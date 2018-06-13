
                                        # {{{ header
## * predictCox (documentation)
#' @title Fast computation of survival probabilities, hazards and cumulative hazards from Cox regression models 
#' @name predictCox
#' 
#' @description Fast routine to get baseline hazards and subject specific hazards
#' as well as survival probabilities from a \code{survival::coxph} or \code{rms::cph} object
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata [data.frame or data.table]  Contain the values of the predictor variables
#' defining subject specific predictions.
#' Should have the same structure as the data set used to fit the \code{object}.
#' @param times [numeric vector] Time points at which to return
#' the estimated hazard/cumulative hazard/survival.
#' @param centered [logical] If \code{TRUE} return prediction at the
#'     mean values of the covariates \code{fit$mean}, if \code{FALSE}
#'     return a prediction for all covariates equal to zero.  in the
#'     linear predictor. Will be ignored if argument \code{newdata} is
#'     used. For internal use.
#' @param type [character vector] the type of predicted value. Choices are \itemize{
#'     \item \code{"hazard"} the baseline hazard function when
#'     argument \code{newdata} is not used and the hazard function
#'     when argument \code{newdata} is used.  \item \code{"cumhazard"}
#'     the cumulative baseline hazard function when argument
#'     \code{newdata} is not used and the cumulative hazard function
#'     when argument \code{newdata} is used.  \item \code{"survival"}
#'     the survival baseline hazard function when argument
#'     \code{newdata} is not used and the cumulative hazard function
#'     when argument \code{newdata} is used.  } Several choices can be
#'     combined in a vector of strings that match (no matter the case)
#'     strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param keep.strata [logical] If \code{TRUE} add the (newdata) strata
#'     to the output. Only if there any.
#' @param keep.times [logical] If \code{TRUE} add the evaluation times
#'     to the output.
#' @param keep.newdata [logical] If \code{TRUE} add the value of the
#'     covariates used to make the prediction in the output list.
#' @param se [logical] If \code{TRUE} compute and add the standard errors to the output.
#' @param band [logical] If \code{TRUE} compute and add the quantiles for the confidence bands to the output.
#' @param iid [logical] If \code{TRUE} compute and add the influence function to the output.
#' @param confint [logical] If \code{TRUE} compute and add the confidence intervals/bands to the output.
#' They are computed applying the \code{confint} function to the output.
#' @param average.iid [logical] If \code{TRUE} add the average of the influence function over \code{newdata} to the output.
#' @param store.iid [character] Implementation used to estimate the influence function and the standard error.
#' Can be \code{"full"} or \code{"minimal"}.
#' @param ... not used.
#' 
#' @details
#' When the argument \code{newdata} is not specified, the function computes the baseline hazard estimate.
#' See (Ozenne et al., 2017) section "Handling of tied event times".
#'
#' Otherwise the function computes survival probabilities with confidence intervals/bands.
#' See (Ozenne et al., 2017) section "Confidence intervals and confidence bands for survival probabilities".
#' The survival is computed using the exponential approximation (equation 3).
#'
#' A detailed explanation about the meaning of the argument \code{store.iid} can be found
#' in (Ozenne et al., 2017) Appendix B "Saving the influence functions".
#' 
#' The function is not compatible with time varying predictor variables.
#' 
#' The centered argument enables us to reproduce the results obtained with the \code{basehaz}
#' function from the survival package but should not be modified by the user.
#'
#' The iid decomposition is output using an array containing the value of the influence
#' of each subject used to fit the object (dim 3),
#' for each subject in newdata (dim 1),
#' and each time (dim 2).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#'
#' @references
#' Brice Ozenne, Anne Lyngholm Sorensen, Thomas Scheike, Christian Torp-Pedersen and Thomas Alexander Gerds.
#' riskRegression: Predicting the Risk of an Event using Cox Regression Models.
#' The R Journal (2017) 9:2, pages 440-460.
#' 
#' @seealso
#' \code{\link{confint.predictCox}} to compute confidence intervals/bands.
#' \code{\link{autoplot.predictCox}} to display the predictions.

## * predictCox (examples)
#' @rdname predictCox
#' @examples 
#' library(survival)
#'
#' #### generate data ####
#' set.seed(10)
#' d <- sampleData(40,outcome="survival") ## training dataset
#' nd <- sampleData(4,outcome="survival") ## validation dataset
#' d$time <- round(d$time,1) ## create tied events
#' # table(duplicated(d$time))
#' 
#' #### stratified Cox model ####
#' fit <- coxph(Surv(time,event)~X1 + strata(X2) + X6,
#'              data=d, ties="breslow", x = TRUE, y = TRUE)
#' 
#' ## compute the baseline cumulative hazard
#' fit.haz <- predictCox(fit)
#' cbind(survival::basehaz(fit), fit.haz$cumhazard)
#'
#' ## compute individual specific cumulative hazard and survival probabilities 
#' fit.pred <- predictCox(fit, newdata=nd, times=c(3,8), se = TRUE, band = TRUE)
#' fit.pred
#'
#' ####  other examples ####
#' # one strata variable
#' fitS <- coxph(Surv(time,event)~strata(X1)+X2,
#'               data=d, ties="breslow", x = TRUE, y = TRUE)
#' 
#' predictCox(fitS)
#' predictCox(fitS, newdata=nd, times = 1)
#'
#' # two strata variables
#' set.seed(1)
#' d$U=sample(letters[1:5],replace=TRUE,size=NROW(d))
#' d$V=sample(letters[4:10],replace=TRUE,size=NROW(d))
#' nd$U=sample(letters[1:5],replace=TRUE,size=NROW(nd))
#' nd$V=sample(letters[4:10],replace=TRUE,size=NROW(nd))
#' fit2S <- coxph(Surv(time,event)~X1+strata(U)+strata(V)+X2,
#'               data=d, ties="breslow", x = TRUE, y = TRUE)
#'
#' cbind(survival::basehaz(fit2S),predictCox(fit2S,type="cumhazard")$cumhazard)
#' predictCox(fit2S)
#' predictCox(fitS, newdata=nd, times = 3)
#'
#' # left truncation
#' test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
#'               stop=c(2,3,6,7,8,9,9,9,14,17), 
#'               event=c(1,1,1,1,1,1,1,0,0,0), 
#'               x=c(1,0,0,1,0,1,1,1,0,0)) 
#' m.cph <- coxph(Surv(start, stop, event) ~ 1, test2, x = TRUE)
#' as.data.table(predictCox(m.cph))
#'
#' basehaz(m.cph)
# }}}

## * predictCox (code)
#' @rdname predictCox
#' @export
predictCox <- function(object,
                       newdata=NULL,
                       times,
                       centered = TRUE,
                       type=c("cumhazard","survival"),
                       keep.strata = TRUE,
                       keep.times = TRUE,
                       keep.newdata = FALSE,
                       se = FALSE,
                       band = FALSE,
                       iid = FALSE,
                       confint = (se+band)>0,
                       average.iid = FALSE,
                       store.iid = "full"){
  status=statusM1=NULL
  
  # {{{ treatment of times and stopping rules
  
#### Extract elements from object ####
                                        # we need:          - the total number of observations, the status and eventtime for each observation
                                        #                   - the strata corresponding to each observation in the training set
                                        #                   - the name of each strata
                                        #                   - the value of the linear predictor for each observation in the training set
  if (se==1L || iid==1L){
      if (missing(newdata)) stop("Argument 'newdata' is missing. Cannot compute standard errors in this case.")
  }
  infoVar <- coxVariableName(object)
  is.strata <- infoVar$is.strata
  if (missing(times)) {
      nTimes <- 0
      times <- numeric(0)
  }else{
      nTimes <- length(times)
  }
  needOrder <- (nTimes>0 && is.unsorted(times))
  if (needOrder) {
      oorder.times <- order(order(times))
      times.sorted <- sort(times)
  }else{
      if (nTimes==0)
          times.sorted <- numeric(0)
      else
          times.sorted <- times
  }

  object.n <- coxN(object)
  object.design <- coxDesign(object)
  object.status <- object.design[["status"]]
  object.start <- object.design[["start"]]
  object.stop <- object.design[["stop"]]
  object.strata <- coxStrata(object,
                             data = NULL,
                             strata.vars = infoVar$strata.vars, 
                             strata.levels = infoVar$strata.levels)
  object.levelStrata <- levels(object.strata)
                                        # if we predict the hazard for newdata then there is no need to center the covariates
  object.eXb <- exp(coxLP(object, data = NULL, center = if(is.null(newdata)){centered}else{FALSE}))

  object.baseEstimator <- coxBaseEstimator(object) 
  nVar <- length(infoVar$lpvars)
  
  #### checks ####
  if(object.baseEstimator == "exact"){
      stop("Prediction with exact handling of ties is not implemented.\n")
  }
  if(nTimes>0 && any(is.na(times))){
      stop("Missing (NA) values in argument \'times\' are not allowed.\n")
  }
  type <- tolower(type)
  if(!is.null(object$weights) && !all(object$weights==1)){
      stop("predictCox does not know how to handle Cox models fitted with weights \n")
  }
  if(any(type %in% c("hazard","cumhazard","survival") == FALSE)){
      stop("type can only be \"hazard\", \"cumhazard\" or/and \"survival\" \n") 
  }
  if(any(object.design[,"start"]!=0)){
      warning("The current version of predictCox was not designed to handle left censoring \n",
              "The function may be used on own risks \n") 
  }    
  if(!is.null(object$naive.var)){
      stop("predictCox does not know how to handle fraitly \n") 
  }
  if(any(is.na(coef(object)))){
      stop("Incorrect object",
           "One or several model parameters have been estimated to be NA \n")
  }
                                        # }}}
                                        # {{{ computation of the baseline hazard
  if(!is.null(newdata)){
      name.regressor <- c(infoVar$lpvars.original, infoVar$strata.vars.original)
      if(length(name.regressor) > 0 && any(name.regressor %in% names(newdata) == FALSE)){
          stop("Missing variables in argument \'newdata\': \"",
               paste0(setdiff(name.regressor,names(newdata)), collapse = "\" \""),
               "\"\n")
      }
    
      new.n <- NROW(newdata)
      newdata <- as.data.table(newdata)
      new.eXb <- exp(coxLP(object, data = newdata, center = FALSE))
      
      new.strata <- coxStrata(object, data = newdata, 
                              sterms = infoVar$sterms, 
                              strata.vars = infoVar$strata.vars, 
                              levels = object.levelStrata, 
                              strata.levels = infoVar$strata.levels)
    
    new.levelStrata <- levels(new.strata)
  }
  
  #### baseline hazard ####
  nStrata <- length(object.levelStrata)
  if(is.strata){etimes.max <- tapply(object.stop, object.strata, max) }else{ etimes.max <- max(object.stop) } # last event time
  
  ## sort the data
  dt.prepare <- data.table::data.table(start = object.start,
                                       stop = object.stop,
                                       status = object.status,
                                       eXb = object.eXb,
                                       strata = as.numeric(object.strata) - 1)
  dt.prepare[, statusM1 := 1-status] # sort by statusM1 such that deaths appear first and then censored events
  data.table::setkeyv(dt.prepare, c("strata","stop","start","statusM1"))
  
  ## compute the baseline hazard
  Lambda0 <- baseHaz_cpp(starttimes = dt.prepare$start,
                         stoptimes = dt.prepare$stop,
                         status = dt.prepare$status,
                         eXb = dt.prepare$eXb,
                         strata = dt.prepare$strata,
                         nPatients = object.n,
                         nStrata = nStrata,
                         emaxtimes = etimes.max,
                         predtimes = times.sorted,
                         cause = 1,
                         Efron = (object.baseEstimator == "efron"))
                                        # }}}
  
#### compute hazard and survival ####        
  if (is.null(newdata)){  
                                        # {{{ results from the training dataset
      if (!("hazard" %in% type)){ Lambda0$hazard <- NULL } 
      if ("survival" %in% type){  # must be before cumhazard
          Lambda0$survival = exp(-Lambda0$cumhazard)
      }
      if (!("cumhazard" %in% type)){ Lambda0$cumhazard <- NULL } 
      if (keep.times==FALSE){
          Lambda0$times <- NULL
      } 
      if (is.strata == TRUE && keep.strata==1L){ ## rename the strata value with the correct levels
          Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
      }else{
          Lambda0$strata <- NULL
      }
      Lambda0$lastEventTime <- etimes.max
      Lambda0$se <- FALSE
      Lambda0$band <- FALSE
      Lambda0$type <- type

      ## Lambda0$time instead of Lambda0$times for compatibility with survival::basehaz
      class(Lambda0) <- "predictCox"
      return(Lambda0)
                                        # }}}
  } else {
    
    # {{{ predictions in new dataset
    out <- list()
    if(missing(times) || nTimes==0){
      stop("Time points at which to evaluate the predictions are missing \n")
    }
    
      ## subject specific hazard
    if (is.strata==FALSE){
      if ("hazard" %in% type){
        out$hazard <- (new.eXb %o% Lambda0$hazard)
        if (needOrder) out$hazard <- out$hazard[,oorder.times,drop=0L]
      }
      if ("cumhazard" %in% type || "survival" %in% type){
        cumhazard <- new.eXb %o% Lambda0$cumhazard
        if ("cumhazard" %in% type){
          if (needOrder)
            out$cumhazard <- cumhazard[,oorder.times,drop=0L]
          else
            out$cumhazard <- cumhazard
        }
        if ("survival" %in% type){
          out$survival <- exp(-cumhazard)
          if (needOrder)
            out$survival <- out$survival[,oorder.times,drop=0L]
        }
      }
      
    }else{ 
      
      ## initialization
      if ("hazard" %in% type){
        out$hazard <- matrix(0, nrow = new.n, ncol = nTimes)
      }
      if ("cumhazard" %in% type){
        out$cumhazard <- matrix(NA, nrow = new.n, ncol = nTimes)                
      }
      if ("survival" %in% type){
        out$survival <- matrix(NA, nrow = new.n, ncol = nTimes)               
      }
      if (is.strata == TRUE){ ## rename the strata value with the correct levels
        Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = object.levelStrata)
      }

      ## loop across strata
      for(S in new.levelStrata){
        id.S <- Lambda0$strata==S
        newid.S <- new.strata==S

        if ("hazard" %in% type){
          out$hazard[newid.S,] <- new.eXb[newid.S] %o% Lambda0$hazard[id.S]
          if (needOrder)
            out$hazard[newid.S,] <- out$hazard[newid.S,oorder.times,drop=0L]
        }
        if ("cumhazard" %in% type || "survival" %in% type){
          cumhazard.S <-  new.eXb[newid.S] %o% Lambda0$cumhazard[id.S]
          if ("cumhazard" %in% type){
            if (needOrder){
              out$cumhazard[newid.S,] <- cumhazard.S[,oorder.times,drop=0L]
            } else{
              out$cumhazard[newid.S,] <- cumhazard.S
            }
          }
          if ("survival" %in% type){
            if (needOrder){
              out$survival[newid.S,] <- exp(-cumhazard.S)[,oorder.times,drop=0L]
            }else{
              out$survival[newid.S,] <- exp(-cumhazard.S)
            }
          }
        }
      }
    }
    # }}}
    # {{{ standard error
    
    if(se || band || iid || average.iid){
      if(se && "hazard" %in% type){
        stop("confidence intervals cannot be computed for the hazard \n")
      }
      if(band && "hazard" %in% type){
        stop("confidence bands cannot be computed for the hazard \n")
      }
      
      if(nVar > 0){
        # remove response variable
        f.object <- stats::reformulate(attr(stats::terms(coxFormula(object)),"term.label"),
                                       response = NULL)
        # use prodlim to get the design matrix
        terms.newdata <- stats::terms(f.object, special = coxSpecialStrata(object), data = newdata)
        new.LPdata <- prodlim::model.design(stats::terms(terms.newdata),
                                            data = newdata,
                                            specialsFactor = TRUE,
                                            dropIntercept = TRUE)$design
        if(NROW(new.LPdata)!=NROW(newdata)){
          stop("NROW of the design matrix and newdata differ \n",
               "maybe because newdata contains NA values \n")
        }
      }else{
        new.LPdata <- matrix(0, ncol = 1, nrow = new.n)
      }
      
      ## Computation of the influence function and/or the standard error
      outSE <- calcSeCox(object, times = times.sorted, nTimes = nTimes, type = type,
                         Lambda0 = Lambda0, object.n = object.n, object.time = object.stop, object.eXb = object.eXb, object.strata = object.strata, nStrata = nStrata,
                         new.eXb = new.eXb, new.LPdata = new.LPdata, new.strata = new.strata,
                         new.survival = out$survival,
                         nVar = nVar, 
                         export = c("iid"[(iid+band)>0],"se"[(se+band)>0],"average.iid"[average.iid==TRUE]), store.iid = store.iid)

      ## restaure orginal time ordering
      if((iid+band)>0){
        if ("hazard" %in% type){
          if (needOrder)
            out$hazard.iid <- outSE$hazard.iid[,oorder.times,,drop=0L]
          else
            out$hazard.iid <- outSE$hazard.iid
        }
        if ("cumhazard" %in% type){
          if (needOrder)
            out$cumhazard.iid <- outSE$cumhazard.iid[,oorder.times,,drop=0L]
          else
            out$cumhazard.iid <- outSE$cumhazard.iid
        }
        if ("survival" %in% type){
          if (needOrder)
            out$survival.iid <- outSE$survival.iid[,oorder.times,,drop=0L]
          else
            out$survival.iid <- outSE$survival.iid
        }
      }
      if(average.iid == TRUE){
        if ("cumhazard" %in% type){
          if (needOrder)
            out$cumhazard.average.iid <- outSE$cumhazard.average.iid[,oorder.times,drop=0L]
          else
            out$cumhazard.average.iid <- outSE$cumhazard.average.iid
        }
        if ("survival" %in% type){
          if (needOrder)
            out$survival.average.iid <- outSE$survival.average.iid[,oorder.times,drop=0L]
          else
            out$survival.average.iid <- outSE$survival.average.iid
        }
      }
      if((se+band)>0){
        if ("cumhazard" %in% type){
          if (needOrder){
            out$cumhazard.se <- outSE$cumhazard.se[,oorder.times,drop=0L]
          }else{
            out$cumhazard.se <- outSE$cumhazard.se
          }          
        }
        if ("survival" %in% type){
            if (needOrder){
                out$survival.se <- outSE$survival.se[,oorder.times,drop=0L]
            } else{
                out$survival.se <- outSE$survival.se
            }          
        }
      }      
    }    
    # }}}
                                        # {{{ export 
      if (keep.times==TRUE) out <- c(out,list(times=times))
      if (is.strata && keep.strata==TRUE) out <- c(out,list(strata=new.strata))
    
      out <- c(out,list(lastEventTime = etimes.max,
                        se = se,
                        band = band,
                        type = type))
      if( keep.newdata==TRUE){
          out$newdata <- newdata[, coxCovars(object), with = FALSE]
      }
      class(out) <- "predictCox"


      if(confint){
          out <- stats::confint(out)
      }
      if(band && se==FALSE){
          out[paste0(type,".se")] <- NULL
      }
      if(band && iid==FALSE){
          out[paste0(type,".iid")] <- NULL
      }
      return(out)
                                        # }}}
  }
  
}

