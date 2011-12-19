CSC <- function (formula,data,cause,survtype="hazard",...){
  # {{{ type
  survtype <- match.arg(survtype,c("hazard","survival"))
  # }}}
  # {{{ formulae
  if (class(formula)=="formula") formula <- list(formula)
  else if (!length(formula)<=4) stop("Formula can either be a single formula (to be used for both cause specific hazards) or a list of formulae, one for each cause.")
  responseFormula <- reformulate("1", formula[[1]][[2]])
  # }}}
  # {{{ response
  require(survival)
  call <- match.call()
  # get information from formula
  mf <- model.frame(responseFormula, data = data, na.action = na.omit)
  response <- model.response(mf)
  time <- response[, "time"]
  status <- response[, "status"]
  event <- getEvent(response)
  # }}}
  # {{{ event times
  eventTimes <- unique(sort(as.numeric(time[status != 0])))
  # }}}
  # {{{ causes
  causes <- getStates(response)
  if (survtype=="hazard")
    NC <- length(causes)
  else
    NC <- 2 # cause of interest and overall survival
  if (length(formula)!=NC && length(formula)>1) stop("Wrong number of formulae. Should be ",NC,".")
  if (length(formula)==1) {
    ## warning("The same formula used for all causes")
    formula <- lapply(1:NC,function(x)formula[[1]])
  }
  # }}}
  # {{{ find the cause of interest
  if (missing(cause)){
    theCause <- causes[1]
  }
  else{
    if ((foundCause <- match(as.character(cause),causes,nomatch=0))==0)
      stop(paste("Requested cause: ",cause," Available causes: ", causes))
    else
      theCause <- foundCause
  }
  otherCauses <- causes[-match(theCause,causes)]
  # }}}
  # {{{ fit Cox models
  if (survtype=="hazard"){
    CoxModels <- lapply(1:NC,function(x){
      if (x==1)
        causeX <- theCause
      else
        causeX <- otherCauses[x-1]
      formulaX <- formula[[x]]
      covData <- model.frame(formulaX,data=data)
      response <- model.response(covData)
      time <- response[, "time"]
      event <- getEvent(response)
      statusX <- event==causeX
      workData <- cbind(time=time,status=statusX,covData[,-1,drop=FALSE])
      formulaXX <- as.formula(paste("Surv(time,status)",as.character(delete.response(terms.formula(formulaX)))[[2]],sep="~"))
      fit <- coxph(formulaXX, data = workData,...)
      fit$call$formula <- fit$formula
      fit
    })
    names(CoxModels) <- paste("Cause",causes)
  }
  else{
    CoxModels <- lapply(1:2,function(x){
      causeX <- theCause
      formulaX <- formula[[x]]
      covData <- model.frame(formulaX,data=data)
      response <- model.response(covData)
      time <- response[, "time"]
      event <- getEvent(response)
      if (x==1){
        statusX <- event==causeX
      }
      else{
        statusX <- response[,"status"]
      }
      workData <- cbind(time=time,status=statusX,covData[,-1,drop=FALSE])
      formulaXX <- as.formula(paste("Surv(time,status)",as.character(delete.response(terms.formula(formulaX)))[[2]],sep="~"))
      coxph(formulaXX, data = workData,...)
    })
    names(CoxModels) <- c(paste("Cause",theCause),"OverallSurvival")
  }
  # }}}
  out <- list(call=call,
              models=CoxModels,
              response=response,
              eventTimes=eventTimes,
              survtype=survtype,
              theCause=theCause,
              causes=c(theCause,otherCauses))
  class(out) <- "CauseSpecificCox"
  out
}


