# -----------------------------cumIncCox-----------------------------
# fits a cox regression for the cause specific hazards for cause of
# interest and a cox regression for total hazard (to the overall
# survival)

cumIncCox <- function (formula,data,cause){
  require(survival)
  call <- match.call()
  # get information from formula
  mf <- model.frame(formula, data = data, na.action = na.omit)
  response <- model.response(mf)
  time <- response[, "time"]
  status <- response[, "status"]
  event <- getEvent(response)
  # get covariates
  covFormula <- delete.response(terms.formula(formula))
  # find the requested cause of interest
  if (missing(cause)){
    cause <- getStates(response)[1]
  }
  if ((foundCause <- match(as.character(cause),getStates(response),nomatch=0))==0)
    stop(paste("Requested cause: ",cause," Available causes: ", getStates(response)))
  else
    cause <- foundCause
  NC <- length(unique(event[status != 0]))
  eventTimes <- unique(sort(as.numeric(time[status != 0])))
  # workData with the cause of interest in dummy code 
  workData <- cbind(time=time,status=status,theCause=as.numeric(event==cause),model.frame(covFormula,data))
  # fit a formula with the cause of interest
  SurvForm <- update(Surv(time, theCause) ~ ., covFormula)
  # fit cox regression model for cause of interest
  csHazCox <- coxph(SurvForm, data = workData)
  # fit formula for overall survival 
  Form <- update(Surv(time, status) ~ ., covFormula)
  # fit cox regression for the overall survival
  survCox <- coxph(Form, data = workData)
  # out
  out <- list(eventTimes = eventTimes,call = call, causes= NC, csHazCox = csHazCox, survCox = survCox)
  class(out) <- "cumIncCox"
  invisible(out)
}







