getPseudoValues <- function(formula,
                            data,
                            times,
                            responseType,
                            censModel="marginal",
                            cause,
                            order=1){
  responseType <- match.arg(responseType,c("survival","competing.risks"))
  censModel <- match.arg(tolower(censModel),tolower(c("marginal","cox","nonpar","aalen")))
  obj <- 1
  class(obj) <- paste(censModel,responseType,"order",order,sep=".")
  UseMethod("getPseudoValues",obj)
}

## getPseudoValues.marginal.competing.risk.order.1 <- function(formula=formula,
                                                            ## data=data,
                                                            ## times=times,
                                                            ## responseType=responseType,
                                                            ## censModel="marginal",
                                                            ## event=event,
                                                            ## order=1){
  ## print("FIXME: pseudo values for competing risks of order 1")
## }
## getPseudoValues.marginal.competing.risk.order.2 <- function(formula=formula,
                                                            ## data=data,
                                                            ## times=times,
                                                            ## responseType=responseType,
                                                            ## censModel="marginal",
                                                            ## event=event,
                                                            ## order=1){
  ## print("FIXME: pseudo values for competing risks of order 2")
## }
## getPseudoValues.marginal.survival.order.1 <- function(formula=formula,
                                                      ## data=data,
                                                      ## times=times,
                                                      ## responseType=responseType,
                                                      ## censModel="marginal",
                                                      ## event=event,
                                                      ## order=1){
  ## print("Pseudo values for survival of order 1 based on independent censoring")
  ## fit <- prodlim(formula,data)
  ## jackknife(fit,times)
## }
## getPseudoValues.marginal.survival.order.2 <- function(formula=formula,
                                                      ## data=data,
                                                      ## times=times,
                                                      ## responseType=responseType,
                                                      ## censModel="marginal",
                                                      ## event=event,
                                                      ## order=1){
  ## print("FIXME: pseudo values for survival of order 2")
  ## fit <- prodlim(formula,data)
  ## pseudoOrder(fit,times)
## }
