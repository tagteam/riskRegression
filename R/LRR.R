LRR <- function(formula,data,times,cause,cens.model,...){
  fit <- riskRegression(formula=formula,data=data,times=times,link="logistic",cause=cause,cens.model=cens.model,...)
  fit$call <- match.call()
  fit
}
