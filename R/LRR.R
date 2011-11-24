LRR <- function(formula,data,times,cause,cens.model,...){
  riskRegression(formula=formula,data=data,times=times,link="logistic",cause=cause,cens.model=cens.model,...)
}
