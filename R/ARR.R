ARR <- function(formula,data,times,cause,cens.model,...){
  riskRegression(formula=formula,data=data,times=times,link="relative",cause=cause,cens.model=cens.model,...)
}
