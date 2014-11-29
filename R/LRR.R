#' @export
LRR <- function(formula,data,times,cause,cens.model,cens.formula,...){
  fit <- riskRegression(formula=formula,
                        data=data,
                        times=times,
                        link="logistic",
                        cause=cause,
                        cens.model=cens.model,
                        cens.formula=cens.formula,
                        ...)
  fit$call <- match.call()
  fit
}
