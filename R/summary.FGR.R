summary.FGR <- function(object,...){
  f <- object$crrFit
  f$call <- object$call
  summary(f,...)
}
