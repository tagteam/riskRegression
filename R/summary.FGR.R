#' Summary of a Fine-Gray regression model
#'
#' Summary of a Fine-Gray regression model
#' @param object Object fitted with function FGR 
#' @param ... passed to cmprsk::summary.crr
#'
#' @method summary FGR
#' @export
summary.FGR <- function(object,...){
  f <- object$crrFit
  f$call <- object$call
  summary(f,...)
}
