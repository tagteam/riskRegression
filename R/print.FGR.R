#' Print of a Fine-Gray regression model
#'
#' Print of a Fine-Gray regression model
#' @param x Object fitted with function FGR 
#' @param ... passed to cmprsk::summary.crr
#'
#' @method print FGR
#' @export
print.FGR <- function(x,...){
  summary(x$response)
  cat("\n\nFine-Gray model: analysis of cause",x$cause,"\n")
  f <- x$crrFit
  cat("\n")
  f$call <- x$call
  print(summary(f,...))
  cat(paste("\nConvergence:",f$converged),"\n\n")
}
