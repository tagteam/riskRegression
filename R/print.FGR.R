#' @S3method print FGR
print.FGR <- function(x,...){
  summary(x$response)
  cat("\n\nFine-Gray model: analysis of cause",x$cause,"\n")
  f <- x$crrFit
  cat("\n")
  f$call <- x$call
  print(summary(f,...))
  cat(paste("\nConvergence:",f$converged),"\n\n")
}
