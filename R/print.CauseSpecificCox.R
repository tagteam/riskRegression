#' Print of a Cause-Specific Cox regression model
#'
#' Print of a Cause-Specific Cox regression model
#' @param x Object obtained with CSC
#' @param ... Passed to print
#'
#' @method print CauseSpecificCox
#' @export
print.CauseSpecificCox <- function(x, ...) {
  print(x$call)
  print(x$response)
  if (x$surv.type == "hazard") {
    lapply(1:length(x$causes), function(c) {
      cat("\n\n----------> Cause: ", x$causes[c], "\n\n")
      xc <- x$models[[c]]
      xc$call$data <- NULL
      if (x$fitter == "coxph") {
        print(summary(xc), ...)
      } else {
        print(xc, ...)
      }
    })
  } else { # surv.type=="survival"
    cat("\n\n----------> Cause: ", x$theCause, "\n\n")
    x1 <- x$models[[1]]
    x1$call$data <- NULL
    print(summary(x1), ...)
    cat("\n\n----------> Event-free survival:\n\n")
    x2 <- x$models[[2]]
    x2$call$data <- NULL
    if (x$fitter == "coxph") {
      print(summary(x2), ...)
    } else {
      print(x2, ...)
    }
  }
}
