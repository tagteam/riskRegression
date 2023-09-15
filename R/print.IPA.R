### print.IPA.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Nov  4 2019 (09:07)
## Version:
## Last-Updated: Jun 25 2020 (13:42)
##           By: Thomas Alexander Gerds
##     Update #: 20
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
##' Print method for IPA
##'
##' @title Print IPA object
##' @param x Object obtained with \code{IPA}
##' @param digits Number of digits
##' @param ... passed to print
##'
##' @method print IPA
##' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
print.IPA <- function(x, digits = 2, ...) {
  Brier <- IPA.drop <- IPA <- NULL
  if (missing(digits)) {
    digits <- 1
  }
  X <- copy(x)
  data.table::setDT(X)
  fmt <- paste0("%1.", digits[[1]], "f")
  X[, Brier := sprintf(fmt = fmt, 100 * Brier)]
  X[, IPA := sprintf(fmt = fmt, 100 * IPA)]
  if (match("IPA.drop", colnames(X), nomatch = 0)) X[, IPA.drop := sprintf(fmt = fmt, 100 * IPA.drop)]
  print(X, ...)
  message("\nNOTE: Values are multiplied by 100 and given in %.")
  message("NOTE: The higher IPA the better.")
  message("NOTE: IPA.drop = IPA(Full model) - IPA. The higher the drop\nthe more important is the variable for the full model.")
}


######################################################################
### print.IPA.R ends here
