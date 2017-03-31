#' Extract coefficients from a Cause-Specific Cox regression model
#' 
#' Extract coefficients from a Cause-Specific Cox regression model
#' @param x Object obtained with CSC
#'
#' @method coef CauseSpecificCox
#' @export
coef.CauseSpecificCox <- function(x){
  return(lapply(x$models,coef))
}
