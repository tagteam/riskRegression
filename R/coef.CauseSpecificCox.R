#' Extract coefficients from a Cause-Specific Cox regression model
#' 
#' Extract coefficients from a Cause-Specific Cox regression model
#' @param object Object obtained with CSC
#' @param ... not used
#'
#' @method coef CauseSpecificCox
#' @export
coef.CauseSpecificCox <- function(object, ...){
  return(lapply(object$models,coef))
}
