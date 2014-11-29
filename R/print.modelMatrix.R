#' @S3method print modelMatrix
print.modelMatrix <- function(x,...){
  levs <- attr(x,"factorLevels")
  refs <- attr(x,"refLevels")
  sargs <- attr(x,"specialArguments")
  varnames <- attr(x,"variableNames")
  cat("\nModel matrix for riskRegression regression models\n")
  numvars <- varnames[!match(varnames,names(levs),nomatch=FALSE)]
  facvars <- names(levs)
  cat("\nNumeric variables:\n")
  if (length(numvars)>0){
    cat(paste(numvars,collapse=", "),"\n")
  }
  else{
    cat("\nNone\n")
  }
  cat("\nFactor variables:\n")
  if (length(facvars)>0){
    fmat <- cbind("Varname"=facvars,"refLevel"=refs[facvars])
    rownames(fmat) <- NULL
    print(fmat)
  }
  else{
    cat("\nNone\n")
  }
  if (!is.null(sargs)){
    cat("\nSpecial arguments:\n")
    print(sargs)
  }
  cat("\nHead of matrix:\n")
  class(x) <- "matrix"
  print(head(x))
}
