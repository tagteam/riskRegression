##' The call is added to an ctree object
##'
##' @title S3-Wrapper for ctree.
##' @param ... passed to ctree
##' @return list with two elements: ctree and call
##' @seealso Cforest
##' @examples
##' if (require("party",quietly=TRUE)){
##' library(prodlim)
##' library(party)
##' library(survival)
##' set.seed(50)
##' d <- SimSurv(50)
##' nd <- data.frame(X1=c(0,1,0),X2=c(-1,0,1))
##' f <- Ctree(Surv(time,status)~X1+X2,data=d)
##' predictRisk(f,newdata=nd,times=c(3,8))
##' }
##'
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
##' @export
Ctree <- function(...) {
  out <- list(ctree = party::ctree(...))
  class(out) <- "Ctree"
  out$call <- match.call()
  out
}

##' @export
predictRisk.Ctree <- function(object, newdata, times, ...) {
  requireNamespace("party")
  # N <- NROW(newdata)
  # NT <- length(times)
  data.table::setDF(newdata)
  survObj <- party::treeresponse(object$ctree, newdata = newdata)
  p <- do.call("rbind", lapply(survObj, function(x) {
    predictRisk(x, newdata = newdata[1, , drop = FALSE], times = times)
  }))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
      stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    }
  }
  p
}

# CFOREST
# --------------------------------------------------------------------
#' S3-wrapper function for cforest from the party package
#'
#' S3-wrapper function for cforest from the party package
#'
#' See \code{cforest} of the \code{party} package.
#'
#' @param formula Passed on as is. See \code{cforest} of the \code{party} package
#' @param data Passed on as is. See \code{cforest} of the \code{party} package
#' @param ... Passed on as they are. See \code{cforest} of the \code{party} package
#' @return list with two elements: cforest and call
#' @references Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23. URL
#' http://www.jstatsoft.org/v50/i11/.
#' @keywords survival
#' @export Cforest
Cforest <- function(formula, data, ...) {
  requireNamespace("party")
  out <- list(forest = party::cforest(formula, data, ...))
  class(out) <- "Cforest"
  out$call <- match.call()
  out
}


##' @export
predictRisk.Cforest <- function(object, newdata, times, ...) {
  requireNamespace("party")
  if (missing(times) || is.null(times)) {
    p <- as.numeric(unlist(party::treeresponse(object$forest, newdata = newdata)))
    return(p)
  } else {
    survObj <- party::treeresponse(object$forest, newdata = newdata)
    p <- do.call("rbind", lapply(survObj, function(x) {
      predictRisk(x, newdata = newdata[1, , drop = FALSE], times = times)
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
      stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ", NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    }
    p
  }
}
