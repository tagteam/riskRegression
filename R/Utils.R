### Utils.R ---
## ----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 19 2019 (15:52)
## Version:
## Last-Updated: sep 23 2019 (22:16)
##           By: Brice Ozenne
##     Update #: 18
## ----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
## ----------------------------------------------------------------------
##
### Code:

## * rowPaste
#' @title Collapse Rows of Characters.
#' @description Collapse rows of characters. Fast alternative to apply(x,1,paste0,collapse="")
#'
#' @param object A matrix/data.frame/list containing the characters.
#'
#' @examples
#' \dontrun{
#' M <- matrix(letters, nrow = 26, ncol = 2)
#' rowPaste(M)
#' }
rowPaste <- function(object) {
  if (is.matrix(object)) {
    return(do.call("paste0", lapply(seq_len(NCOL(object)), function(iC) {
      object[, iC]
    })))
  } else if (is.list(object) || is.data.frame(object)) {
    return(do.call("paste0", object))
  } else {
    stop("Arugment \'object\' must be a matrix, data.frame, or list \n")
  }
}

## * subsetIndex
#' @title Extract Specific Elements From An Object
#' @description Extract specific elements from an object.
#' @name subsetIndex
#'
#' @param object A vector or a matrix.
#' @param index index of the elements to be extracted.
#' 0 indicates that the column should be set to the default value.
#' NA indicates that the column should be set to NA.
#' @param default the default value.
#' @param col If object is a matrix, \code{TRUE} lead to extract the columns and \code{FALSE} the rows.
#' @param ... Only used by the generic method.
#'
#' @examples
#' M <- matrix(rnorm(50), 5, 10)
#' subsetIndex(M, index = c(0, 0, 1), default = 0)
#' subsetIndex(M, index = c(0, 2, 3, NA), default = 0)
#' subsetIndex(M, index = c(0, NA, 2, 3, NA), default = 0)
#'
#' C <- 1:10
#' subsetIndex(C, index = c(0, 0, 1, 5, NA), default = 0)
#' @export
`subsetIndex` <- function(object, index, default, ...) UseMethod("subsetIndex")

## ** subsetIndex.default
#' @rdname subsetIndex
#' @export
subsetIndex.default <- function(object, index, default, ...) {
  out <- c(default, object)[index + 1]
  return(out)
}

## ** subsetIndex.matrix
#' @rdname subsetIndex
#' @export
subsetIndex.matrix <- function(object, index, default, col = TRUE, ...) {
  if (col) {
    out <- cbind(default, object)[, index + 1, drop = FALSE]
  } else {
    out <- rbind(default, object)[index + 1, , drop = FALSE]
  }
  dimnames(out) <- NULL
  return(out)
}

######################################################################
### Utils.R ends here
