### print.synth_code.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 14 2025 (06:44) 
## Version: 
## Last-Updated: May 14 2025 (16:08) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Print method for synthesized code
##'
##' @title Print synthesized code
##' @param x Object obtained with \code{synthesize}
##' @param digits Number of digits. Not used.
##' @param ... Not used.
#'
#' method print synth_code
#' @export
print.synth_code <- function(x,digits,...){
    cat(x)
}


######################################################################
### print.synth_code.R ends here
