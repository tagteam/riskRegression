### as.data.frame.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  9 2026 (16:13) 
## Version: 
## Last-Updated: nov  9 2026 (16:23) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * as.data.frame.ate (code)
#' @export
as.data.frame.ate <- function(x, row.names = NULL, optional = FALSE, ...){

    return(as.data.frame(data.table::as.data.table(x = x, ...)))
    
}


##----------------------------------------------------------------------
### as.data.frame.ate.R ends here
