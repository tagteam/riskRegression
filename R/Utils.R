### Utils.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 19 2019 (15:52) 
## Version: 
## Last-Updated: sep 19 2019 (16:08) 
##           By: Brice Ozenne
##     Update #: 11
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
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
#' M <- matrix(letters,nrow = 26, ncol = 2)
#' rowPaste(M)
#' }
rowPaste <- function(object){
    if(is.matrix(object)){
        return(    do.call("paste0",lapply(1:NCOL(object), function(iC){object[,iC]}))  )
    }else if(is.list(x) || is.data.frame(x)){
        return( do.call("paste0",object ) )
    }else{
        stop("Arugment \'object\' must be a matrix, data.frame, or list \n")
    }
    

}


## * subsetCols
#' @title Extract Columns From Matrix
#' @description Extract columns from matrix.
#'
#' @param object A matrix
#' @param index index of the columns to be extracted.
#' 0 indicates that the column should be set to the default value.
#' NA indicates that the column should be set to NA.
#' @param default the default value.
#'
#' @examples
#' M <- matrix(rnorm(50),5,10)
#' subsetCols(M, index = c(0,0,1), default = 0)
#' subsetCols(M, index = c(0,2,3,NA), default = 0)
#' subsetCols(M, index = c(0,NA,2,3,NA), default = 0)
subsetCols <- function(object, index, default){

    ## if(!is.matrix(object)){
    ##     stop("Argument \'object\' must be a matrix \n")
    ## }
    out <- cbind(default,object)[,index+1,drop = FALSE]
    dimnames(out) <- NULL
    return(out)

    
    ## n.index <- length(index)
    ## n.row <- NROW(object)

    ## pos0 <- which(index==0)
    ## posNA <- which(is.na(index))
    ## posValid <- setdiff(1:n.index,c(pos0,posNA))

    ## new.object <- matrix(default, nrow = n.row, ncol = n.index)
    ## if(length(posNA)>0){
    ##     new.object[,posNA] <- NA
    ## }
    ## if(length(posValid)>0){
    ##     new.object[,posValid] <- object[,index[posValid]]
    ## }
    ## return(new.object)
}

######################################################################
### Utils.R ends here
