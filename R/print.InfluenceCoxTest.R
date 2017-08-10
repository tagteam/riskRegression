### print.InfluenceCoxTest.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jul 20 2017 (11:05) 
## Version: 
## last-updated: jul 20 2017 (11:08) 
##           By: Brice Ozenne
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Print the results of the influence test
#' @description Print the results of the influence test
#'
#' @param x Object obtained with influenceCoxTest
#' @param ... Passed to print
#'
#' @method print influenceCoxTest
#' @export
print.influenceCoxTest <- function(x,...){
  
  ## transformation
  # if(class(x$transformation)=="function"){
  #   message <- c(paste0("\nWARNING: values (delta,inf,sup",if(!is.na(x$pBand)){",infBand,supBand"},") are not given on the original scale \n"),
  #                paste0("transformation: ",attr(x$transformation,"srcref"),"\n"))
  #   
  # }else{
  #   message <- "\nValues are given on the original scale \n"
  # }
  message <- NULL
  
  ## band
  keep.name <- c("time","delta","inf","sup")
  if(!is.na(x$pBand)){
    keep.name <- c(keep.name,"infBand","supBand")
    message <- c("No difference in estimate across all times : p = ",x$pBand,"\n",
                 message)
  }else{
    messageFinal <- NULL
  }
  keep.name <- c(keep.name,"p")
  
  ## display
  print(x$table[,keep.name])
  cat(message, sep = "")
  
  ## export
  return(invisible(x))
}

#----------------------------------------------------------------------
### print.InfluenceCoxTest.R ends here
