### reconstructData.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:39) 
## Version: 
## Last-Updated: Apr 27 2025 (07:39) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
## * reconstructData
#' @title Reconstruct the original dataset
#' @description Reconstruct the original dataset from the elements stored in the coxph object
#' @name reconstructData
#' 
#' @param object a coxph object.
#'
#' @author Brice Ozenne broz@@sund.ku.dk and Thomas A. Gerds tag@@biostat.ku.dk
#'
reconstructData <- function(object){
  
  
  ## combine response variable, regressors and strata variable
  newdata <- as.data.frame(cbind(coxModelFrame(object),splitStrataVar(object)))
  
  ## set response variable to their original names
  infoVar <- coxVariableName(object, model.frame = newdata)
  if(!is.null(infoVar$entry)){
    names(newdata)[names(newdata) == "start"] <- infoVar$entry
  }
  names(newdata)[names(newdata) == "stop"] <- infoVar$time
  names(newdata)[names(newdata) == "status"] <-  infoVar$status
  
  ## reconstruct categorical variables (from binary indicators to factors)
  if(length(object$contrasts)>0){ 
    name.categorical <- names(object$contrasts)
    n.categorical <- length(name.categorical)
    for(iCat in 1:n.categorical){ # iCat <- 1
      iName.categorical <- name.categorical[iCat]
      iContrast <- paste0(iName.categorical,object$xlevels[[iName.categorical]])
      
      ## add reference level
      newdata[[iContrast[1]]] <- 1-rowSums(newdata[,iContrast[-1],drop=FALSE])
      iIndex.level <- apply(newdata[,iContrast,drop=FALSE],1,function(x){which(x==1)})
      newdata[[iName.categorical]] <- object$xlevels[[iName.categorical]][iIndex.level]
    }
  }
  
  return(newdata)
}


######################################################################
### reconstructData.R ends here
