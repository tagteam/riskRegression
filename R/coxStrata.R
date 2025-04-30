### coxStrata.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:35) 
## Version: 
## Last-Updated: Apr 29 2025 (06:51) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Define the strata for a new dataset
#' @description Define the strata in a dataset to match those of a stratified Cox model
#' @name coxStrata
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param sterms terms in the formula corresponding to the strata variables
#' @param strata.vars the name of the variables used to define the strata
#' @param levels the strata levels that have been used to fit the Cox model
#' @param strata.levels a named list containing for each variable used to form the strata all its possible levels
#'
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @details if no strata variables returns a vector of \code{"1"} (factor).

#' @rdname coxStrata
#' @export
coxStrata <- function(object, data, sterms, strata.vars, strata.levels) UseMethod("coxStrata")

## ** coxStrata.cph
#' @rdname coxStrata
#' @method coxStrata cph
#' @export
coxStrata.cph <- function(object, data, sterms, strata.vars, strata.levels){
  
  if(length(strata.vars)==0){ ## no strata variables
    
    n <- if(is.null(data)) coxN(object) else NROW(data)
    strata <- as.factor(rep("1", n))
    
  }else{  ## strata variables

      if(is.null(data)){ ## training dataset
          strata <- object$strata
      }else { ## new dataset
          tmp <- model.frame(sterms,data)
          colnames(tmp) <- names(prodlim::parseSpecialNames(names(tmp),"strat"))
          tmp <- data.frame(lapply(1:NCOL(tmp),function(j){factor(paste0(names(tmp)[j],"=",tmp[,j,drop=TRUE]))}))
          strata <- apply(tmp,1,paste,collapse=".")
          if (any(unique(strata) %in% strata.levels == FALSE)){
              stop("unknown strata: ",paste(unique(strata[strata %in% strata.levels == FALSE]), collapse = " | "),"\n")
          }
          strata <- factor(strata, levels = strata.levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert strata to numeric
      }
    
  }
  
  return(strata)
}

## ** coxStrata.coxph
#' @rdname coxStrata
#' @method coxStrata coxph
#' @export
coxStrata.coxph <- function(object, data, sterms, strata.vars, strata.levels){
  
  if(length(strata.vars)==0){ ## no strata variables
    
    n <- if(is.null(data)) coxN(object) else NROW(data)
    strata <- as.factor(rep("1", n))
    
  }else{  ## strata variables
    
    if(is.null(data)){ ## training dataset
      strata <- object$strata
    }else { ## new dataset
      strata <- prodlim::model.design(sterms,data=data,xlev=strata.levels,specialsFactor=TRUE)$strata[[1]]
      if (any(unique(strata) %in% strata.levels == FALSE)){
        stop("unknown strata: ",paste(unique(strata[strata %in% strata.levels == FALSE]), collapse = " | "),"\n")
      }
      strata <- factor(strata, levels = strata.levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert strata to numeric
    }
    
  }
  return(strata)
}

## ** coxStrata.phreg
#' @rdname coxStrata
#' @method coxStrata phreg
#' @export
coxStrata.phreg <- function(object, data, sterms, strata.vars, strata.levels){
  
  if(length(strata.vars)==0){ ## no strata variables
    
    n <- if(is.null(data)) coxN(object) else NROW(data)
    strata <- as.factor(rep("1", n))
    
  }else{  ## strata variables
      if(is.null(data)){ ## training dataset
          strata <- object$model.frame[[object$strata.name]]
      }else { ## new dataset
          strata <- prodlim::model.design(sterms,data=data,xlev=strata.levels,specialsFactor=TRUE)$strata[[1]]
          if (any(unique(strata) %in% strata.levels == FALSE)){
              stop("unknown strata: ",paste(unique(strata[strata %in% strata.levels == FALSE]), collapse = " | "),"\n")
          }
      }
    
  }
  return(strata)
}
# }}}

## ** coxStrata.prodlim
#' @rdname coxStrata
#' @method coxStrata prodlim
#' @export
coxStrata.prodlim <- function(object, data, sterms, strata.vars, strata.levels){

  if(length(strata.vars)==0){ ## no strata variables
    
    n <- if(is.null(data)) coxN(object) else NROW(data)
    strata <- as.factor(rep("1", n))
    
  }else{  ## strata variables

      if(is.null(data)){ ## training dataset
          dt.X <- as.data.table(object$model.matrix[object$originalDataOrder,,drop=FALSE])
          strata <- interaction(dt.X, sep = ", ")
      }else { ## new dataset
          strata <- interaction(as.data.table(data)[,.SD,.SDcols = attr(sterms,"term.labels")], sep = ", ")
          if (any(levels(strata) %in% strata.levels == FALSE)){
              stop("unknown strata: ",paste(unique(strata[strata %in% strata.levels == FALSE]), collapse = " | "),"\n")
          }
          strata <- factor(strata, levels = strata.levels)
      }
    
  }
    return(strata)
}

## ** coxStrata.coxnet
#' @rdname coxStrata
#' @method coxStrata coxnet
#' @export
coxStrata.coxnet <- function(object, data, sterms, strata.vars, strata.levels){
  
  if(length(strata.vars)==0){ ## no strata variables
    
    n <- if(is.null(data)) coxN(object) else NROW(data)
    strata <- as.factor(rep("1", n))
    
  }else{  ## strata variables
    
    if(is.null(data)){ ## training dataset
      strata <- object$strata
    }else { ## new dataset
      strata <- prodlim::model.design(sterms,data=data,xlev=strata.levels,specialsFactor=TRUE)$strata[[1]]
      if (any(unique(strata) %in% strata.levels == FALSE)){
        stop("unknown strata: ",paste(unique(strata[strata %in% strata.levels == FALSE]), collapse = " | "),"\n")
      }
      strata <- factor(strata, levels = strata.levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert strata to numeric
    }
    
  }
  return(strata)
}

######################################################################
### coxStrata.R ends here
