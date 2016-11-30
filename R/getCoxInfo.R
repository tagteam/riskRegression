#' @title Extract the event times used to fit a Cox model
#' @description Extract the event times used to fit a Cox model
#' @rdname CoxEventtime 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxEventtime(mCox)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxEventtime(mCox)
#' }

#' @rdname CoxEventtime
CoxEventtime <- function(object){
  UseMethod("CoxEventtime") 
} 

#' @rdname CoxEventtime
#' @method CoxEventtime cph
CoxEventtime.cph <- function(object){
  return(object$y[,"time"])
}

#' @rdname CoxEventtime
#' @method CoxEventtime coxph
CoxEventtime.coxph <- function(object){
  return(object$y[,"time"])
}

#' @title Extract data
#' @description Extract the columns from a dataset involved by a model. Optional argument to center the columns values.
#' @name CoxData
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param lpVars the name of the covariates 
#' @param stratavars the name of the strata variables
#' @param center should the covariates be centered
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @details if no variables for the linear predictor or the strata, returns a data.table with one column filled with 0.
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCoxS <- coxph(Surv(time, event) ~ X1, data = d) 
#' CoxData(mCoxS, lpVars = "X1", stratavars = NULL, center = FALSE)
#' mCoxS <- coxph(Surv(time, event) ~ X1*X2, data = d) 
#' CoxData(mCoxS, lpVars = c("X1","X2"), stratavars = NULL, center = FALSE)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d) 
#' CoxData(mCoxS, lpVars = "X2", stratavars = "strata(X1)", center = FALSE)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d) 
#' CoxData(mCoxS, lpVars = NULL, stratavars = c("strata(X1)","strata(X2)"), center = FALSE)
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ X1, data = d) 
#' CoxData(mCoxS, lpVars = "X1", data = d, stratavars = NULL, center = FALSE)
#' mCoxS <- cph(Surv(time, event) ~ X1*X2, data = d) 
#' CoxData(mCoxS, lpVars = c("X1","X2"), data = d, stratavars = NULL, center = FALSE)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d) 
#' xx <- CoxData(mCoxS, lpVars = "X2", data = d, stratavars = "X1", center = FALSE)
#' print(xx)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d) 
#' xx <- CoxData(mCoxS, lpVars = NULL, data = d, stratavars = c("X1","X2"), center = FALSE)
#' print(xx)
#' }
CoxData <- function(object, data = model.frame(object), lpVars, stratavars, center){
 
  if(!is.data.table(data)){
    data <- as.data.table(data[,c(lpVars,stratavars), drop = FALSE])
  }
  modeldata <- NULL
  
  # variables for the linear predictor
  if(length(lpVars)>0){
    modeldata <- data[,lpVars,with=FALSE]
    
    if(center){
      modeldata <- sweep(modeldata, FUN = "-", MARGIN = 2, STATS = object$means[lpVars])
      modeldata <- as.data.table(modeldata)
    }
    
  }
  
  # variables for the strata
  if(length(stratavars)>0){
     
    if(is.null(modeldata)){
      modeldata <- data[,stratavars,with=FALSE]
    }else{
      modeldata[,stratavars := data[,stratavars,with=FALSE], with = FALSE]
    }
    
  }
  
  if(is.null(modeldata)){
    
    n <- if(NROW(data)==0) CoxN(object) else NROW(data)
    modeldata <- data.table(matrix(0, nrow = n, ncol = 1))
    
  }
  
  return(modeldata)
  
} 

#' @title Compute the linear predictor of a Cox model
#' @description Compute the linear predictor of a Cox model
#' @rdname CoxLP 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param lpVars the name of all the covariates
#' @param stratavars the variables used to define the strata
#' @param center should the linear predictor be computed after centering the covariates
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @details In case of empty linear predictor returns a vector of 0 with the same length as the number of rows of the dataset
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxLP(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL, center = FALSE)
#' CoxLP(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL, center = TRUE)  
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d)
#' CoxLP(mCoxS, data = d, xvars = c("X1","X2"), stratavars = "strata(X1)", center = FALSE) 
#' mCoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' CoxLP(mCoxS2, data = d, xvars = c("X1","X2"), stratavars = c("X1","X2"), center = FALSE) 
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxLP(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL, center = FALSE) 
#' CoxLP(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL, center = TRUE) 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d)
#' CoxLP(mCoxS, data = d, xvars = c("X1","X2"), stratavars = "strat(X1)", center = FALSE) 
#' mCoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' CoxLP(mCoxS2, data = d, xvars = c("X1","X2"), stratavars = c("X1","X2"), center = FALSE) 
#' }

#' @rdname CoxLP
CoxLP <- function(object, data, lpVars, stratavars, center){
  UseMethod("CoxLP") 
} 

#' @rdname CoxLP
#' @method CoxLP cph
CoxLP.cph <- function(object, data = NULL, lpVars, stratavars, center){
  
  if(is.null(data)){ ## training dataset
    
    Xb <- object$linear.predictors
    if(center == FALSE){
      Xb <- Xb + sum(object$means*stats::coef(object))
    }
    
  }else{ ## new dataset
    if(length(lpVars)>0){
      data <- CoxData(object, data = data, lpVars = lpVars, stratavars = stratavars, center = center) # here the center argument will enable to center the X variables
      Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
    }else{ 
      Xb <- rep(0, NROW(data)) 
    } 
  }
  
  return(Xb)
}

#' @rdname CoxLP
#' @method CoxLP coxph
CoxLP.coxph <- function(object, data = NULL, lpVars, stratavars, center){
  
  if(is.null(data)){ ## training dataset
    
    Xb <- object$linear.predictors
    if(center == FALSE){
      Xb <- Xb + sum(object$means*stats::coef(object))
    }
    
  }else{ ## new dataset
    
    if(length(lpVars) == 0){ 
      Xb <- rep(0, NROW(data))
    } else{
      data <- CoxData(object, data = data, lpVars = lpVars, stratavars = stratavars, center = center) # here the center argument will enable to center the X variables
      if(length(stratavars)>0){
        Xb <- rowSums(stats::predict(object, newdata = as.data.frame(data), type = "terms"))
      }else { 
        Xb <- stats::predict(object, newdata = data.frame(data), type = "lp")
      }
    } 
    
  }
  
  return(Xb)
}

#' @title Extract the number of observations from a Cox model
#' @description Extract the number of observations from a Cox model
#' @rdname CoxN 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxN(mCox)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxN(mCox)
#' }

#' @rdname CoxN
CoxN <- function(object){
  UseMethod("CoxN") 
} 

#' @rdname CoxN
#' @method CoxN cph
CoxN.cph <- function(object){
  return(sum(object$n))
}

#' @rdname CoxN
#' @method CoxN coxph
CoxN.coxph <- function(object){
  return(object$n)
}


#' @title Extract the time and event variable from a Cox mode
#' @description Extract the time and event variable from a Cox mode
#' @rdname CoxResponseVar 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxResponseVar(mCox)
#'
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxResponseVar(mCox) 
#' }
CoxResponseVar <- function(object){

  f.surv <- object$call$formula[[2]]
  if("time" %in% names(f.surv) == FALSE && "event" %in% names(f.surv) == FALSE ){
    names(f.surv)[2:3] <- c("time","event")
  }else if("time" %in% names(f.surv) == FALSE){
    names(f.surv)[setdiff(2:3,which(names(f.surv)=="event"))] <- "time"
  }else if("event" %in% names(f.surv) == FALSE){
    names(f.surv)[setdiff(2:3,which(names(f.surv)=="time"))] <- "event"
  }
  
  return(list(time = deparse(f.surv$time),
              event = deparse(f.surv$event))
  )
  
}

#' @title Extract the event indicator used to fit a Cox model
#' @description Extract the event indicator used to fit a Cox model
#' @rdname CoxStatus 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxStatus(mCox)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxStatus(mCox)
#' }

#' @rdname CoxStatus
CoxStatus <- function(object){
  UseMethod("CoxStatus") 
} 

#' @rdname CoxStatus
#' @method CoxStatus cph
CoxStatus.cph <- function(object){
  return(object$y[,"status"])
}

#' @rdname CoxStatus
#' @method CoxStatus coxph
CoxStatus.coxph <- function(object){
  return(object$y[,"status"])
}

#' @title Define the strata for a new dataset
#' @description Define the strata in a dataset to match those of a stratified Cox model
#' @name CoxStrata
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param sterms terms in the formula corresponding to the strata variables
#' @param stratavars the name of the variables used to define the strata
#' @param levels the strata levels that have been used to fit the Cox model
#' @param stratalevels a named list containing for each variable used to form the strata all its possible levels
#'
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @details if no strata variables returns a vector of \code{"1"} (factor).
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' resInfo <- CoxStrataVar(mCoxS)
#' Ostrata <- CoxStrata(mCoxS, stratavars = resInfo$strata$stratavars)
#' CoxStrata(mCoxS, data = d, sterms = resInfo$strata$sterms, stratavars = resInfo$strata$stratavars, 
#'           levels = levels(Ostrata), stratalevels = resInfo$strata$stratalevels)
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' resInfo <- CoxStrataVar(mCoxS)
#' Ostrata <- CoxStrata(mCoxS, stratavars = resInfo$strata$stratavars)
#' CoxStrata(mCoxS, data = d, sterms = resInfo$strata$sterms, stratavars = resInfo$strata$stratavars, 
#'           levels = levels(Ostrata), stratalevels = resInfo$strata$stratalevels)
#' }

#' @rdname CoxStrata
CoxStrata <- function(object, data, sterms, stratavars, levels, stratalevels) UseMethod("CoxStrata")

#' @rdname CoxStrata
#' @method CoxStrata coxph
CoxStrata.cph <- function(object, data = NULL, sterms, stratavars, levels, stratalevels){
  
  if(length(stratavars)==0){ ## no strata variables
    
    n <- if(is.null(data)) CoxN(object) else NROW(data)
    strata <- as.factor(rep("1", n))
    
  }else{  ## strata variables
    
    if(is.null(data)){ ## training dataset
      strata <- object$strata
    }else { ## new dataset
      tmp <- model.frame(sterms,data)
      colnames(tmp) <- names(prodlim::parseSpecialNames(names(tmp),"strat"))
      tmp <- data.frame(lapply(1:NCOL(tmp),function(j){factor(paste0(names(tmp)[j],"=",tmp[,j,drop=TRUE]))}))
      strata <- apply(tmp,1,paste,collapse=".")
      if (any(unique(strata) %in% levels == FALSE)){
        stop("unknown strata: ",paste(unique(strata[strata %in% levels == FALSE]), collapse = " | "),"\n")
      }
      strata <- factor(strata, levels = levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert strata to numeric
    }
    
  }
  
  return(strata)
}

#' @rdname CoxStrata
CoxStrata.coxph <- function(object, data = NULL, sterms, stratavars, levels, stratalevels){
  
  if(length(stratavars)==0){ ## no strata variables
    
    n <- if(is.null(data)) CoxN(object) else NROW(data)
    strata <- as.factor(rep("1", n))
    
  }else{  ## strata variables
    
    if(is.null(data)){ ## training dataset
      strata <- interaction(stats::model.frame(object)[,stratavars], drop = TRUE, sep = ", ", lex.order = TRUE) 
    }else { ## new dataset
      strata <- prodlim::model.design(sterms,data=data,xlev=stratalevels,specialsFactor=TRUE)$strata[[1]]
      if (any(unique(strata) %in% levels == FALSE)){
        stop("unknown strata: ",paste(unique(strata[strata %in% levels == FALSE]), collapse = " | "),"\n")
      }
      strata <- factor(strata, levels = levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert strata to numeric
    }
    
  }
  
  return(strata)
}


#' @title Extract variable names from a model
#' @description Extract the name of the variables belonging to the linear predictor or used to form the strata
#' @rdname CoxStrataVar 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxStrataVar(mCox)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d)
#' CoxStrataVar(mCoxS)
#' mCoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' CoxStrataVar(mCoxS2)
#' mCoxI <- coxph(Surv(time, event) ~ X1*X2, data = d)
#' CoxStrataVar(mCoxI)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxStrataVar(mCox) 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d)
#' CoxStrataVar(mCoxS) 
#' mCoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' CoxStrataVar(mCoxS2) 
#' mCoxI <- cph(Surv(time, event) ~ X1*X2, data = d)
#' CoxStrataVar(mCoxI)
#' }

#' @rdname CoxStrataVar
CoxStrataVar <- function(object){
  UseMethod("CoxStrataVar") 
} 

#' @rdname CoxStrataVar
#' @method CoxStrataVar coxph
CoxStrataVar.coxph <- function(object){
  
  xterms <- delete.response(object$terms)
  
  # renamed variables (e.g. strata(X1))
  xvars <- attr(xterms,"term.labels")
  strataspecials <- attr(xterms,"specials")$strata
  stratavars <- xvars[strataspecials]
  
  # original variables (e.g. X1)
  allVars.X <- all.vars(xterms)
  stratavars.original <- allVars.X[strataspecials]
  
  is.strata <- length(strataspecials)>0
  
  if(is.strata){
    if (length(xvars)>length(strataspecials)) 
      sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
    else 
      sterms <- xterms
    stratalevels <- object$xlevels[stratavars]
  }else{
    sterms <- NULL
    stratalevels <- NULL
  }
  
  ## export
  return(list(stratavars = stratavars,
              stratavars.original = stratavars.original,
              strataspecials = strataspecials,
              stratalevels = stratalevels,
              is.strata = is.strata,
              sterms = sterms,
              lpvars = setdiff(allVars.X,stratavars.original))
         )
}

#' @rdname CoxStrataVar
#' @method CoxStrataVar cph
CoxStrataVar.cph <- function(object){
  
  xterms <- delete.response(object$terms)
  
  # renamed variables (e.g. strat(X1))
  xvars <- attr(xterms,"term.labels")
  strataspecials <- attr(xterms,"specials")$strat
  stratavars <- xvars[strataspecials]
  
  # original variables (e.g. X1)
  allVars.X <- all.vars(xterms)
  stratavars.original <- allVars.X[strataspecials]
  
  is.strata <- length(strataspecials)>0
  
  if(is.strata){ ## cph:strata for estimation of the baseline hazard
    if (length(xvars)>length(strataspecials)) 
      sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
    else 
      sterms <- xterms
    stratavars <- xvars[strataspecials]
  }else{
    sterms <- NULL
  }
  
  ## export
  return(list(stratavars = stratavars,
              stratavars.original = stratavars.original,
              strataspecials = strataspecials,
              stratalevels = NULL,
              is.strata = is.strata,
              sterms = sterms,
              lpvars = setdiff(allVars.X,stratavars.original))
  )
}

