#### functions #####

#' @title Extract the design matrix for a Cox model
#' @description Extract the design matrix for a Cox model
#' @name CoxDesign
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param lpvars the variables used to define the linear predictor
#' @param stratavars the special variables used to define the strata
#' @param rm.intercept should the intercept be removed
#' @param center should the covariates be centered
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
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
#' info <- CoxStrataVar(mCoxS)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = FALSE)
#' CoxDesign(mCoxS, data = d, lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = TRUE)
#' 
#' mCoxS <- coxph(Surv(time, event) ~ X1*X2, data = d) 
#' info <- CoxStrataVar(mCoxS)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = FALSE)
#' 
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d) 
#' info <- CoxStrataVar(mCoxS)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = FALSE)
#' 
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d) 
#' info <- CoxStrataVar(mCoxS)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = FALSE)
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ X1, data = d) 
#' info <- CoxStrataVar(mCoxS)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = FALSE)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = TRUE)
#' 
#' mCoxS <- cph(Surv(time, event) ~ X1*X2, data = d) 
#' info <- CoxStrataVar(mCoxS)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = FALSE)
#' 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d) 
#' info <- CoxStrataVar(mCoxS)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = FALSE)
#' 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d) 
#' info <- CoxStrataVar(mCoxS)
#' CoxDesign(mCoxS, data = CoxData(mCoxS), lpvars = info$lpvars, stratavars = info$stratavars, rm.intercept = TRUE, center = FALSE)
#' }
CoxDesign <- function(object, data, lpvars, stratavars, rm.intercept, center){
  
  f <- CoxFormula(object)
  
  # remove strata
  if(length(stratavars)>0){
    f.strata <- as.formula(paste0(".~ .-",paste(stratavars, collapse = "-")))
    f <- update(f,
                f.strata)
  }
  
  # remove Surv
  terms.f <- delete.response(terms(f))
  
  # remove intercept
  if(rm.intercept == TRUE){attr(terms.f, "intercept") <- 0}
  
  # add lp variables
  modeldata <- as.data.frame(model.matrix(terms.f, 
                                          data = data, 
                                          contrasts = NULL))
  
  # center data
  if(center && length(lpvars)>0){
    modeldata[,lpvars] <- rowCenter_cpp(as.matrix(modeldata[,lpvars,drop = FALSE]), center = CoxCenter(object))
  }
  
  # # add strata variables
  # if(rm.stratavars == FALSE && length(stratavars)>0){
  #   f.strata <- as.formula(paste0("~0+",paste(stratavars.original, collapse = "+")))
  #   
  #   modeldata <- cbind(modeldata,
  #                      model.frame(formula = f.strata, 
  #                                  data = data, contrasts = NULL)
  #   )
  #   
  # }
  
  if(NCOL(modeldata) == 0){
    modeldata <- cbind(modeldata, 0)
  }
  
  return(modeldata)
} 

#' @title Extract the time and event variable from a Cox model
#' @description Extract the time and event variable from a Cox model
#' @rdname CoxResponseVar 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
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
  
  f.surv <- CoxFormula(object)[[2]]
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

#### methods #####


#' @title Extract the mean value of the covariates
#' @description Extract the mean value of the covariates
#' @rdname CoxCenter 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$Xcat <- paste0(d$X3,d$X4)
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxCenter(mCox)
#' mCox <- coxph(Surv(time, event) ~ X1+X2+Xcat, data = d)
#' CoxCenter(mCox)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxCenter(mCox)
#' mCox <- cph(Surv(time, event) ~ X1+X2+Xcat, data = d)
#' CoxCenter(mCox)
#' }

#' @rdname CoxCenter
CoxCenter <- function(object){
  UseMethod("CoxCenter") 
} 

#' @rdname CoxCenter
#' @method CoxCenter cph
CoxCenter.cph <- function(object){
  return(setNames(object$means, object$mmcolnames))
}

#' @rdname CoxCenter
#' @method CoxCenter coxph
CoxCenter.coxph <- function(object){
  return(object$means)
}

#' @title Extract the dataset used to train a Cox model
#' @description Extract the dataset used to train a Cox model
#' @rdname CoxData 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxData(mCox)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxData(mCox)
#' }

#' @rdname CoxData
CoxData <- function(object){
  UseMethod("CoxData") 
} 

#' @rdname CoxData
#' @method CoxData cph
CoxData.cph <- function(object){
  return(eval(object$call$data))
}

#' @rdname CoxData
#' @method CoxData coxph
CoxData.coxph <- function(object){
  return(eval(object$call$data))
}

#' @title Extract the event times used to fit a Cox model
#' @description Extract the event times used to fit a Cox model
#' @rdname CoxEventtime 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
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

#' @title Extract the formula from a Cox model
#' @description Extract the formula from a Cox model
#' @rdname CoxFormula 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxFormula(mCox)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxFormula(mCox)
#' }

#' @rdname CoxFormula
CoxFormula <- function(object){
  UseMethod("CoxFormula") 
} 

#' @rdname CoxFormula
#' @method CoxFormula cph
CoxFormula.cph <- function(object){
  return(object$sformula)
}

#' @rdname CoxFormula
#' @method CoxFormula coxph
CoxFormula.coxph <- function(object){
  return(object$formula)
}

#' @title Compute the linear predictor of a Cox model
#' @description Compute the linear predictor of a Cox model
#' @rdname CoxLP 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param center should the linear predictor be computed after centering the covariates
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
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
#' CoxLP(mCox, data = d, center = FALSE)
#' CoxLP(mCox, data = d, center = TRUE)  
#' 
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d)
#' CoxLP(mCoxS, data = d, center = FALSE) 
#' 
#' mCoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' CoxLP(mCoxS2, data = d, center = FALSE) 
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxLP(mCox, data = d, center = FALSE) 
#' CoxLP(mCox, data = d, center = TRUE) 
#' 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d)
#' CoxLP(mCoxS, data = d, center = FALSE) 
#' 
#' mCoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' info <- CoxStrataVar(mCoxS2)
#' CoxLP(mCoxS2, data = d, center = FALSE) 
#' }

#' @rdname CoxLP
CoxLP <- function(object, data, center){
  UseMethod("CoxLP") 
} 

#' @rdname CoxLP
#' @method CoxLP cph
CoxLP.cph <- function(object, data, center){
  
  coef <- stats::coef(object)
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    
    Xb <- object$linear.predictors
    
    if(center == FALSE && n.varLP != 0){
      Xb <- Xb + sum(CoxCenter(object)*coef)
    }
    
  }else{ ## new dataset
    
    if(n.varLP>0){
      
      Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
      
      if(center && n.varLP != 0){
        Xb <- Xb - sum(CoxCenter(object)*coef)
      }
      
    }else{ 
      Xb <- rep(0, NROW(data)) 
    } 
  }
  
  return(Xb)
}

#' @rdname CoxLP
#' @method CoxLP coxph
CoxLP.coxph <- function(object, data, center){
  
  coef <- stats::coef(object)
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    
    Xb <- object$linear.predictors
    
    if(center == FALSE && n.varLP != 0){
      Xb <- Xb + sum(CoxCenter(object)*coef)
    }
    
  }else{ ## new dataset
    
    if(n.varLP>0){
      is.strata <- attr(object$terms, "special")$strata
      
      if(length(is.strata)>0){
        Xb <- rowSums(stats::predict(object, newdata = as.data.frame(data), 
                                     type = "terms"))
      }else{ 
        Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
      }
      
      if(center && n.varLP != 0){
        Xb <- Xb - sum(CoxCenter(object)*coef)
      }
      
    }else{ 
      Xb <- rep(0, NROW(data)) 
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
#' @author Brice Ozenne broz@@sund.ku.dk
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

#' @title Extract the event indicator used to fit a Cox model
#' @description Extract the event indicator used to fit a Cox model
#' @rdname CoxStatus 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
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
#' @author Brice Ozenne broz@@sund.ku.dk
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
#' @author Brice Ozenne broz@@sund.ku.dk
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
              lpvars = setdiff(xvars, stratavars))
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
              lpvars = setdiff(xvars, stratavars)
  ))
}

#' @title Extract the variance covariance matrix of the beta from a Cox model
#' @description Extract the variance covariance matrix of the beta from a Cox model
#' @rdname CoxVarCov 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @details Should return \code{NULL} if the Cox model has no covariate. 
#' The rows and columns of the variance covariance matrix must be named with the names used in the design matrix.
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxVarCov(mCox)
#' mCox <- coxph(Surv(time, event) ~ 1, data = d, y = TRUE)
#' CoxVarCov(mCox)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxVarCov(mCox)
#' mCox <- cph(Surv(time, event) ~ 1, data = d, y = TRUE)
#' CoxVarCov(mCox)
#' 
#' }

#' @rdname CoxVarCov
CoxVarCov <- function(object){
  UseMethod("CoxVarCov") 
} 

#' @rdname CoxVarCov
#' @method CoxVarCov cph
CoxVarCov.cph <- function(object){
  
  Sigma <- object$var
  if(!is.null(Sigma)){
    colnames(Sigma) <- object$mmcolnames
    rownames(Sigma) <- object$mmcolnames
  }
  
  return(Sigma)
}

#' @rdname CoxVarCov
#' @method CoxVarCov coxph
CoxVarCov.coxph <- function(object){
  
  Sigma <- object$var
  if(!is.null(Sigma)){
    coefName <- names(coef(object))
    colnames(Sigma) <- coefName
    rownames(Sigma) <- coefName
  }
  return(Sigma)
}
