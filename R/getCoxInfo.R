#### functions ####

# {{{ CoxVariableName
#' @title Extract variable names from a model
#' @description Extract the name of the variables belonging to the linear predictor or used to form the strata
#' @rdname CoxVariableName 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
#' d$id <- 1:NROW(d)
#' 
#' ##
#' library(survival)
#' d$X1=factor(d$X1)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d, x = TRUE, y = TRUE)
#' CoxVariableName(mCox)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d, x = TRUE, y = TRUE)
#' CoxVariableName(mCoxS)
#' mCoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d, x = TRUE, y = TRUE)
#' CoxVariableName(mCoxS2)
#' mCoxI <- coxph(Surv(time, event) ~ X1*X2, data = d, x = TRUE, y = TRUE)
#' CoxVariableName(mCoxI)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, x = TRUE, y = TRUE)
#' CoxVariableName(mCox) 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d, x = TRUE, y = TRUE)
#' CoxVariableName(mCoxS) 
#' mCoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, x = TRUE, y = TRUE)
#' CoxVariableName(mCoxS2) 
#' mCoxI <- cph(Surv(time, event) ~ X1*X2, data = d, x = TRUE, y = TRUE)
#' CoxVariableName(mCoxI)
#' 
#' ##
#' library(mets)
#' mCox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' CoxVariableName(mCox) 
#' mCoxS <- phreg(Surv(entry, time, event) ~ strata(X1)+cluster(id)+X2, data = d)
#' CoxVariableName(mCoxS) 
#' mCoxI <- phreg(Surv(entry, time, event) ~ X1*X2, data = d)
#' CoxVariableName(mCoxI)
#' }

#' @rdname CoxVariableName
#' @export
CoxVariableName <- function(object){

    f <- CoxFormula(object)
    special.object <- CoxSpecialStrata(object)
    ## response
    ls.SurvVar <- SurvResponseVar(f)
    ## strata
    xterms <- delete.response(terms(f,
                                    special = special.object,
                                    data = CoxDesign(object)))
    ls.StrataInfo <- extractStrata(xterms,
                                   special = CoxSpecialStrata(object))
  
    ## export
    return(c(list(entry = ls.SurvVar$entry),
             list(time = ls.SurvVar$time),
             list(status = ls.SurvVar$status),
             ls.StrataInfo,
             list(lpvars = names(coef(object)))
             ))
} 
# }}}

CoxCovars <- function(object){
    ttobj <- stats::terms(object)
    ## colnames(attr(ttobj,"factors"))
    all.vars(attr(delete.response(ttobj),"variables"))
}

#### methods #####

# {{{ CoxBaseEstimator
#' @title Extract the type of estimator for the baseline hazard
#' @description Extract the type of estimator for the baseline hazard
#' @rdname CoxBaseEstimator 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' CoxBaseEstimator(mCox)
#'
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' CoxBaseEstimator(mCox)
#'
#' ##
#' library(mets)
#' mCox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' CoxBaseEstimator(mCox)
#' }

#' @rdname CoxBaseEstimator
#' @export
CoxBaseEstimator <- function(object){
  UseMethod("CoxBaseEstimator") 
} 

#' @rdname CoxBaseEstimator
#' @method CoxBaseEstimator coxph
#' @export
CoxBaseEstimator.coxph <- function(object){
  return(object$method)
}

#' @rdname CoxBaseEstimator
#' @method CoxBaseEstimator phreg
#' @export
CoxBaseEstimator.phreg <- function(object){
  return("breslow")
}
# }}}

# {{{ CoxCenter
#' @title Extract the mean value of the covariates
#' @description Extract the mean value of the covariates
#' @rdname CoxCenter 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$Xcat <- paste0(d$X3,d$X4)
#' d$entry <- 0
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
#'
#' ##
#' library(mets)
#' mCox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' CoxCenter(mCox)
#' mCox <- cph(Surv(time, event) ~ X1+X2+Xcat, data = d)
#' CoxCenter(mCox)
#' }

#' @rdname CoxCenter
#' @export
CoxCenter <- function(object){
  UseMethod("CoxCenter") 
} 

#' @rdname CoxCenter
#' @method CoxCenter cph
#' @export
CoxCenter.cph <- function(object){
  return(setNames(object$means, object$mmcolnames))
}

#' @rdname CoxCenter
#' @method CoxCenter coxph
#' @export
CoxCenter.coxph <- function(object){
  return(setNames(object$means, names(coef(object))))
}

#' @rdname CoxCenter
#' @method CoxCenter phreg
#' @export
CoxCenter.phreg <- function(object){
  data <- model.matrix(object = eval(object$call$formula),  object$model.frame)[,names(coef(object)), drop = FALSE]
  return(apply(data,2,mean))
}
# }}}

# {{{ CoxDesign
#' @title Extract the design matrix used to train a Cox model
#' @description Extract the design matrix used to train a Cox model. Should contain the time of event, the type of event, 
#' the variable for the linear predictor, the strata variables and the date of entry (in case of delayed entry).
#' @rdname CoxDesign 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' @param center logical. Should the variable of the linear predictor be centered ?
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' set.seed(10)
#' d <- sampleData(2e1, outcome = "survival")
#' d$entry <- 0
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(entry, time, event) ~ X1+X2, data = d, x = TRUE, y = TRUE)
#' CoxDesign(mCox)
#' mCox <- coxph(Surv(time, event) ~ X1*X2+strata(X3), data = d, x = TRUE, y = TRUE)
#' CoxDesign(mCox)
#' mCox <- coxph(Surv(time, event) ~ X1+X2+strata(X3)+strata(X4), data = d, x = TRUE, y = TRUE)
#' CoxDesign(mCox, center = TRUE)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, x = TRUE, y = TRUE)
#' CoxDesign(mCox)
#' mCox <- cph(Surv(time, event) ~ X1*X2+strat(X3), data = d, x = TRUE, y = TRUE)
#' CoxDesign(mCox, center = TRUE)
#' mCox <- cph(Surv(time, event) ~ X1*X2+strat(X3)+strat(X4), data = d, x = TRUE, y = TRUE)
#' CoxDesign(mCox)
#' 
#' ##
#' library(mets) 
#' d$id <- 1:NROW(d)
#' 
#' mCox <- phreg(Surv(entry, time, event) ~ X1*X2, data = d)
#' CoxDesign(mCox)
#' mCox <- phreg(Surv(entry, time, event) ~ X1+X2+strata(X3)+cluster(id), data = d)
#' CoxDesign(mCox)
#' # mCox <- phreg(Surv(entry, time, event) ~ X1+X2+strata(X3)+strata(X4)+cluster(id), data = d)
#' }

#' @rdname CoxDesign
#' @export
CoxDesign <- function(object, center){
  UseMethod("CoxDesign") 
} 

#' @rdname CoxDesign
#' @method CoxDesign coxph
#' @export
CoxDesign.coxph <- function(object, center = FALSE){

    default.start <- 0
    
        
    if("x" %in% names(object) == FALSE){
        stop("invalid object \n",
             "set x=TRUE in the call to ",class(object)[1]," \n")
    }
        
    if("y" %in% names(object) == FALSE){
        stop("invalid object \n",
             "set y=TRUE in the call to ",class(object)[1]," \n")
    }
 
    ## set to null to be able to use cbind after 
    # (otherwise $x may have dimension <0,0> that cannot be bind with $y)
    if(NCOL(object[["x"]])==0){
        object[["x"]] <- NULL
    }else if(center){
        object[["x"]][] <- rowCenter_cpp(object[["x"]], center = CoxCenter(object))
    }

    if("strata" %in% names(object) == FALSE){
        object[["strata"]] <- NULL
    }
    if("start" %in% colnames(object$y) == FALSE){
        object[["y"]] <- cbind(start = default.start, 
                               stop = object[["y"]][,"time"],
                               status = object[["y"]][,"status"])
    }
    
    return(as.data.frame(cbind(object[["y"]],
                               object[["x"]],
                               strata=object[["strata"]])))
    
}

#' @rdname CoxDesign
#' @method CoxDesign phreg
#' @export
CoxDesign.phreg <- function(object, center = FALSE){
   M.outcome <- as.matrix(object$model.frame[,1])
   if("entry" %in% names(M.outcome) == FALSE){
     M.outcome <- cbind(entry = 0, M.outcome)
   }
     
   # normalize names
   name.default <- colnames(M.outcome)
   name.default<- gsub("entry","start",gsub("time","stop",name.default))
   colnames(M.outcome) <- name.default
   
   # get covariates
  M.X <- model.matrix(CoxFormula(object), data = object$model.frame)[,names(coef(object)),drop=FALSE]
  if(center){
    M.X <- rowCenter_cpp(M.X, center = CoxCenter(object))
  }
  
  if("strata" %in% names(object) == FALSE){
    object[["strata"]] <- NULL
  }
  
  return(as.data.frame(cbind(M.outcome, M.X, strata = object[["strata"]])))
}
# }}}

# {{{ CoxFormula
#' @title Extract the formula from a Cox model
#' @description Extract the formula from a Cox model
#' @rdname CoxFormula 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
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
#' f <- Surv(time, event) ~ X1+X2
#' mCox <- cph(f, data = d, y = TRUE)
#' CoxFormula(mCox)
#' }

#' @rdname CoxFormula
#' @export
CoxFormula <- function(object){
  UseMethod("CoxFormula") 
} 

#' @rdname CoxFormula
#' @method CoxFormula cph
#' @export
CoxFormula.cph <- function(object){
  return(object$sformula)
}

#' @rdname CoxFormula
#' @method CoxFormula coxph
#' @export
CoxFormula.coxph <- function(object){
  return(object$formula)
}

#' @rdname CoxFormula
#' @method CoxFormula phreg
#' @export
CoxFormula.phreg <- function(object){
  return(eval(object$call$formula))
}
# }}}

# {{{ CoxLP

#' @title Compute the linear predictor of a Cox model
#' @description Compute the linear predictor of a Cox model
#' @rdname CoxLP 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param center should the linear predictor be computed after centering the covariates
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @details In case of empty linear predictor returns a vector of 0 with the same length as the number of rows of the dataset
#'
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
#' d$id <- 1:NROW(d)
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
#' CoxLP(mCox, data = NULL, center = FALSE) 
#' CoxLP(mCox, data = NULL, center = TRUE) 
#' 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d)
#' CoxLP(mCoxS, data = d, center = FALSE) 
#' 
#' mCoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' CoxLP(mCoxS2, data = d, center = FALSE) 
#' 
#' ##
#' library(mets)
#' mCox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' CoxLP(mCox, data = NULL, center = FALSE) 
#' CoxLP(mCox, data = d, center = TRUE) 
#' 
#' mCoxS <- phreg(Surv(entry, time, event) ~ strata(X1)+X2+cluster(id), data = d)
#' CoxLP(mCoxS, data = d, center = FALSE) 
#' 
#' mCoxS2 <- phreg(Surv(entry, time, event) ~ X1*X2, data = d)
#' CoxLP(mCoxS2, data = d, center = FALSE) 
#
#' }

#' @rdname CoxLP
#' @export
CoxLP <- function(object, data, center){
  UseMethod("CoxLP") 
} 

#' @rdname CoxLP
#' @method CoxLP cph
#' @export
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
      
      if(center == FALSE){
        Xb <- Xb + sum(CoxCenter(object)*coef)
      }
      
    }else{ 
      Xb <- rep(0, NROW(data)) 
    } 
  }
  
  return(unname(Xb))
}

#' @rdname CoxLP
#' @method CoxLP coxph
#' @export
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
      
      object[["strata"]] <- NULL # solve a bug in survival:::predict.coxph when fitting the model with x = TRUE
      if(length(is.strata)>0){
        Xb <- rowSums(stats::predict(object, newdata = as.data.frame(data), 
                                     type = "terms"))
      }else{ 
        Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
      }
      
      if(center == FALSE){
        Xb <- Xb + sum(CoxCenter(object)*coef)
      }
      
    }else{ 
      Xb <- rep(0, NROW(data)) 
    } 
  }
  
  return(unname(Xb))
}

#' @rdname CoxLP
#' @method CoxLP phreg
#' @export
CoxLP.phreg <- function(object, data, center){
  coef <- stats::coef(object)
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    data <- object$model.frame
  }
  
  f.term <- delete.response(terms(eval(object$call$formula)))
  attr(f.term,"intercept") <- 0
  data <- model.matrix(object = f.term,  data)[,names(coef), drop = FALSE]
  
  if(n.varLP>0){
    Xb <- as.vector(as.matrix(data) %*% coef)
    
    if(center){
      Xb <- Xb - sum(CoxCenter(object)*coef)
    }
    
  }else{ 
    Xb <- rep(0, NROW(data)) 
  } 
  
  
  return(Xb)
}

# }}}

# {{{ CoxN
#' @title Extract the number of observations from a Cox model
#' @description Extract the number of observations from a Cox model
#' @rdname CoxN 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
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
#'
#' ##
#' library(mets)
#' mCox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' CoxN(mCox)
#' }

#' @rdname CoxN
#' @export
CoxN <- function(object){
  UseMethod("CoxN") 
} 

#' @rdname CoxN
#' @method CoxN cph
#' @export
CoxN.cph <- function(object){
  return(sum(object$n))
}

#' @rdname CoxN
#' @method CoxN coxph
#' @export
CoxN.coxph <- function(object){
  return(object$n)
}

#' @rdname CoxN
#' @method CoxN phreg
#' @export
CoxN.phreg <- function(object){
  return(NROW(object$model.frame))
}
# }}}

# {{{ CoxSpecialStrata
#' @title Special character for strata in Cox model
#' @description Return the special character used to indicate the strata variables of the Cox model
#' @name CoxSpecialStrata
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#'
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
#' 
#' ##
#' library(survival)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' CoxSpecialStrata(mCoxS)
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' CoxSpecialStrata(mCoxS)
#' 
#' ##
#' library(mets)
#' mCoxS <- phreg(Surv(entry, time, event) ~ strat(X1)+strat(X2), data = d)
#' CoxSpecialStrata(mCoxS)
#' }

#' @rdname CoxSpecialStrata
#' @export
CoxSpecialStrata <- function(object) UseMethod("CoxSpecialStrata")

#' @rdname CoxSpecialStrata
#' @method CoxSpecialStrata coxph
#' @export
CoxSpecialStrata.coxph <- function(object){
  return("strata")
}

#' @rdname CoxSpecialStrata
#' @method CoxSpecialStrata cph
#' @export
CoxSpecialStrata.cph <- function(object){
  return("strat")
}

#' @rdname CoxSpecialStrata
#' @method CoxSpecialStrata phreg
#' @export
CoxSpecialStrata.phreg <- function(object){
  return("strata")
}
# }}}

# {{{ CoxStrata
#' @title Define the strata for a new dataset
#' @description Define the strata in a dataset to match those of a stratified Cox model
#' @name CoxStrata
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
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
#' d$entry <- 0
#' d$id <- 1:NROW(d)
#' 
#' ##
#' library(survival)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d, x = TRUE, y = TRUE)
#' resInfo <- CoxVariableName(mCoxS)
#' Ostrata <- CoxStrata(mCoxS, stratavars = resInfo$stratavars)
#' CoxStrata(mCoxS, data = d, sterms = resInfo$sterms, stratavars = resInfo$stratavars, 
#'           levels = levels(Ostrata), stratalevels = resInfo$stratalevels)
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' resInfo <- CoxVariableName(mCoxS)
#' Ostrata <- CoxStrata(mCoxS, stratavars = resInfo$stratavars)
#' CoxStrata(mCoxS, data = d, sterms = resInfo$sterms, stratavars = resInfo$stratavars, 
#'           levels = levels(Ostrata), stratalevels = resInfo$stratalevels)
#'           
#' ##
#' library(mets)
#' mCoxS <- phreg(Surv(entry, time, event) ~ strata(X1)+X2+cluster(id), data = d)
#' resInfo <- CoxVariableName(mCoxS)
#' Ostrata <- CoxStrata(mCoxS, stratavars = resInfo$stratavars)
#' CoxStrata(mCoxS, data = d, sterms = resInfo$sterms, stratavars = resInfo$stratavars, 
#'           levels = levels(Ostrata), stratalevels = resInfo$stratalevels)
#' }

#' @rdname CoxStrata
#' @export
CoxStrata <- function(object, data, sterms, stratavars, levels, stratalevels) UseMethod("CoxStrata")

#' @rdname CoxStrata
#' @method CoxStrata coxph
#' @export
CoxStrata.cph <- function(object, data, sterms, stratavars, levels, stratalevels){
  
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
#' @method CoxStrata coxph
#' @export
CoxStrata.coxph <- function(object, data, sterms, stratavars, levels, stratalevels){
  
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

#' @rdname CoxStrata
#' @method CoxStrata phreg
# '@export
CoxStrata.phreg <- CoxStrata.coxph
# }}}

# {{{ CoxVarCov
#' @title Extract the variance covariance matrix of the beta from a Cox model
#' @description Extract the variance covariance matrix of the beta from a Cox model
#' @rdname CoxVarCov 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @details Should return \code{NULL} if the Cox model has no covariate. 
#' The rows and columns of the variance covariance matrix must be named with the names used in the design matrix.
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
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
#' ##
#' library(mets)
#' mCox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' CoxVarCov(mCox)
#' mCox <- cph(Surv(time, event) ~ 1, data = d, y = TRUE)
#' CoxVarCov(mCox)
#' }

#' @rdname CoxVarCov
#' @export
CoxVarCov <- function(object){
  UseMethod("CoxVarCov") 
} 

#' @rdname CoxVarCov
#' @method CoxVarCov cph
#' @export
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
#' @export
CoxVarCov.coxph <- function(object){
  
  Sigma <- object$var
  if(!is.null(Sigma)){
    coefName <- names(coef(object))
    colnames(Sigma) <- coefName
    rownames(Sigma) <- coefName
  }
  return(Sigma)
}

#' @rdname CoxVarCov
#' @method CoxVarCov phreg
#' @export
CoxVarCov.phreg <- function(object){
  
  Sigma <- -solve(object$hessian)
  if(!is.null(Sigma)){
    coefName <- names(coef(object))
    colnames(Sigma) <- coefName
    rownames(Sigma) <- coefName
  }
  return(Sigma)
}
# }}}


#### Auxiliary function #### 

# {{{ SurvResponseVar
#' @title Extract the time and event variable from a Cox model
#' @description Extract the time and event variable from a Cox model
#' @rdname SurvResponseVar
#' @param formula a formula
#'
#' @author Brice Ozenne broz@@sund.ku.dk
#'
#' @examples
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' SurvResponseVar(mCox$formula)
#' 
#' mCox <- coxph(Surv(entry, time, event) ~ X1+X2, data = d)
#' SurvResponseVar(mCox$formula)
#' 
#' mCox <- coxph(Surv(event, time = entry, time2 = time) ~ X1+X2, data = d)
#' SurvResponseVar(mCox$formula)
#'
#' mCox <- coxph(Surv(time, event) ~ 1, data = d)
#' SurvResponseVar(mCox$formula)
#'
#' }
SurvResponseVar <- function(formula){
  
  formula <- formula[[2]]
  length.f <- length(formula)
  
  if(length.f==4){
    all.names <- c("time","time2","event")
    all.positions <- 2:4
  }else{ # length.f==3
    all.names <- c("time","event")
    all.positions <- 2:3
  }  
  
  ## find named argument
  all.names_copy <- all.names
  for(iNames in all.names_copy){
    if(iNames %in% names(formula)){
      all.names <- setdiff(all.names, iNames)
      all.positions <- setdiff(all.positions, which(iNames %in% names(formula)))
    }
  }
  
  ## set names to unnamed arguments
  nRemaining <- length(all.names)
  if(nRemaining>0){
    for(iNames in 1:nRemaining){
      names(formula)[all.positions[iNames]] <- all.names[iNames]
    } 
  }
  
  ## export
  return(list(entry = if(length.f==4){deparse(formula$time)}else{NULL},
              time = if(length.f==4){deparse(formula$time2)}else{deparse(formula$time)},
              status = deparse(formula$event)
  )
  )
  
}
# }}}

# {{{ extractStrata
#' @title Extract the information about the strata
#' @description Extract the information about the strata stored in a Cox model
#' @rdname extractStrata
#' 
#' @param xterms ...
#' @param xlevels ...
#' @param special the special character indicating the strata variables
#'
#' @author Brice Ozenne broz@@sund.ku.dk and Thomas A. Gerds tag@@biostat.ku.dk
#'
#' @examples
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
#' 
#' ##
#' library(survival)
#' mCox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' extractStrata(delete.response(mCox$terms), special = CoxSpecialStrata(mCox))
#' 
#' mCox <- coxph(Surv(event, time = entry, time2 = time) ~ strata(X1)+strata(X2), data = d)
#' extractStrata(delete.response(mCox$terms),
#'               xlevels =  mCox$xlevels, special = CoxSpecialStrata(mCox))
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d)
#' extractStrata(delete.response(mCox$terms), special = CoxSpecialStrata(mCox))
#' 
#' mCox <- cph(Surv(event, time = entry, time2 = time) ~ strat(X1)+strat(X2), data = d)
#' extractStrata(delete.response(mCox$terms),
#'               xlevels =  mCox$xlevels, special = CoxSpecialStrata(mCox))
#'
#' }
extractStrata <- function(xterms, xlevels = NULL, special){
  
  # renamed variables (e.g. strata(X1))
  xvars <- attr(xterms,"term.labels")
  strataspecials <- attr(xterms,"specials")[[special]]
  stratavars <- xvars[strataspecials]
  
  # original variables (e.g. X1)
  allVars.X <- all.vars(xterms)
  stratavars.original <- allVars.X[strataspecials]
  
  is.strata <- length(strataspecials)>0
  
  if(is.strata){
    if (length(xvars)>length(strataspecials)){
      sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
    } else {
      sterms <- xterms
    }
    
  }else{
    sterms <- NULL
  }
  
  # levels of the strata variables [only useful for coxph]
  if(!is.null(xlevels)){
    stratalevels <- xlevels[stratavars]
  }else{
    stratalevels <- NULL
  }
  
  return(list(stratavars = stratavars,
              stratavars.original = stratavars.original,
              strataspecials = strataspecials,
              stratalevels = stratalevels,
              is.strata = is.strata,
              sterms = sterms))
}
# }}}

# {{{ model.matrix.phreg
#' @title Extract design matrix for phreg objects
#' @description Extract design matrix for phreg objects
#' @param object a phreg object.
#' @param data a dataset.
#' 
#' @details mainly a copy paste of the begining of the \code{phreg} function.
#' 
#' @examples 
#' \dontrun{
#' library(mets)
#' 
#' n <- 10;
#' d <- mets:::simCox(n); d$id <- seq(nrow(d)); d$group <- factor(rbinom(nrow(d),1,0.5))
#' m1 <- phreg(Surv(entry, time,status)~X1+X2+cluster(id)+strata(group),data=d)
#' riskRegression:::model.matrix(m1, d) 
#' }
#' @method model.matrix phreg
model.matrix.phreg <- function(object, data){
  special <- c("strata", "cluster")
  Terms <- terms(eval(object$call$formula), special, data = data)

  ## remove specials
  if (!is.null(attributes(Terms)$specials$cluster)) {
    ts <- survival::untangle.specials(Terms, "cluster")
    Terms <- Terms[-ts$terms]
  }
  if (!is.null(stratapos <- attributes(Terms)$specials$strata)) {
    ts <- survival::untangle.specials(Terms, "strata")
    Terms <- Terms[-ts$terms]
  }
  attr(Terms,"intercept") <- 0
    
  return(model.matrix(Terms, data))
}

# }}}
