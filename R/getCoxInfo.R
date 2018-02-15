#### functions ####

# {{{ coxVariableName
#' @title Extract variable names from a model
#' @description Extract the name of the variables belonging to the linear predictor or used to form the strata
#' @rdname coxVariableName 
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d, x = TRUE, y = TRUE)
#' coxVariableName(mcox)
#' mcoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d, x = TRUE, y = TRUE)
#' coxVariableName(mcoxS)
#' mcoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d, x = TRUE, y = TRUE)
#' coxVariableName(mcoxS2)
#' mcoxI <- coxph(Surv(time, event) ~ X1*X2, data = d, x = TRUE, y = TRUE)
#' coxVariableName(mcoxI)
#' 
#' ##
#' library(rms)
#' mcox <- cph(Surv(time, event) ~ X1+X2, data = d, x = TRUE, y = TRUE)
#' coxVariableName(mcox) 
#' mcoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d, x = TRUE, y = TRUE)
#' coxVariableName(mcoxS) 
#' mcoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, x = TRUE, y = TRUE)
#' coxVariableName(mcoxS2) 
#' mcoxI <- cph(Surv(time, event) ~ X1*X2, data = d, x = TRUE, y = TRUE)
#' coxVariableName(mcoxI)
#' 
#' ##
#' library(mets)
#' mcox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' coxVariableName(mcox) 
#' mcoxS <- phreg(Surv(entry, time, event) ~ strata(X1)+cluster(id)+X2, data = d)
#' coxVariableName(mcoxS) 
#' mcoxI <- phreg(Surv(entry, time, event) ~ X1*X2, data = d)
#' coxVariableName(mcoxI)
#' }

#' @rdname coxVariableName
#' @export
coxVariableName <- function(object){
  
  f <- coxFormula(object)
  special.object <- coxSpecialStrata(object)
  
  ## response
  ls.SurvVar <- SurvResponseVar(f)
  
  ## strata
  xterms <- delete.response(terms(f,
                                  special = special.object,
                                  data = coxDesign(object)))
  ls.StrataInfo <- extractStrata(xterms,
                                 special = coxSpecialStrata(object))
  
  ## regressor
  lpvars.original <- setdiff(all.vars(f),c(unlist(ls.SurvVar), ls.StrataInfo$strata.vars.original))
  
  ## export
  return(c(list(entry = ls.SurvVar$entry),
           list(time = ls.SurvVar$time),
           list(status = ls.SurvVar$status),
           ls.StrataInfo,
           list(lpvars = names(coef(object))),
           list(lpvars.original = lpvars.original)
  ))
} 
# }}}

coxCovars <- function(object){
  ttobj <- stats::terms(object)
  ## colnames(attr(ttobj,"factors"))
  all.vars(attr(delete.response(ttobj),"variables"))
}

#### methods #####

# {{{ coxBaseEstimator
#' @title Extract the type of estimator for the baseline hazard
#' @description Extract the type of estimator for the baseline hazard
#' @rdname coxBaseEstimator 
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' coxBaseEstimator(mcox)
#'
#' ##
#' library(rms)
#' mcox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' coxBaseEstimator(mcox)
#'
#' ##
#' library(mets)
#' mcox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' coxBaseEstimator(mcox)
#' }

#' @rdname coxBaseEstimator
#' @export
coxBaseEstimator <- function(object){
  UseMethod("coxBaseEstimator") 
} 

#' @rdname coxBaseEstimator
#' @method coxBaseEstimator coxph
#' @export
coxBaseEstimator.coxph <- function(object){
  return(object$method)
}

#' @rdname coxBaseEstimator
#' @method coxBaseEstimator phreg
#' @export
coxBaseEstimator.phreg <- function(object){
  return("breslow")
}
# }}}

# {{{ coxCenter
#' @title Extract the mean value of the covariates
#' @description Extract the mean value of the covariates
#' @rdname coxCenter 
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' coxCenter(mcox)
#' mcox <- coxph(Surv(time, event) ~ X1+X2+Xcat, data = d)
#' coxCenter(mcox)
#' 
#' ##
#' library(rms)
#' mcox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' coxCenter(mcox)
#' mcox <- cph(Surv(time, event) ~ X1+X2+Xcat, data = d)
#' coxCenter(mcox)
#'
#' ##
#' library(mets)
#' mcox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' coxCenter(mcox)
#' mcox <- cph(Surv(time, event) ~ X1+X2+Xcat, data = d)
#' coxCenter(mcox)
#' }

#' @rdname coxCenter
#' @export
coxCenter <- function(object){
  UseMethod("coxCenter") 
} 

#' @rdname coxCenter
#' @method coxCenter cph
#' @export
coxCenter.cph <- function(object){
  return(setNames(object$means, object$mmcolnames))
}

#' @rdname coxCenter
#' @method coxCenter coxph
#' @export
coxCenter.coxph <- function(object){
  return(setNames(object$means, names(coef(object))))
}

#' @rdname coxCenter
#' @method coxCenter phreg
#' @export
coxCenter.phreg <- function(object){
  data <- model.matrix(object = eval(object$call$formula),  object$model.frame)[,names(coef(object)), drop = FALSE]
  return(apply(data,2,mean))
}
# }}}

# {{{ coxDesign
#' @title Extract the design matrix used to train a Cox model
#' @description Extract the design matrix used to train a Cox model. Should contain the time of event, the type of event, 
#' the variable for the linear predictor, the strata variables and the date of entry (in case of delayed entry).
#' @rdname coxDesign 
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
#' mcox <- coxph(Surv(entry, time, event) ~ X1+X2, data = d, x = TRUE, y = TRUE)
#' coxDesign(mcox)
#' mcox <- coxph(Surv(time, event) ~ X1*X2+strata(X3), data = d, x = TRUE, y = TRUE)
#' coxDesign(mcox)
#' mcox <- coxph(Surv(time, event) ~ X1+X2+strata(X3)+strata(X4), data = d, x = TRUE, y = TRUE)
#' coxDesign(mcox, center = TRUE)
#' 
#' ##
#' library(rms)
#' mcox <- cph(Surv(time, event) ~ X1+X2, data = d, x = TRUE, y = TRUE)
#' coxDesign(mcox)
#' mcox <- cph(Surv(time, event) ~ X1*X2+strat(X3), data = d, x = TRUE, y = TRUE)
#' coxDesign(mcox, center = TRUE)
#' mcox <- cph(Surv(time, event) ~ X1*X2+strat(X3)+strat(X4), data = d, x = TRUE, y = TRUE)
#' coxDesign(mcox)
#' 
#' ##
#' library(mets) 
#' d$id <- 1:NROW(d)
#' 
#' mcox <- phreg(Surv(entry, time, event) ~ X1*X2, data = d)
#' coxDesign(mcox)
#' mcox <- phreg(Surv(entry, time, event) ~ X1+X2+strata(X3)+cluster(id), data = d)
#' coxDesign(mcox)
#' # mcox <- phreg(Surv(entry, time, event) ~ X1+X2+strata(X3)+strata(X4)+cluster(id), data = d)
#' }

#' @rdname coxDesign
#' @export
coxDesign <- function(object, center){
  UseMethod("coxDesign") 
} 

#' @rdname coxDesign
#' @method coxDesign coxph
#' @export
coxDesign.coxph <- function(object, center = FALSE){
  
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
    object[["x"]][] <- rowCenter_cpp(object[["x"]], center = coxCenter(object))
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

#' @rdname coxDesign
#' @method coxDesign phreg
#' @export
coxDesign.phreg <- function(object, center = FALSE){
  M.outcome <- as.matrix(object$model.frame[,1])
  if("entry" %in% names(M.outcome) == FALSE){
    M.outcome <- cbind(entry = 0, M.outcome)
  }
  
  # normalize names
  name.default <- colnames(M.outcome)
  name.default<- gsub("entry","start",gsub("time","stop",name.default))
  colnames(M.outcome) <- name.default
  
  # get covariates
  M.X <- model.matrix(coxFormula(object), data = object$model.frame)[,names(coef(object)),drop=FALSE]
  
  if(center){
    M.X <- rowCenter_cpp(M.X, center = coxCenter(object))
  }
  
  if("strata" %in% names(object) == FALSE){
    object[["strata"]] <- NULL
  }
  
  return(as.data.frame(cbind(M.outcome, M.X, strata = object[["strata"]])))
}
# }}}

# {{{ coxFormula
#' @title Extract the formula from a Cox model
#' @description Extract the formula from a Cox model
#' @rdname coxFormula 
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' coxFormula(mcox)
#' 
#' ##
#' library(rms)
#' f <- Surv(time, event) ~ X1+X2
#' mcox <- cph(f, data = d, y = TRUE)
#' coxFormula(mcox)
#' }

#' @rdname coxFormula
#' @export
coxFormula <- function(object){
  UseMethod("coxFormula") 
} 

#' @rdname coxFormula
#' @method coxFormula cph
#' @export
coxFormula.cph <- function(object){
  return(object$sformula)
}

#' @rdname coxFormula
#' @method coxFormula coxph
#' @export
coxFormula.coxph <- function(object){
  return(object$formula)
}

#' @rdname coxFormula
#' @method coxFormula phreg
#' @export
coxFormula.phreg <- function(object){
  return(eval(object$call$formula))
}
# }}}

# {{{ coxLP

#' @title Compute the linear predictor of a Cox model
#' @description Compute the linear predictor of a Cox model
#' @rdname coxLP 
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' coxLP(mcox, data = d, center = FALSE)
#' coxLP(mcox, data = d, center = TRUE)  
#' 
#' mcoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d)
#' coxLP(mcoxS, data = d, center = FALSE) 
#' 
#' mcoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' coxLP(mcoxS2, data = d, center = FALSE) 
#' 
#' ##
#' library(rms)
#' mcox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' coxLP(mcox, data = d, center = FALSE) 
#' coxLP(mcox, data = d, center = TRUE) 
#' coxLP(mcox, data = NULL, center = FALSE) 
#' coxLP(mcox, data = NULL, center = TRUE) 
#' 
#' mcoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d)
#' coxLP(mcoxS, data = d, center = FALSE) 
#' 
#' mcoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' coxLP(mcoxS2, data = d, center = FALSE) 
#' 
#' ##
#' library(mets)
#' mcox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' coxLP(mcox, data = NULL, center = FALSE) 
#' coxLP(mcox, data = d, center = TRUE) 
#' 
#' mcoxS <- phreg(Surv(entry, time, event) ~ strata(X1)+X2+cluster(id), data = d)
#' coxLP(mcoxS, data = d, center = FALSE) 
#' 
#' mcoxS2 <- phreg(Surv(entry, time, event) ~ X1*X2, data = d)
#' coxLP(mcoxS2, data = d, center = FALSE) 
#
#' }

#' @rdname coxLP
#' @export
coxLP <- function(object, data, center){
  UseMethod("coxLP") 
} 

#' @rdname coxLP
#' @method coxLP cph
#' @export
coxLP.cph <- function(object, data, center){
  
  coef <- stats::coef(object)
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    
    Xb <- object$linear.predictors
    
    if(center == FALSE && n.varLP != 0){
      Xb <- Xb + sum(coxCenter(object)*coef)
    }
    
  }else{ ## new dataset
    
    if(n.varLP>0){
      
      Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
      
      if(center == FALSE){
        Xb <- Xb + sum(coxCenter(object)*coef)
      }
      
    }else{ 
      Xb <- rep(0, NROW(data)) 
    } 
  }
  
  return(unname(Xb))
}

#' @rdname coxLP
#' @method coxLP coxph
#' @export
coxLP.coxph <- function(object, data, center){
  
  coef <- stats::coef(object)
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    
    Xb <- object$linear.predictors
    
    if(center == FALSE && n.varLP != 0){
      Xb <- Xb + sum(coxCenter(object)*coef)
    }
    
  }else{ ## new dataset
    if(n.varLP>0){
      is.strata <- attr(object$terms, "special")$strata
      
      
      if(length(is.strata)>0){
        object.strata <- object[["strata"]]
        object[["strata"]] <- NULL # solve a bug in survival:::predict.coxph when fitting the model with x = TRUE
        
        Xb <- try(rowSums(stats::predict(object, newdata = as.data.frame(data), 
                                         type = "terms")), silent = TRUE)
        if("try-error" %in% class(Xb)){ ## Fix an error when the dataset used to fit the object is removed from the global environment
          ## survival:::predict.coxph search for it and read (at least) the status variable
          txt <- paste0("survival::predict.coxph returns the following error:\n",
                        as.character(Xb),
                        "It seems that the dataset used to fit the model is no more compatible with the model,\n",
                        "probably because it has been modified afterwards.\n",
                        "coxLP.coxph will try to reconstruct the original dataset and continue the execution.\n")
          warning(txt)
          ## So avoid an error, the following code re-create the original dataset
          object[["strata"]] <- object.strata
          object$call$data <- reconstructData(object)
          object[["strata"]] <- NULL
          
          Xb <- rowSums(stats::predict(object, newdata = as.data.frame(data), 
                                       type = "terms"))
        }
      }else{ 
        Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
      }
      
      if(center == FALSE){
        Xb <- Xb + sum(coxCenter(object)*coef)
      }
      
    }else{ 
      Xb <- rep(0, NROW(data)) 
    } 
  }
  
  return(unname(Xb))
}

#' @rdname coxLP
#' @method coxLP phreg
#' @export
coxLP.phreg <- function(object, data, center){
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
      Xb <- Xb - sum(coxCenter(object)*coef)
    }
    
  }else{ 
    Xb <- rep(0, NROW(data)) 
  } 
  
  
  return(Xb)
}

# }}}

# {{{ coxN
#' @title Extract the number of observations from a Cox model
#' @description Extract the number of observations from a Cox model
#' @rdname coxN 
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' coxN(mcox)
#' 
#' ##
#' library(rms)
#' mcox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' coxN(mcox)
#'
#' ##
#' library(mets)
#' mcox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' coxN(mcox)
#' }

#' @rdname coxN
#' @export
coxN <- function(object){
  UseMethod("coxN") 
} 

#' @rdname coxN
#' @method coxN cph
#' @export
coxN.cph <- function(object){
  return(sum(object$n))
}

#' @rdname coxN
#' @method coxN coxph
#' @export
coxN.coxph <- function(object){
  return(object$n)
}

#' @rdname coxN
#' @method coxN phreg
#' @export
coxN.phreg <- function(object){
  return(NROW(object$model.frame))
}
# }}}

# {{{ coxSpecialStrata
#' @title Special character for strata in Cox model
#' @description Return the special character used to indicate the strata variables of the Cox model
#' @name coxSpecialStrata
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
#' mcoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' coxSpecialStrata(mcoxS)
#' 
#' ##
#' library(rms)
#' mcoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' coxSpecialStrata(mcoxS)
#' 
#' ##
#' library(mets)
#' mcoxS <- phreg(Surv(entry, time, event) ~ strat(X1)+strat(X2), data = d)
#' coxSpecialStrata(mcoxS)
#' }

#' @rdname coxSpecialStrata
#' @export
coxSpecialStrata <- function(object) UseMethod("coxSpecialStrata")

#' @rdname coxSpecialStrata
#' @method coxSpecialStrata coxph
#' @export
coxSpecialStrata.coxph <- function(object){
  return("strata")
}

#' @rdname coxSpecialStrata
#' @method coxSpecialStrata cph
#' @export
coxSpecialStrata.cph <- function(object){
  return("strat")
}

#' @rdname coxSpecialStrata
#' @method coxSpecialStrata phreg
#' @export
coxSpecialStrata.phreg <- function(object){
  return("strata")
}
# }}}

# {{{ coxStrata
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
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' d$entry <- 0
#' d$id <- 1:NROW(d)
#' 
#' ##
#' library(survival)
#' mcoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d, x = TRUE, y = TRUE)
#' resInfo <- coxVariableName(mcoxS)
#' Ostrata <- coxStrata(mcoxS, strata.vars = resInfo$strata.vars)
#' coxStrata(mcoxS, data = d, sterms = resInfo$sterms, strata.vars = resInfo$strata.vars, 
#'           levels = levels(Ostrata), strata.levels = resInfo$strata.levels)
#' 
#' ##
#' library(rms)
#' mcoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' resInfo <- coxVariableName(mcoxS)
#' Ostrata <- coxStrata(mcoxS, strata.vars = resInfo$strata.vars)
#' coxStrata(mcoxS, data = d, sterms = resInfo$sterms, strata.vars = resInfo$strata.vars, 
#'           levels = levels(Ostrata), strata.levels = resInfo$strata.levels)
#'           
#' ##
#' library(mets)
#' mcoxS <- phreg(Surv(entry, time, event) ~ strata(X1)+X2+cluster(id), data = d)
#' resInfo <- coxVariableName(mcoxS)
#' Ostrata <- coxStrata(mcoxS, strata.vars = resInfo$strata.vars)
#' coxStrata(mcoxS, data = d, sterms = resInfo$sterms, strata.vars = resInfo$strata.vars, 
#'           levels = levels(Ostrata), strata.levels = resInfo$strata.levels)
#' }

#' @rdname coxStrata
#' @export
coxStrata <- function(object, data, sterms, strata.vars, levels, strata.levels) UseMethod("coxStrata")

#' @rdname coxStrata
#' @method coxStrata coxph
#' @export
coxStrata.cph <- function(object, data, sterms, strata.vars, levels, strata.levels){
  
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
      if (any(unique(strata) %in% levels == FALSE)){
        stop("unknown strata: ",paste(unique(strata[strata %in% levels == FALSE]), collapse = " | "),"\n")
      }
      strata <- factor(strata, levels = levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert strata to numeric
    }
    
  }
  
  return(strata)
}


#' @rdname coxStrata
#' @method coxStrata coxph
#' @export
coxStrata.coxph <- function(object, data, sterms, strata.vars, levels, strata.levels){
  
  if(length(strata.vars)==0){ ## no strata variables
    
    n <- if(is.null(data)) coxN(object) else NROW(data)
    strata <- as.factor(rep("1", n))
    
  }else{  ## strata variables
    
    if(is.null(data)){ ## training dataset
      strata <- object$strata
    }else { ## new dataset
      strata <- prodlim::model.design(sterms,data=data,xlev=strata.levels,specialsFactor=TRUE)$strata[[1]]
      if (any(unique(strata) %in% levels == FALSE)){
        stop("unknown strata: ",paste(unique(strata[strata %in% levels == FALSE]), collapse = " | "),"\n")
      }
      strata <- factor(strata, levels = levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert strata to numeric
    }
    
  }
  return(strata)
}

#' @rdname coxStrata
#' @method coxStrata phreg
# '@export
coxStrata.phreg <- coxStrata.coxph
# }}}

# {{{ coxVarCov
#' @title Extract the variance covariance matrix of the beta from a Cox model
#' @description Extract the variance covariance matrix of the beta from a Cox model
#' @rdname coxVarCov 
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' coxVarCov(mcox)
#' mcox <- coxph(Surv(time, event) ~ 1, data = d, y = TRUE)
#' coxVarCov(mcox)
#' 
#' ##
#' library(rms)
#' mcox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' coxVarCov(mcox)
#' mcox <- cph(Surv(time, event) ~ 1, data = d, y = TRUE)
#' coxVarCov(mcox)
#' 
#' ##
#' library(mets)
#' mcox <- phreg(Surv(entry, time, event) ~ X1+X2, data = d)
#' coxVarCov(mcox)
#' mcox <- cph(Surv(time, event) ~ 1, data = d, y = TRUE)
#' coxVarCov(mcox)
#' }

#' @rdname coxVarCov
#' @export
coxVarCov <- function(object){
  UseMethod("coxVarCov") 
} 

#' @rdname coxVarCov
#' @method coxVarCov cph
#' @export
coxVarCov.cph <- function(object){
  
  Sigma <- object$var
  if(!is.null(Sigma)){
    colnames(Sigma) <- object$mmcolnames
    rownames(Sigma) <- object$mmcolnames
  }
  
  return(Sigma)
}

#' @rdname coxVarCov
#' @method coxVarCov coxph
#' @export
coxVarCov.coxph <- function(object){
  
  Sigma <- object$var
  if(!is.null(Sigma)){
    coefName <- names(coef(object))
    colnames(Sigma) <- coefName
    rownames(Sigma) <- coefName
  }
  return(Sigma)
}

#' @rdname coxVarCov
#' @method coxVarCov phreg
#' @export
coxVarCov.phreg <- function(object){
  
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' SurvResponseVar(mcox$formula)
#' 
#' mcox <- coxph(Surv(entry, time, event) ~ X1+X2, data = d)
#' SurvResponseVar(mcox$formula)
#' 
#' mcox <- coxph(Surv(event, time = entry, time2 = time) ~ X1+X2, data = d)
#' SurvResponseVar(mcox$formula)
#'
#' mcox <- coxph(Surv(time, event) ~ 1, data = d)
#' SurvResponseVar(mcox$formula)
#' 
#' mcox <- coxph(Surv(time, event>0) ~ 1, data = d)
#' SurvResponseVar(mcox$formula)
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
  
  ## extract variable names
  name.vars <- all.vars(formula)
  names(name.vars) <- names(formula)[-1]
  
  ## export
  return(list(entry = if(length.f==4){unname(name.vars["time"])}else{NULL},
              time = if(length.f==4){unname(name.vars["time2"])}else{unname(name.vars["time"])},
              status = unname(name.vars["event"])
  ))
  
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
#' mcox <- coxph(Surv(time, event) ~ X1+X2, data = d)
#' extractStrata(delete.response(mcox$terms), special = coxSpecialStrata(mcox))
#' 
#' mcox <- coxph(Surv(event, time = entry, time2 = time) ~ strata(X1)+strata(X2), data = d)
#' extractStrata(delete.response(mcox$terms),
#'               xlevels =  mcox$xlevels, special = coxSpecialStrata(mcox))
#' 
#' ##
#' library(rms)
#' mcox <- cph(Surv(time, event) ~ X1+X2, data = d)
#' extractStrata(delete.response(mcox$terms), special = coxSpecialStrata(mcox))
#' 
#' mcox <- cph(Surv(event, time = entry, time2 = time) ~ strat(X1)+strat(X2), data = d)
#' extractStrata(delete.response(mcox$terms),
#'               xlevels =  mcox$xlevels, special = coxSpecialStrata(mcox))
#'
#' }
extractStrata <- function(xterms, xlevels = NULL, special){
  
  # renamed variables (e.g. strata(X1))
  xvars <- attr(xterms,"term.labels")
  strataspecials <- attr(xterms,"specials")[[special]]
  strata.vars <- xvars[strataspecials]
  
  # original variables (e.g. X1)
  allVars.X <- all.vars(xterms)
  strata.vars.original <- allVars.X[strataspecials]
  
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
    strata.levels <- xlevels[strata.vars]
  }else{
    strata.levels <- NULL
  }
  
  return(list(strata.vars = strata.vars,
              strata.vars.original = strata.vars.original,
              strataspecials = strataspecials,
              strata.levels = strata.levels,
              is.strata = is.strata,
              sterms = sterms))
}
# }}}

# {{{ splitStrataVar
#' @title Reconstruct each of the strata variables
#' @description Reconstruct each of the strata variables from the strata variable stored in the coxph object.
#' @rdname splitStrataVar
#' 
#' @param object a coxph object.
#'
#' @author Brice Ozenne broz@@sund.ku.dk and Thomas A. Gerds tag@@biostat.ku.dk
#'
#' @examples
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' library(survival)
#' mcox <- coxph(Surv(time, event) ~ strata(X1)+strata(X2)+X3, data = d, x = TRUE, y = TRUE)
#' splitStrataVar(mcox)
splitStrataVar <- function(object){
  
  grid.levels <- expand.grid(object$xlevels)
  vec.levels <- apply(grid.levels,1,paste,collapse=", ")
  
  ## find original name of the strata variables
  xterms <- delete.response(object$terms)
  name.strataVar.original <- all.vars(xterms)[attr(xterms,"specials")[["strata"]]]
  name.strataVar <- attr(xterms,"term.labels")[attr(xterms,"specials")[["strata"]]]
  n.strataVar <- length(name.strataVar)
  
  ## identify value relative to each strata variable
  df.strata <- NULL
  for(iStrataVar in 1:n.strataVar){ # iStrataVar <- 1
    value2level <- setNames(grid.levels[,name.strataVar[iStrataVar]], vec.levels)
    if(iStrataVar==1){
      df.strata <- data.frame(unname(value2level[object$strata]))
    }else{
      df.strata <- cbind(df.strata,
                         unname(value2level[object$strata]))
    }
  }
  
  ## export
  names(df.strata) <- name.strataVar.original
  return(df.strata)
}
# }}}

# {{{ 
#' @title Reconstruct the original dataset
#' @description Reconstruct the original dataset from the elements stored in the coxph object
#' @rdname splitStrataVar
#' 
#' @param object a coxph object.
#'
#' @author Brice Ozenne broz@@sund.ku.dk and Thomas A. Gerds tag@@biostat.ku.dk
#'
#' @examples
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' library(survival)
#' mcox <- coxph(Surv(time, event) ~ strata(X1)+strata(X2)+X3, data = d, x = TRUE, y = TRUE)
#' reconstructData(mcox)
reconstructData <- function(object){
  
  
  ## combine response variable, regressors and strata variable
  newdata <- as.data.frame(cbind(coxDesign(object),splitStrataVar(object)))
  
  ## set response variable to their original names
  infoVar <- coxVariableName(object)
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
  Terms <- terms(coxFormula(object), special, data = data)
  
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
