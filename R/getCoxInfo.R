#### functions ####

                                        # {{{ coxVariableName
## * coxVariableName
#' @title Extract variable names from a model
#' @description Extract the name of the variables belonging to the linear predictor or used to form the strata
#' @name coxVariableName
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param model.frame [data.frame] dataset containing all the relevant variables (entry, time to event, type of event, variables in the linear predictor, strata).
#' Output from \code{coxModelFrame}.
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxVariableName
#' @export
coxVariableName <- function(object, model.frame){

    ## ** get formula
    ff <- coxFormula(object)

    ## ** get special operators
    special.object <- coxSpecial(object)

    ## ** get names of the start/stop/status variables
    out <- SurvResponseVar(ff)

    ## ** identify special terms
    xterms <- delete.response(terms(ff,
                                    special = unlist(special.object),
                                    data = model.frame))    
    ls.specials <- attr(xterms,"specials")
    if(any(attr(xterms,"order")>1) && length(ls.specials$strata)>0){
        ## the presence of an interaction makes attr(xterms,"specials")$strata consider the 'main' effects which may not automatically be there
        ## e.g. coxph(Surv(time,event)~X1*X3+strata(X5),data=train, y=TRUE, x = TRUE) lead to design matrix [X1,X3,strata(X5),X1:X3] and attr(xterms,"specials")$strata being 3 (as expected)
        ## e.g. coxph(Surv(time,event)~X1:X3+strata(X5),data=train, y=TRUE, x = TRUE) lead to design matrix [strata(X5),X1:X3] and attr(xterms,"specials")$strata being 3 (incorrect)
        ## This is a fix
        term.strata <- names(which(colSums(attr(xterms,"factors")[attr(xterms,"specials")$strata,,drop=FALSE])>0))
        ls.specials$strata <- intersect(which(attr(xterms,"order")==1), ## strata terms are always main terms (handle X2*strata(X1)--> strata(X1), X2:strata(X1))
                                        which(attr(xterms,"term.labels") %in% term.strata))
    }
    n.specials <- length(unlist(ls.specials))
    
    n.xterms <- length(attr(xterms,"term.labels"))

    ## ** linear predictor
    if(n.xterms>n.specials){
        out$lpvars <- names(coef(object)) ##attr(xterms.lp, "term.labels") not ok for categorical variables with coxph
        if(n.specials>0){
            out$lp.sterms <- stats::drop.terms(xterms,unlist(ls.specials))
        }else{
            out$lp.sterms <- xterms
        }
        out$lpvars.original <- all.vars(out$lp.sterms)
        
    }else{
        out["lpvars"] <- list(NULL)
        out["lpvars.original"] <- list(NULL)
        out["lp.sterms"] <- list(NULL)
    }

    ## ** strata variables
    out$strataspecials <- special.object$strata
    out$is.strata <- length(ls.specials[[special.object$strata]])>0

    if(out$is.strata){
        if(length(ls.specials[[special.object$strata]])!=n.xterms){
            out$strata.sterms <- stats::drop.terms(xterms,setdiff(1:n.xterms,ls.specials[[special.object$strata]]))
        }else{
            out$strata.sterms <- xterms
        }
        out$stratavars <- attr(out$strata.sterms, "term.labels")
        out$stratavars.original <- all.vars(out$strata.sterms)
        out$strata.levels <- coxStrataLevel(object)

        if(length(out$strata.levels)==1){
            out$is.strata <- FALSE
        }
    }else{
        out["stratavars"] <- list(NULL)
        out["stratavars.original"] <- list(NULL)
        out["strata.sterms"] <- list(NULL)
        out$strata.levels <- factor("1") ## not NULL - coherent with coxStrata
    }

    ## ** export
    return(out)
} 
# }}}


#### methods #####

                                        # {{{ coxBaseEstimator


## * coxBaseEstimator
#' @title Extract the type of estimator for the baseline hazard
#' @description Extract the type of estimator for the baseline hazard
#' @name coxBaseEstimator 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxBaseEstimator
#' @export
coxBaseEstimator <- function(object){
  UseMethod("coxBaseEstimator") 
} 

## ** coxBaseEstimator.coxph
#' @rdname coxBaseEstimator
#' @method coxBaseEstimator coxph
#' @export
coxBaseEstimator.coxph <- function(object){
  return(object$method)
}

## ** coxBaseEstimator.phreg
#' @rdname coxBaseEstimator
#' @method coxBaseEstimator phreg
#' @export
coxBaseEstimator.phreg <- function(object){
  return("breslow")
}
# }}}

                                        # {{{ coxCenter
## * coxCenter
#' @title Extract the mean value of the covariates
#' @description Extract the mean value of the covariates
#' @name coxCenter 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxCenter
#' @export
coxCenter <- function(object){
  UseMethod("coxCenter") 
} 

## ** coxCenter.cph
#' @rdname coxCenter
#' @method coxCenter cph
#' @export
coxCenter.cph <- function(object){
  return(setNames(object$means, object$mmcolnames))
}

## ** coxCenter.coxph
#' @rdname coxCenter
#' @method coxCenter coxph
#' @export
coxCenter.coxph <- function(object){
  return(setNames(object$means, names(coef(object))))
}

## ** coxCenter.phreg
#' @rdname coxCenter
#' @method coxCenter phreg
#' @export
coxCenter.phreg <- function(object){
    return(apply(object$X,2,mean))
}
                                        # }}}

                                        # {{{ coxModelFrame
## * coxModelFrame
#' @title Extract the design matrix used to train a Cox model
#' @description Extract the design matrix used to train a Cox model. Should contain the time of event, the type of event, 
#' the variable for the linear predictor, the strata variables and the date of entry (in case of delayed entry).
#' @name coxModelFrame 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' @param center [logical] Should the variables of the linear predictor be added ?
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxModelFrame
#' @export
coxModelFrame <- function(object, center){
  UseMethod("coxModelFrame") 
} 

## ** coxModelFrame.coxph
#' @rdname coxModelFrame
#' @method coxModelFrame coxph
#' @export
coxModelFrame.coxph <- function(object, center = FALSE){

    default.start <- 0
    
    if("y" %in% names(object) == FALSE){
        stop("invalid object \n",
             "set y=TRUE in the call to ",class(object)[1]," \n")
    }
    
    if("x" %in% names(object) == FALSE){
        rhs.formula <- update(formula(object),"0~.")
        if(length(all.vars(rhs.formula))==0 || (!is.null(object$nevent) && object$nevent==0)){
            ## name.var <- setdiff(all.vars(object$terms), unlist(SurvResponseVar(coxFormula(object))[-1]))
            name.var <- names(object$coef)
            object$x <- matrix(0, nrow = NROW(object$y), ncol = length(name.var),
                               dimnames = list(NULL, name.var))
        }else{
            stop("invalid object \n",
                 "set x=TRUE in the call to ",class(object)[1]," \n")
        }
    }
 
    ## ** add x
    if(NCOL(object[["x"]])!=0){
        if(center){
            dt <- as.data.table(rowCenter_cpp(object[["x"]], center = coxCenter(object)))
        }else{
            dt <- as.data.table(object[["x"]])
        }        
        if(any(names(dt) != names(coef(object))) && NCOL(dt) == length(coef(object))){
            colnames(dt) <- names(coef(object))
        }
    }else{
        dt <- NULL
    }
    if("strata" %in% names(dt)){
        stop("The variables in the linear predictor should be named \"strata\" \n")
    }

    ## ** add y
    if(is.null(dt)){
        dt <- data.table(status = object[["y"]][,"status"])
    }else{
        dt[,c("status") := object[["y"]][,"status"]]
    }
    if("start" %in% colnames(object$y) == FALSE){
        dt[,c("start","stop") := list(default.start, object[["y"]][,"time"])]
    }else{        
        dt[,c("start","stop") := list(object[["y"]][,"start"], object[["y"]][,"stop"])]
    }
  
    ## ** add strata
    if("strata" %in% names(object)){
        dt[,c("strata") := object[["strata"]]]
    }else{
        dt[,c("strata") := factor(1)]
    }
    
    ## ** export
    first.col <- c("start","stop","status")
    data.table::setcolorder(dt, c(first.col,setdiff(names(dt),first.col)))
    return(dt)
}

## ** coxModelFrame.cph
#' @rdname coxModelFrame
#' @method coxModelFrame cph
#' @export
coxModelFrame.cph <- coxModelFrame.coxph

## ** coxModelFrame.phreg
#' @rdname coxModelFrame
#' @method coxModelFrame phreg
#' @export
coxModelFrame.phreg <- function(object, center = FALSE){

    default.start <- 0

    ## ** add y
    dt <- as.data.table(unclass(object$model.frame[,1]))
    if("start" %in% names(dt) == FALSE){
        dt[, c("start") := default.start]
    }
  
    ## normalize names
    name.old <- names(dt)
    name.new <- gsub("entry","start",gsub("time","stop",name.old))
    setnames(dt, old = name.old, new = name.new)
  
    ## ** add x
    if(center){
        M.X <- rowCenter_cpp(object$X, center = coxCenter(object))
        colnames(M.X) <- colnames(object$X)
        dt <- cbind(dt, as.data.table(M.X))
    }else{
        dt <- cbind(dt, as.data.table(object$X))
    }
    if("strata" %in% names(dt)){
        stop("The variables in the linear predictor should be named \"strata\" \n")
    }

    ## ** add strata
    if(!is.null(object$strata.name)){
        dt[,c("strata") := as.factor(object$model.frame[[object$strata.name]])]
    }else{
        dt[,c("strata") := factor(1)]
    }

    ## ** export
    first.col <- c("start","stop","status")
    data.table::setcolorder(dt, c(first.col,setdiff(names(dt),first.col)))
    return(dt)
    
}
# }}}

                                        # {{{ coxFormula
## * coxFormula
#' @title Extract the formula from a Cox model
#' @description Extract the formula from a Cox model
#' @name coxFormula 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxFormula
#' @export
coxFormula <- function(object){
  UseMethod("coxFormula") 
} 

## ** coxFormula.cph
#' @rdname coxFormula
#' @method coxFormula cph
#' @export
coxFormula.cph <- function(object){
  return(object$sformula)
}

## ** coxFormula.coxph
#' @rdname coxFormula
#' @method coxFormula coxph
#' @export
coxFormula.coxph <- function(object){
    if(object$nevent>0){
        return(object$formula)
    }else{
        out <- object$terms
        extra.attr <- names(attributes(out))
        for(iAttr in extra.attr){
            attr(out, iAttr) <- NULL
        }
        return(as.formula(out))
    }
}

## ** coxFormula.phreg
#' @rdname coxFormula
#' @method coxFormula phreg
#' @export
coxFormula.phreg <- function(object){
  return(object$formula)
}
## ** coxFormula.glm
#' @rdname coxFormula
#' @method coxFormula glm
#' @export
coxFormula.glm <- function(object){
  return(stats::formula(object))
}

# }}}

# {{{ coxLP
## * coxLP
#' @title Compute the linear predictor of a Cox model
#' @description Compute the linear predictor of a Cox model
#' @name coxLP 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param center should the linear predictor be computed after centering the covariates
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @details In case of empty linear predictor returns a vector of 0 with the same length as the number of rows of the dataset

#' @rdname coxLP
#' @export
coxLP <- function(object, data, center){
  UseMethod("coxLP") 
} 

## ** coxLP.cph
#' @rdname coxLP
#' @method coxLP cph
#' @export
coxLP.cph <- function(object, data, center){
  
    coef <- stats::coef(object)
    n.varLP <- length(coef)

    if(n.varLP==0){
        if(is.null(data)){
            Xb <- rep(0, coxN(object))
        }else{
            Xb <- rep(0, NROW(data))
        }
    }else if(is.null(data)){ ## training dataset
    
        Xb <- object$linear.predictors
    
        if(center[[1]] == FALSE && n.varLP != 0){
            Xb <- Xb + sum(coxCenter(object)*coef)
        }
    
    }else{ ## new dataset
        ## if(all(object$Design$assume %in% c("category","asis","interaction","strata"))){
        ## if(!("rcspline" %in% object$Design$assume)){
        Xb <- stats::predict(object, newdata = as.data.frame(data), type = "lp")
        ## }else{ ## does not work with splines
        ## X <- model.matrix(object, data = data)
        ## Xb <- as.vector(X %*% coef)
        ## }
        
      
      if(center == FALSE){
        Xb <- Xb + sum(coxCenter(object)*coef)
      }
      
  }
  
  return(unname(Xb))
}

## ** coxLP.coxph
#' @rdname coxLP
#' @method coxLP coxph
#' @export
coxLP.coxph <- function(object, data, center){
  
  coef <- stats::coef(object)
  n.varLP <- length(coef)
  
  if(is.null(data)){ ## training dataset
    
    Xb <- object$linear.predictors
    
    if(center[[1]] == FALSE && n.varLP != 0){
      Xb <- Xb + sum(coxCenter(object)*coef)
    }
    
  }else{ ## new dataset
      if(n.varLP>0){
      is.strata <- attr(object$terms, "special")$strata

      
      if(length(is.strata)>0){
        object.strata <- object[["strata"]]
        object[["strata"]] <- NULL # solve a bug in survival:::predict.coxph when fitting the model with strata

        Xb <- try(rowSums(stats::predict(object, newdata = as.data.frame(data), 
                                         type = "terms")), silent = TRUE)
        if(inherits(x=Xb,what="try-error")){ ## Fix an error when the dataset used to fit the object is removed from the global environment
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

## ** coxLP.phreg
#' @rdname coxLP
#' @method coxLP phreg
#' @export
coxLP.phreg <- function(object, data, center){

    coef <- stats::coef(object)
    n.varLP <- length(coef)

    if(n.varLP>0){
        if(is.null(data)){ ## training dataset
            X <- object$X
        }else{
            X <- model.matrix(object, data = data)
        }
  
        Xb <- as.vector(X[,names(coef), drop = FALSE] %*% coef)
    
        if(center){
            Xb <- Xb - sum(coxCenter(object)*coef)
        }
    
    }else{
        if(is.null(data)){
            Xb <- rep(0, coxN(object))
        }else{
            Xb <- rep(0, NROW(data))
        }
    } 
  
  
  return(Xb)
}

# }}}

                                        # {{{ coxN
## * coxN
#' @title Extract the number of observations from a Cox model
#' @description Extract the number of observations from a Cox model
#' @name coxN 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#'

#' @rdname coxN
#' @export
coxN <- function(object){
  UseMethod("coxN") 
} 

## ** coxN.cph
#' @rdname coxN
#' @method coxN cph
#' @export
coxN.cph <- function(object){
    return(sum(object$n))
}

## ** coxN.coxph
#' @rdname coxN
#' @method coxN coxph
#' @export
coxN.coxph <- function(object){
  return(object$n)
}

## ** coxN.phreg
#' @rdname coxN
#' @method coxN phreg
#' @export
coxN.phreg <- function(object){
  return(NROW(object$model.frame))
}

## ** coxN.CSC
#' @rdname coxN
#' @method coxN CauseSpecificCox
#' @export
coxN.CauseSpecificCox <- function(object){
  return(sapply(object$models,coxN))
}

## ** coxN.CSC
#' @rdname coxN
#' @method coxN glm
#' @export
coxN.glm <- function(object){
  return(stats::nobs(object))
}

# }}}

                                        # {{{ coxSpecialStrata
## * coxSpecial
#' @title Special characters in Cox model
#' @description Return the special character(s) of the Cox model, e.g. used to indicate the strata variables.
#' @name coxSpecial
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#'
#' @details Must return a list with at least one element strata
#' indicating the character in the formula marking the variable(s) defining the strata.
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxSpecial
#' @export
coxSpecial <- function(object) UseMethod("coxSpecial")

## ** coxSpecial.coxph
#' @rdname coxSpecial
#' @method coxSpecial coxph
#' @export
coxSpecial.coxph <- function(object){
    return(list(strata = "strata",
                cluster = "cluster"))
}

## ** coxSpecial.cph
#' @rdname coxSpecial
#' @method coxSpecial cph
#' @export
coxSpecial.cph <- function(object){
  return(list(strata = "strat"))
}

## ** coxSpecial.phreg
#' @rdname coxSpecial
#' @method coxSpecial phreg
#' @export
coxSpecial.phreg <- function(object){
    return(list(strata = "strata",
                cluster = "cluster"))
}
# }}}

## * coxStrataLevel
#' @title Returns the name of the strata in Cox model
#' @description Return the name of the strata in Cox model
#' @name coxStrataLevel
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#'
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxStrataLevel
#' @export
coxStrataLevel <- function(object) UseMethod("coxStrataLevel")

## ** coxStrataLevel.coxph
#' @rdname coxStrataLevel
#' @method coxStrataLevel coxph
#' @export
coxStrataLevel.coxph <- function(object){
    if(!is.null(object$strata)){
        return(levels(object$strata))
    }else{
        return(NULL)
    }
}

## ** coxStrataLevel.cph
#' @rdname coxStrataLevel
#' @method coxStrataLevel cph
#' @export
coxStrataLevel.cph <- function(object){
    if(!is.null(object$strata)){
        return(levels(object$strata))
    }else{
        return(NULL)
    }
}

## ** coxStrataLevel.phreg
#' @rdname coxStrataLevel
#' @method coxStrataLevel phreg
#' @export
coxStrataLevel.phreg <- function(object){
    if(!is.null(object$strata.name)){
       return( levels(object$model.frame[[object$strata.name]]) )
    }else{
        return( NULL )
    }
}
                                        # }}}

                                        # {{{ coxStrata
## * coxStrata
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
# '@export
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

                                        # {{{ coxVarCov
## * coxVarCov
#' @title Extract the variance covariance matrix of the beta from a Cox model
#' @description Extract the variance covariance matrix of the beta from a Cox model
#' @name coxVarCov 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' 
#' @author Brice Ozenne broz@@sund.ku.dk
#' 
#' @details Should return \code{NULL} if the Cox model has no covariate. 
#' The rows and columns of the variance covariance matrix must be named with the names used in the design matrix.

#' @rdname coxVarCov
#' @export
coxVarCov <- function(object){
  UseMethod("coxVarCov") 
} 

## ** coxVarCov.cph
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

## ** coxVarCov.coxph
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

## ** coxVarCov.phreg
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

## * SurvResponseVar
#' @title Extract the time and event variable from a Cox model
#' @description Extract the time and event variable from a Cox model
#' @name SurvResponseVar
#' 
#' @param formula a formula
#'
#' @author Brice Ozenne broz@@sund.ku.dk
#'
#' @examples
#' \dontrun{
#' SurvResponseVar(Surv(time,event)~X1+X2)
#' SurvResponseVar(Hist(time,event==0)~X1+X2)
#' SurvResponseVar(Surv(start,time, status,type="counting") ~ X3+X5)
#' SurvResponseVar(Surv(start,event=status, time2=time,type="counting") ~ X3+X5)
#'
#' SurvResponseVar(survival::Surv(start,event=status, time2=time,type="counting") ~ X3+X5)
#' SurvResponseVar(status ~ X3+X5)
#' SurvResponseVar(I(status == 1) ~ X3+X5)
#' SurvResponseVar(list(Hist(time, event) ~ X1+X6,Hist(time, event) ~ X6))
#' }
SurvResponseVar <- function(formula){

    if(inherits(x=formula,what="call")){
        formula <- eval(formula)
    }
    
    if(inherits(x=formula,what="list")){
        if(any(unlist(lapply(formula,class))!="formula")){
            stop("When a list, the elements of argument \'formula\' must all be a formula \n")
        }
        
        ls.formula <- lapply(formula, SurvResponseVar)
        if(length(ls.formula)>1){
            test.identical <- sapply(ls.formula[-1], function(x){identical(ls.formula[[1]],x)})
            if(test.identical){
                return(ls.formula[[1]])
            }else{
                stop("Different left hand side in the formulae \n")
            }
        }else{
            return(ls.formula[[1]])
        }
        
    }else if(!inherits(formula,"formula")){
        stop("Argument \'formula\' can either be a formula or a list of formula \n")
    }

    ls.formula <- as.list(formula[[2]])
    if(length(ls.formula)==1){ ## deal with the case of simple formula, e.g. Y~X
        return(list(status = all.vars(formula)[1]))
    }
    
    ## extra names of the inputs
    n.input <- length(ls.formula) - 1
    all.input <- names(ls.formula)[-1]
    if(is.null(all.input)){
        all.input <- rep("",n.input)
    }
    

    ## extract arguments
    operator <- deparse(ls.formula[[1]])
    if(operator %in% c("survival::Surv","Surv")){
        ## names(as.list(args(survival::Surv)))
        if( ("time2" %in% all.input) || ( sum(all.input=="") >= (1 + ("time" %in% all.input == FALSE) + ("event" %in% all.input == FALSE)) )){
            all.args <- c(entry = "time", time = "time2", status = "event", type = "type", origin = "origin")
        }else{
            all.args <- c(time = "time", status = "event", type = "type", origin = "origin")
        }
        
    }else if(operator %in% c("riskRegression::Hist","Hist")){
        ## names(as.list(args(Hist)))
        all.args <- c(time = "time", status = "event", entry = "entry", id = "id", cens.code = "cens.code", addInitialState = "addInitialState")
    }else if(operator %in% "I"){
        xx <- all.vars(formula)[1]
        attr(xx, "call") <- deparse(formula[[2]])
        return(list(status = xx))
    }else{
        stop("The left side of the formula must contain Hist() or Surv() \n") 
    }

    ## associate named inputs
    remaining.input <- all.input
    remaining.args <- all.args
    index.input <- 1:n.input
    if(any(all.input!="")){
        for(iNames in all.input[all.input!=""]){
            remaining.args <- remaining.args[remaining.args != iNames]
            index.input <- index.input[remaining.input != iNames]
            remaining.input <- remaining.input[remaining.input != iNames]
        }
    }
    ## set names to unnamed inputs
    all.input[index.input] <- remaining.args[1:length(index.input)]

    all.args <- all.args[all.args %in% all.input]
    for(iNewname in 1:length(all.args)){ ## iNewname <- 1
        all.input[all.input == all.args[iNewname]] <- names(all.args)[iNewname]
    }
    
    ## output
    names(ls.formula) <- c("operator",all.input)
    
    for(iArg in 1:length(ls.formula)){ ## iArg <- 3
        if(inherits(x=ls.formula[[iArg]],what="name")){
            ls.formula[[iArg]] <- deparse(ls.formula[[iArg]])
        }else if(inherits(x=ls.formula[[iArg]],what="call")){
            xx <- as.character(ls.formula[[iArg]])[2]
            attr(xx,"call") <- deparse(ls.formula[[iArg]])
            ls.formula[[iArg]] <- xx
        }
    }

  ## export
  return(ls.formula)
  
}
# }}}

                                        # {{{ splitStrataVar
## * splitStrataVar
#' @title Reconstruct each of the strata variables
#' @description Reconstruct each of the strata variables from the strata variable stored in the coxph object.
#' @name splitStrataVar
#' 
#' @param object a coxph object.
#'
#' @author Brice Ozenne broz@@sund.ku.dk and Thomas A. Gerds tag@@biostat.ku.dk
#'
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

# }}}

                                        # {{{ model.matrix.phreg

## * model.matrix
## ** model.matrix.cph
#' @title Extract design matrix for cph objects
#' @description Extract design matrix for cph objects
#' @param object a cph object.
#' @param data a dataset.
#' 
#' @method model.matrix cph
model.matrix.cph <- function(object, data){

    if(all(object$Design$assume %in% c("category","asis","interaction","strata"))){
        out <- predict(object, newdata = data, type = "x")

        ## put back rms names
        colnames(out) <- object$Design$colnames[match(object$Design$mmcolnames,colnames(out))]
    }else{
        type.special <- setdiff(object$Design$assume,c("category","asis","interaction","strata"))
        name.special <- attr(object$terms,"term.labels")[object$Design$assume %in% type.special]
        ## warning("Special operator(s) in the formula: ",paste(name.special,collapse = ", "),".\n",
        ##         "Associated type: \"",paste(type.special,collapse = "\" \""),"\". \n",
        ##         "If it is a spline term, consider using survival::coxph instead of rms::cph. \n")

        ## set categorical variables to their appropriate level
        if(any(object$Design$assume=="category")){ 
            for(iVar in object$Design$name[object$Design$assume=="category"]){
                if(!is.factor(data[[iVar]]) || any(sort(levels(data[[iVar]])) != sort(object$Design$parms[[iVar]])) ){
                    data[[iVar]] <- factor(data[[iVar]], levels=object$Design$parms[[iVar]])
                }
            }
        }

        ## borrowed from survival::coxph
        mf <- stats::model.frame(object)
        Terms <- stats::delete.response(terms(mf))
        out <- model.matrix(Terms, data)

        ## remove intercept
        out <- out[,attr(out,"assign")!=0,drop=FALSE]

        ## put back rms names
        if(!is.null(attr(object$Design$mmcolnames,"alt"))){
            colnames(out) <- object$Design$colnames[match(attr(object$Design$mmcolnames,"alt"),colnames(out))]
        }else{
            colnames(out) <- object$Design$colnames[match(object$Design$mmcolnames,colnames(out))]
        }
        ## add strata if possible
        try(attr(out,"strata") <- attr(predict(object, newdata = data, type = "x"),"strata"), silent = TRUE)
    }
    return(out)
    
}

## ** model.matrix.phreg
#' @title Extract design matrix for phreg objects
#' @description Extract design matrix for phreg objects
#' @param object a phreg object.
#' @param data a dataset.
#' 
#' @details mainly a copy paste of the begining of the \code{phreg} function.
#' 
#' @method model.matrix phreg
model.matrix.phreg <- function(object, data){
    special <- c("strata", "cluster")
    Terms <- stats::delete.response(attr(object$model.frame, "terms"))

    ## remove specials
    if (length(attributes(Terms)$specials$cluster)>0) {
        ts <- survival::untangle.specials(Terms, "cluster")
        Terms <- Terms[-ts$terms]
    }
    if (length(attributes(Terms)$specials$strata)>0) {
        ts <- survival::untangle.specials(Terms, "strata")
        Terms <- Terms[-ts$terms]
    }
    attr(Terms,"intercept") <- 1 ## keep intercept to have the same behavior with and without categorical variables 
    X <- model.matrix(Terms, data)
    return(X[,setdiff(colnames(X),"(Intercept)"),drop=FALSE])
}

# }}}

## * terms
## ** terms.phreg
#' @title Extract terms for phreg objects
#' @description Extract terms for phreg objects
#' @param x a phreg object.
#' @param ... not used.
#' 
#' @method terms phreg
terms.phreg <- function(x, ...){
    stats::terms(x$formula)
}


