#' @title Extract information from a Cox model
#' @description Extract information from a Cox model that are necessary for computing the baseline hazard, performing predictions or estimating the influence function associated to the Cox model
#' @name getCoxInfo
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param design should the design matrix be returned? 
#' @param center should the covariates be centered
#' 
#' @return A named list containing the following elements:
#' \itemize{
#'  \item{"nPatients"}{the number of observations}
#'  \item{"xvars"}{the name of all the regressors (including those used to form the strata)}
#'  \item{"modeldata"}{the design matrix for the regressors}
#'  \item{"stratavars"}{the name of the variables used to define the strata}
#'  \item{"is.strata"}{is there any strata?}
#'  \item{"sterms"}{terms corresponding to the strata variables}
#'  \item{"strataF"}{a vector contain the strata factor for each observation}
#'  \item{"stratalevels"}{a named list containing for each variable used to form the strata all its possible levels}
#' }
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
#' getCoxInfo(mCox, design = FALSE) 
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' getCoxInfo(mCoxS, design = FALSE) 
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' getCoxInfo(mCox, design = FALSE) 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' getCoxInfo(mCoxS, design = FALSE) 
#' }

#' @rdname getCoxInfo
getCoxInfo <- function(object, design, center) UseMethod("getCoxInfo")

#' @rdname getCoxInfo
#' @method getCoxInfo coxph
getCoxInfo.cph <- function(object, design, center){
  
  if(is.null(object$y)){
    stop("Argument \'y\' must be set to TRUE in cph \n")
  }
  
  nPatients <- sum(object$n)
  infoVar <- getVar(object)
  is.strata <- infoVar$strata$is.strata
  
  #### strata
  if(is.strata){
    strataF <- object$strata
  }else{
    strataF <- factor("1")
  }
  
  if(design){ ## cph:design matrix for standard error
    modeldata <- getData(object, data = model.frame(object), lpVars = infoVar$lp$vars, stratavars = stratavars, center = center)
  }else{
    modeldata <- matrix(nrow = 0, ncol = 0)
  }
  
  return(list(nPatients = nPatients,
              lpvars = infoVar$lp$vars,
              stratavars = infoVar$strata$stratavars,
              strataOvars = infoVar$strata$stratavars.original,
              timevars = infoVar$response$time,
              eventvars = infoVar$response$event,
              is.strata = infoVar$strata$is.strata,
              stratalevels = infoVar$strata$stratalevels,
              sterms = infoVar$strata$sterms,
              modeldata = modeldata,
              strataF = strataF))
}

#' @rdname getCoxInfo
#' @method getCoxInfo coxph
getCoxInfo.coxph <- function(object, design, center){
  
  nPatients <- object$n
  infoVar <- getVar(object)
  
  is.strata <- infoVar$strata$is.strata
  stratavars <- infoVar$strata$stratavars
    
  #### strata
  if(is.strata){
    strataF <- interaction(stats::model.frame(object)[,stratavars], drop = TRUE, sep = ", ", lex.order = TRUE) 
  }else{
    strataF <- factor("1")
  }
  
  if(design){ ## coxph:design matrix for standard error
    modeldata <- centerData(object, data = model.frame(object), stratavars = stratavars, center = center)
  }else{
    modeldata <- matrix(nrow = 0, ncol = 0)
  }
  
  return(list(nPatients = nPatients,
              lpvars = infoVar$lp$vars,
              stratavars = infoVar$strata$stratavars,
              strataOvars = infoVar$strata$stratavars.original,
              timevars = infoVar$response$time,
              eventvars = infoVar$response$event,
              is.strata = infoVar$strata$is.strata,
              stratalevels = infoVar$strata$stratalevels,
              sterms = infoVar$strata$sterms,
              modeldata = modeldata,
              strataF = strataF
  ))
}

#' @title Extract variable names from a model
#' @description Extract the name of the variables belonging to the linear predictor or used to form the strata
#' @rdname getVar 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param lp should the name of the variable contained in the linear predictor be output
#' @param strata should the variables used to form the strata be output
#' @param response should the name of the response variable be output
#' 
#' @return a list with three elements, one for the linear predictor, one for the strata variables and one for the response variables
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
#' getVar(mCox)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d)
#' getVar(mCoxS, response = FALSE)
#' mCoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' getVar(mCoxS2, response = FALSE)
#' mCoxI <- coxph(Surv(time, event) ~ X1*X2, data = d)
#' getVar(mCoxI, response = FALSE)
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' getVar(mCox) 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d)
#' getVar(mCoxS, response = FALSE) 
#' mCoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' getVar(mCoxS2, response = FALSE) 
#' }

#' @rdname getVar
getVar <- function(object, lp, strata, response){
  UseMethod("getVar") 
} 

#' @rdname getVar
#' @method getVar coxph
getVar.coxph <- function(object, lp = TRUE, strata = TRUE, response = TRUE){
  
  if(strata || lp){
    xterms <- delete.response(object$terms)
    allVars.X <- all.vars(xterms)
    
    strataspecials <- attr(xterms,"specials")$strata
    stratavars.original <- allVars.X[strataspecials]
  }
  
  ## strata
  ls.strata <- list()
  
  if(strata){
    xvars <- attr(xterms,"term.labels")
    stratavars <- xvars[strataspecials]
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
    
    ls.strata <- list(stratavars = stratavars,
                      stratavars.original = stratavars.original,
                      strataspecials = strataspecials,
                      stratalevels = stratalevels,
                      is.strata = is.strata,
                      sterms = sterms)
  }
  
  ## LP
  ls.lp <- list()
  if(lp){
    
    ls.lp <- list(vars = setdiff(allVars.X,stratavars.original),
                  test = !all(allVars.X %in% stratavars.original))
    
  }
  
  ## Response
  ls.response <- list()
  if(response){
    
    f.surv <- object$call$formula[[2]]
    if("time" %in% names(f.surv) == FALSE && "event" %in% names(f.surv) == FALSE ){
      names(f.surv)[2:3] <- c("time","event")
    }else if("time" %in% names(f.surv) == FALSE){
      names(f.surv)[setdiff(2:3,which(names(f.surv)=="event"))] <- "time"
    }else if("event" %in% names(f.surv) == FALSE){
      names(f.surv)[setdiff(2:3,which(names(f.surv)=="time"))] <- "event"
    }
    
    ls.response <- list(time = deparse(f.surv$time),
                        event = deparse(f.surv$event))
    
  }
  
  ## export
  return(list(strata = ls.strata,
              lp = ls.lp,
              response = ls.response
  ))
}

#' @rdname getVar
#' @method getVar cph
getVar.cph <- function(object, lp = TRUE, strata = TRUE, response = TRUE){
  
  if(strata || lp){
    xterms <- delete.response(object$terms)
    allVars.X <- all.vars(xterms)
    
    strataspecials <- attr(xterms,"specials")$strat
    stratavars.original <- allVars.X[strataspecials]
  }
  
  ## strata
  ls.strata <- list()
  
  if(strata){
    xvars <- attr(xterms,"term.labels")
    stratavars <- xvars[strataspecials]
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
    
    ls.strata <- list(stratavars = stratavars,
                      stratavars.original = stratavars.original,
                      strataspecials = strataspecials,
                      stratalevels = NULL,
                      is.strata = is.strata,
                      sterms = sterms)
  }
  
  ## LP
  ls.lp <- list()
  if(lp){
    
    ls.lp <- list(vars = setdiff(allVars.X,stratavars.original),
                  test = !all(allVars.X %in% stratavars.original))
    
  }
  
  ## Response
  ls.response <- list()
  if(response){
    
    f.surv <- object$call$formula[[2]]
    if("time" %in% names(f.surv) == FALSE && "event" %in% names(f.surv) == FALSE ){
      names(f.surv)[2:3] <- c("time","event")
    }else if("time" %in% names(f.surv) == FALSE){
      names(f.surv)[setdiff(2:3,which(names(f.surv)=="event"))] <- "time"
    }else if("event" %in% names(f.surv) == FALSE){
      names(f.surv)[setdiff(2:3,which(names(f.surv)=="time"))] <- "event"
    }
    
    ls.response <- list(time = deparse(f.surv$time),
                        event = deparse(f.surv$event))
    
  }
  
  ## export
  return(list(strata = ls.strata,
              lp = ls.lp,
              response = ls.response
  ))
}

#' @title Extract data
#' @description Extract the columns from a dataset involved by a model. Optional argument to center the columns values.
#' @name getData
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param stratavars the name of the strata variables
#' @param center should the covariates be centered
#' 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCoxS <- coxph(Surv(time, event) ~ X1, data = d) 
#' getData(mCoxS, lpVars = "X1", stratavars = NULL, center = FALSE)
#' mCoxS <- coxph(Surv(time, event) ~ X1*X2, data = d) 
#' getData(mCoxS, lpVars = c("X1","X2"), stratavars = NULL, center = FALSE)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d) 
#' getData(mCoxS, lpVars = "X2", stratavars = "strata(X1)", center = FALSE)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d) 
#' getData(mCoxS, lpVars = NULL, stratavars = c("strata(X1)","strata(X2)"), center = FALSE)
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ X1, data = d) 
#' getData(mCoxS, lpVars = "X1", data = d, stratavars = NULL, center = FALSE)
#' mCoxS <- cph(Surv(time, event) ~ X1*X2, data = d) 
#' getData(mCoxS, lpVars = c("X1","X2"), data = d, stratavars = NULL, center = FALSE)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d) 
#' xx <- getData(mCoxS, lpVars = "X2", data = d, stratavars = "X1", center = FALSE)
#' print(xx)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d) 
#' xx <- getData(mCoxS, lpVars = NULL, data = d, stratavars = c("X1","X2"), center = FALSE)
#' print(xx)
#' }
getData <- function(object, data = model.frame(object), lpVars, stratavars, center){
  
  if(is.data.table(data)){
    modeldata <- data.table()
  }else{
    modeldata <- data.frame(matrix(ncol = 0, nrow = NROW(data)))
  }
  
  if(length(lpVars)>0){
    
    if(is.data.table(data)){
      modeldata <- data[,lpVars,with=FALSE]
    }else{
      modeldata <- data[,lpVars,drop=FALSE]
    }
    if(center){
      modeldata <- sweep(modeldata, FUN = "-", MARGIN = 2, STATS = object$means[lpVars])
    }
    
  }
  
  if(length(stratavars)>0){
    if(is.data.table(data)){
      if(NROW(modeldata)==0){
        modeldata <- data[,stratavars,with=FALSE]
      }else{
        modeldata[,stratavars := data[,stratavars,with=FALSE], with = FALSE]
      }
    }else{
      modeldata <- cbind(modeldata, 
                         data[,stratavars,drop=FALSE])
    }
  }
  
  return(modeldata)
  
} 

#' @title Compute the linear predictor of a Cox model
#' @description Compute the linear predictor of a Cox model
#' @rdname getLP 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param xvars the name of all the regressors (including those used to form the strata)
#' @param stratavars the variables used to define the strata
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
#' lpCox(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL, center = FALSE)
#' lpCox(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL, center = TRUE)  
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+X2, data = d)
#' lpCox(mCoxS, data = d, xvars = c("X1","X2"), stratavars = "strata(X1)", center = FALSE) 
#' mCoxS2 <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' lpCox(mCoxS2, data = d, xvars = c("X1","X2"), stratavars = c("X1","X2"), center = FALSE) 
#' 
#' ##
#' library(rms)
#' mCox <- cph(Surv(time, event) ~ X1+X2, data = d, y = TRUE)
#' lpCox(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL, center = FALSE) 
#' lpCox(mCox, data = d, xvars = c("X1","X2"), stratavars = NULL, center = TRUE) 
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+X2, data = d)
#' lpCox(mCoxS, data = d, xvars = c("X1","X2"), stratavars = "strat(X1)", center = FALSE) 
#' mCoxS2 <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' lpCox(mCoxS2, data = d, xvars = c("X1","X2"), stratavars = c("X1","X2"), center = FALSE) 
#' }

#' @rdname getLP
getLP <- function(object, data, lpVars, stratavars, center){
  UseMethod("getLP") 
} 

#' @rdname getLP
#' @method getLP cph
getLP.cph <- function(object, data, lpVars, stratavars, center){
  
  if(length(lpVars)>0){
    data <- getData(object, data = data, lpVars = lpVars, stratavars = stratavars, center = center)
    Xb <- stats::predict(object, data, type = "lp")
  }else{ 
    Xb <- rep(0, NROW(data)) 
  }
  
  return(Xb)
}

#' @rdname getLP
#' @method getLP coxph
getLP.coxph <- function(object, data, lpVars, stratavars, center){
  
  if(length(lpVars) == 0){ 
    Xb <- rep(0, NROW(data))
  } else{
    data <- getData(object, data = data, lpVars = lpVars, stratavars = stratavars, center = center)
    if(length(stratavars)>0){
      Xb <- rowSums(stats::predict(object, newdata = data, type = "terms"))
    }else { 
      Xb <- stats::predict(object, data, type = "lp")
    }
  } 
  
  return(Xb)
}

#' @title Define the strata for a new dataset
#' @description Define the strata in a dataset to match those of a stratified Cox model
#' @name getStrata
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param data a \code{data.frame} or a \code{data.table}
#' @param sterms terms in the formula corresponding to the strata variables
#' @param stratavars the name of the variables used to define the strata
#' @param levelsStrata the strata levels that have been used to fit the Cox model
#' @param stratalevels a named list containing for each variable used to form the strata all its possible levels
#'
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' 
#' @examples 
#' \dontrun{
#' d <- sampleData(1e2, outcome = "survival")
#' 
#' ##
#' library(survival)
#' mCoxS <- coxph(Surv(time, event) ~ strata(X1)+strata(X2), data = d)
#' info <- getCoxInfo(mCoxS, design = FALSE)
#' getStrata(mCoxS, data = d, sterms = info$sterms, stratavars = info$stratavars, 
#'              levels = levels(info$strataF), stratalevels = info$stratalevels) 
#' 
#' ##
#' library(rms)
#' mCoxS <- cph(Surv(time, event) ~ strat(X1)+strat(X2), data = d, y = TRUE)
#' info <- getCoxInfo(mCoxS, design = FALSE)
#' getStrata(mCoxS, data = d, sterms = info$sterms, stratavars = info$stratavars, 
#'              levels = levels(info$strataF), stratalevels = info$stratalevels) 
#' }

#' @rdname getStrata
getStrata <- function(object, data, sterms, stratavars, levels, stratalevels) UseMethod("getStrata")

#' @rdname getStrata
#' @method getStrata coxph
getStrata.cph <- function(object, data, sterms, stratavars, levels, stratalevels){
  
  if(length(stratavars)>0){
    tmp <- model.frame(sterms,data)
    colnames(tmp) <- names(prodlim::parseSpecialNames(names(tmp),"strat"))
    tmp <- data.frame(lapply(1:NCOL(tmp),function(j){factor(paste0(names(tmp)[j],"=",tmp[,j,drop=TRUE]))}))
    newstrata <- apply(tmp,1,paste,collapse=".")
    newstrata <- factor(newstrata, levels = levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert newstrata to numeric
  }else{
    newstrata <- NULL
  } 
  
  return(newstrata)
}

#' @rdname getStrata
getStrata.coxph <- function(object, data, sterms, stratavars, levels, stratalevels){
  
  if(length(stratavars)>0){
    newstrata <- prodlim::model.design(sterms,data=data,xlev=stratalevels,specialsFactor=TRUE)$strata[[1]]
    newstrata <- factor(newstrata, levels = levels) # add all levels - necessary for predict.CauseSpecificCox to able to correctly convert newstrata to numeric
  }else{
    newstrata <- NULL
  } 
  
  return(newstrata)
}



