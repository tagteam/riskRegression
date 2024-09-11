### compress.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  9 2024 (14:04) 
## Version: 
## Last-Updated: sep 10 2024 (13:41) 
##           By: Brice Ozenne
##     Update #: 50
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * compressData
##' @description Find unique covariate sets and restrict the dataset to those.
##' Used by predictCox and predict.CauseSpecificCox
##' @noRd
compressData <- function(object, newdata, times, diag, average.iid,
                         oorder.times, times.sorted, nStrata, infoVar, level){


    ## ** extract information
    if(is.null(newdata) || (!is.null(level) && level == "full")){
        return(NULL)
    }else if(inherits(object,"prodlim")){
        n.cov <- 1 ## unique covariate profile (as there is no covariate so all predictions are equals up to the strata)
        if((!is.null(level) && level == "minimal")){
            nStrata <- 1
        }
        test.allCategorical <- TRUE
    }else if(inherits(object,"CauseSpecificCox")){
        ls.infoVar <- lapply(object$model, function(iO){coxVariableName(iO, model.frame = coxModelFrame(iO))})
        test.allCategorical <- all(sapply(ls.infoVar, function(iInfo){length(intersect(iInfo$lpvars, iInfo$lpvars.original))==0}))

        if(test.allCategorical || (!is.null(level) && level == "minimal")){
            ## *** covariates
            keep.col <- !duplicated(do.call("c",lapply(ls.infoVar, "[[","lpvars.original")))
            newdata.X <- do.call(cbind,lapply(object$model, FUN = stats::model.matrix, data = newdata))[,keep.col,drop=FALSE]
            if((!is.null(level) && level == "minimal")){
                n.cov <- 1
            }else{
                n.cov <- 2^NCOL(newdata.X)
            }
            ## *** strata variables
            allStrata.var <- unique(unlist(lapply(ls.infoVar,"[[","stratavars.original")))
            if(length(setdiff(allStrata.var,names(newdata.X)))>0){
                newdata.X <- cbind(newdata.X,as.data.frame(newdata)[setdiff(allStrata.var,names(newdata.X))])
            }
            if((!is.null(level) && level == "minimal")){
                nStrata <- 1
            }else if(sum(!duplicated(lapply(ls.infoVar,"[[","stratavars.original")))==1){
                nStrata <- length(ls.infoVar[[1]][["strata.levels"]])
            }else{ ## conservative approximation
                nStrata <- prod(lengths(lapply(ls.infoVar,"[[","strata.levels")))
            }            
        }else{
            n.cov <- Inf
            nStrata <- Inf
        }
    }else{
        test.allCategorical <- length(intersect(infoVar$lpvars, infoVar$lpvars.original))==0

        if(test.allCategorical || (!is.null(level) && level == "minimal")){
            ## *** covariates
            newdata.X <- model.matrix(object, data = newdata)
            if((!is.null(level) && level == "minimal")){
                n.cov <- 1
                nStrata <- 1
            }else{
                n.cov <- 2^NCOL(newdata.X)
            }
            ## *** strata variables
            if(length(setdiff(infoVar$stratavars.original,names(newdata.X)))>0){
                newdata.X <- cbind(newdata.X,as.data.frame(newdata)[setdiff(infoVar$stratavars.original,names(newdata.X))])
            }
            
        }else{
            n.cov <- Inf
        }
    }

    ## ** find unique patient profiles
    if((!test.allCategorical && is.null(level)) || NROW(newdata) <= nStrata*n.cov){
        return(NULL)
    }else if(inherits(object,"prodlim")){ ## Kaplan Meier
        newdata.save <- newdata

        ## find unique combinations of predictors
        newdata.strata <- coxStrata(object, data = newdata, 
                                    sterms = infoVar$strata.sterms, 
                                    strata.vars = infoVar$stratavars, 
                                    strata.levels = infoVar$strata.levels)
        ## associate each observation to a unique combination of predictors
        newdata.index <- tapply(1:NROW(newdata), INDEX = newdata.strata, FUN = identity, simplify = FALSE)
        attr(newdata.index,"vectorwise") <- as.numeric(newdata.strata)

    }else if(nStrata == 1 && n.cov == 1){ ## No covariate nor strata
        newdata.save <- newdata

        ## find unique combinations of predictors
        newdata.index <- list(1:NROW(newdata))

        ## associate each observation to a unique combination of predictors
        attr(newdata.index,"vectorwise") <- rep(1, NROW(newdata))

    }else{ ## Cox or Cause specific Cox

        ## associate each observation to a unique combination of predictors
        if(NCOL(newdata.X)==1){
            newdata.index <- tapply(1:NROW(newdata), INDEX = droplevels(newdata.X[,1]), FUN = identity, simplify = FALSE)
        }else if(all(newdata.X %in% 0:1) && NCOL(newdata.X)<=10){
            newdata.index <- tapply(1:NROW(newdata), INDEX = rowSums(newdata.X * matrix(10^(1:NCOL(newdata.X)), nrow = NROW(newdata.X), ncol = NCOL(newdata.X), byrow = TRUE)), FUN = identity, simplify = FALSE)
        }else{
            newdata.index <- tapply(1:NROW(newdata), INDEX = do.call(paste, args = c(as.data.frame(newdata.X), sep = ".")), FUN = identity, simplify = FALSE)
        }

        ## library(microbenchmark)
        ## microbenchmark(A = rowSums(newdata.X * 10^(1:NCOL(newdata.X))),
        ##                B = do.call(paste, args = c(as.data.frame(newdata.X), sep = ".")),
        ##                C = interaction(as.data.frame(newdata.X)))
        
        ## associate each observation to a unique combination of predictors
        attr(newdata.index,"vectorwise") <- unlist(lapply(1:length(newdata.index), function(iN){rep(iN,length(newdata.index[[iN]]))}))[order(unlist(newdata.index))]

    }

    ## ** compress data
    ## save full data
    newdata.save <- newdata
    ## subset data
    if(data.table::is.data.table(newdata)){
        newdata <- newdata.save[sapply(newdata.index,"[",1)]
    }else{
        newdata <- as.data.frame(newdata.save)[sapply(newdata.index,"[",1),,drop=FALSE]
    }
    ## do not do diag
    diag.save <- diag
    diag <- FALSE
    ## update times (mostly relevant when diag=TRUE to avoid duplicated times)
    oorder.times.save <- oorder.times
    times.sorted.save <- times.sorted

    times.sorted <- unique(times.sorted)
    order.times <- order(times.sorted)
    oorder.times <- order(order.times)
    nTimes <- length(times.sorted)

    if(average.iid==TRUE){
        if(is.null(attr(average.iid,"factor"))){
            attr(average.iid,"rm.list") <- TRUE
            if(diag.save){
                newdata.table <- table(attr(newdata.index,"vectorwise"),times.sorted.save[oorder.times.save])
                attr(average.iid,"factor") <- list(NROW(newdata)*matrix(as.double(newdata.table), nrow = NROW(newdata.table), ncol = NCOL(newdata.table))/NROW(newdata.save))
            }else{
                attr(average.iid,"factor") <- list(matrix(NROW(newdata)*lengths(newdata.index)/NROW(newdata.save), byrow = FALSE, nrow = NROW(newdata), ncol = nTimes))
            }
        }else{
            attr(average.iid,"rm.list") <- FALSE
            for(iF in 1:length(attr(average.iid,"factor"))){ ## iF <- 1
                attr(average.iid,"factor")[[iF]] <- do.call(rbind,lapply(newdata.index, function(iIndex){
                    if(diag.save){
                        iFactor <- tapply(attr(average.iid,"factor")[[iF]][iIndex,1],
                                          INDEX = factor(times.sorted.save[oorder.times.save][iIndex], levels = times.sorted),
                                          FUN = sum, default = 0)
                    }else{
                        iFactor <- colSums(attr(average.iid,"factor")[[iF]][iIndex,,drop=FALSE])
                    }
                    return(NROW(newdata)*iFactor/NROW(newdata.save))
                }))
            }
        }
    }

    ## ** export
    return(list(newdata = newdata, newdata.save = newdata.save, newdata.index = newdata.index,
                times.sorted = times.sorted, times.sorted.save = times.sorted.save,
                diag = diag, diag.save = diag.save,
                order.times = order.times,
                oorder.times = oorder.times, oorder.times.save = oorder.times.save,
                nTimes = nTimes,
                average.iid = average.iid))
}

## * decompressData
##' @description Expand the prediction and iid to the original data based on the results for the unique covariate sets.
##' Used by predictCox and predict.CauseSpecificCox
##' @noRd
decompressData <- function(object, newdata, type, diag, times, se, iid, average.iid,
                           newdata.index, times.sorted, needOrder){

    ## ** recover original profiles
    if(diag){
        for(iType in type){ ## iType <- "cumhazard"                
            object[[iType]] <- cbind(object[[iType]][attr(newdata.index,"vectorwise") + (match(times, times.sorted)-1) * NROW(newdata)])
            if(se){
                object[[paste(iType,"se",sep=".")]] <- cbind(object[[paste(iType,"se",sep=".")]][attr(newdata.index,"vectorwise") + (match(times, times.sorted)-1) * NROW(newdata)])
            }
            if(iid){
                iLS.iid <- mapply(slice = attr(newdata.index,"vectorwise"), column = match(times, times.sorted), function(slice,column){
                    object[[paste(iType,"iid",sep=".")]][,column,slice]
                }, SIMPLIFY = FALSE)
                object[[paste(iType,"iid",sep=".")]] <- array(unlist(iLS.iid), dim = c(NROW(object[[paste(iType,"iid",sep=".")]]),1,length(iLS.iid)))
            }
            if(average.iid){
                if(attr(average.iid,"rm.list")){
                    object[[paste(iType,"average","iid",sep=".")]] <- cbind(rowSums(object[[paste(iType,"average","iid",sep=".")]][[1]]))
                }else{
                    object[[paste(iType,"average","iid",sep=".")]] <- lapply(object[[paste(iType,"average","iid",sep=".")]], function(iM){cbind(rowSums(iM))})
                }
            }
        }
        object$diag <- TRUE
    }else{        
        for(iType in type){ ## iType <- "cumhazard"
            object[[iType]] <- object[[iType]][attr(newdata.index,"vectorwise"),,drop=FALSE]
            if(se){
                object[[paste(iType,"se",sep=".")]] <- object[[paste(iType,"se",sep=".")]][attr(newdata.index,"vectorwise"),,drop=FALSE]
            }
            if(iid){
                object[[paste(iType,"iid",sep=".")]] <- object[[paste(iType,"iid",sep=".")]][,,attr(newdata.index,"vectorwise"),drop=FALSE]                    
            }
            if(average.iid && attr(average.iid,"rm.list")){
                object[[paste(iType,"average","iid",sep=".")]] <- object[[paste(iType,"average","iid",sep=".")]][[1]]
            }
            if(needOrder || length(times) != length(times.sorted)){
                col.reorder <- match(times,times.sorted) ## restaure original times and original time ordering (has been reduce to ascending unique times)
                object[[iType]] <- object[[iType]][,col.reorder,drop=FALSE]
                if(iid){
                    object[[paste(iType,"iid",sep=".")]] <- object[[paste(iType,"iid",sep=".")]][,col.reorder,,drop=FALSE]                    
                }
                if(average.iid){
                    if(attr(average.iid,"rm.list")){
                        object[[paste(iType,"average","iid",sep=".")]] <- object[[paste(iType,"average","iid",sep=".")]][,col.reorder,drop=FALSE]
                    }else{
                        object[[paste(iType,"average","iid",sep=".")]] <- lapply(object[[paste(iType,"average","iid",sep=".")]], function(iM){iM[,col.reorder,drop=FALSE]})
                    }
                        
                }
            }
            
        }
    }

    ## ** export
    return(object)
}


##----------------------------------------------------------------------
### compress.R ends here
