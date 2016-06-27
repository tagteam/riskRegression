#' Fast prediction of survival, hazard and cumulative hazard from Cox regression model 
#'
#' Fast routine to get baseline hazards and subject specific hazards
#' as well as survival probabilities from a coxph or cph object
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata  A data frame or table containing the values of the
#'     predictor variables defining subject specific predictions. Should have
#'     the same structure as the data set used to fit the \code{object}.
#' @param times Time points at which to evaluate the predictions. 
#' @param centered If TRUE remove the centering factor used by \code{coxph}
#'     in the linear predictor.
#' @param maxtime Baseline hazard will be computed for each event before maxtime
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}
#' @param keep.strata Logical. If \code{TRUE} add the (newdata) strata to the output. Only if there any. 
#' @param keep.times Logical. If \code{TRUE} add the evaluation times to the output. 
#' @details Not working with time varying predictor variables or
#'     delayed entry.
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#' @return A list optionally containing the time, the strata (if any), the hazard, the
#'         cumulative hazard and survival probabilities.
#' @examples 
#' library(survival)
#' 
#' set.seed(10)
#' d <- SimSurv(1e2)
#' nd <- SimSurv(10)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d, ties="breslow")
#' # table(duplicated(d$time))
#' 
#' predictCox(fit)
#' predictCox(fit,newdata=nd)
#' cbind(survival::basehaz(fit),predictCox(fit,type="cumHazard"))
#' predictCox(fit, maxtime = 5,keep.times=TRUE)
#' predictCox(fit, maxtime = 1,keep.times=TRUE)
#' 
#' # one strata variable
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' predictCox(fitS)
#' predictCox(fitS, maxtime = 5)
#' predictCox(fitS, maxtime = 5,newdata=nd)
#'
#' # two strata variables
#' set.seed(1)
#' d$U=sample(letters[1:5],replace=TRUE,size=NROW(d))
#' d$V=sample(letters[4:10],replace=TRUE,size=NROW(d))
#' nd$U=sample(letters[1:5],replace=TRUE,size=NROW(nd))
#' nd$V=sample(letters[4:10],replace=TRUE,size=NROW(nd))
#' fit2S <- coxph(Surv(time,status)~X1+strata(U)+strata(V)+X2,data=d, ties="breslow")
#'
#' cbind(survival::basehaz(fit2S),predictCox(fit2S,type="cumHazard"))
#' predictCox(fit2S)
#' predictCox(fitS, maxtime = 5)
#' predictCox(fitS, maxtime = 5,newdata=nd)
#' 
#' 
#' @export
predictCox <- function(object,
                       newdata=NULL,
                       times,
                       centered = TRUE,
                       maxtime = Inf,
                       type=c("hazard","cumHazard","survival"),
                       keep.strata = FALSE,
                       keep.times = FALSE){
    ## extract elements from objects
    xterms <- delete.response(object$terms)
    xvars <- attr(xterms,"term.labels")
    if ("cph" %in% class(object)){
        nPatients <- sum(object$n)
        if(is.null(object$y)){
            stop("Argument \'y\' must be set to TRUE in cph \n")
        }
        strataspecials <- attr(xterms,"specials")$strat
        stratavars <- xvars[strataspecials]
        is.strata <- length(strataspecials)>0
        if(is.strata){
            if (length(xvars)>length(strataspecials)) 
                sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
            else 
                sterms <- xterms
            stratavars <- xvars[strataspecials]
            strataF <- object$Strata
            ## stratalevels <- levels(factor(object$strata))
            ## names(stratalevels) <- stratavars
        }else{
            strataF <- factor("1")
        }
    } else if ("coxph" %in% class(object)){
        nPatients <- object$n
        strataspecials <- attr(xterms,"specials")$strata
        stratavars <- xvars[strataspecials]
        is.strata <- length(strataspecials)>0
        if(is.strata){
            if (length(xvars)>length(strataspecials)) 
                sterms <- stats::drop.terms(xterms,(1:length(xvars))[-strataspecials])
            else 
                sterms <- xterms
            stratalevels <- object$xlevels[stratavars]
            strataF <- interaction(stats::model.frame(object)[,stratavars], drop = TRUE, sep = ", ", lex.order = TRUE) 
        }else{
            strataF <- factor("1")
        }
    } else {
        stop("Only implemented for \"coxph\" and \"cph\" objects \n")
    }
    ## method
    if(object$method == "exact"){
        stop("Prediction with exact correction for ties is not implemented \n")
    }
    if (is.na(maxtime)) maxtime <- Inf
    ## main
    levelsStrata <- levels(strataF)
    nStrata <- length(levelsStrata)
    ytimes <- object$y[,"time"]
    status <- object$y[,"status"]
    Lambda0 <- BaseHazStrata_cpp(alltimes = ytimes,
                                 status = status,
                                 Xb = if(centered == FALSE){object$linear.predictors + sum(object$means*stats::coef(object))}else{object$linear.predictors},
                                 strata = as.numeric(strataF) - 1,
                                 nPatients = nPatients,
                                 nStrata = nStrata,
                                 maxtime = maxtime,
                                 cause = 1,
                                 Efron = (object$method == "efron"))
    if (is.strata)
        etimes <- sort(unique(Lambda0$time))
    else
        etimes <- Lambda0$time
    if (is.strata == TRUE){ ## rename the strata value with the correct levels
        Lambda0$strata <- factor(Lambda0$strata, levels = 0:(nStrata-1), labels = levelsStrata)
    }
    if (is.null(newdata)){
        hazard <- Lambda0$hazard
        cumHazard <- Lambda0$cumHazard
        survival <- exp(-Lambda0$cumHazard)
        
        if(is.strata && keep.strata==TRUE){
          newstrata <- Lambda0$strata
        }
        
    } else{
        ## linear predictor
        if  ("cph" %in% class(object)){
            if(length(xvars) > length(stratavars)){
                Xb <- stats::predict(object, newdata, type = "lp")}
            else{ 
                Xb <- rep(0, NROW(newdata)) 
            }
        } else{ ## coxph
            if(length(xvars) == length(stratavars)){ 
                Xb <- rep(0, NROW(newdata))
            } else if(is.strata){ 
                Xb <- rowSums(stats::predict(object, newdata = newdata, type = "terms")) 
            }else { 
                Xb <- stats::predict(object, newdata, type = "lp") 
            }
        }
        ## subject specific hazard
        if (is.strata==FALSE){
            if ("hazard" %in% type){
                hazard <- exp(Xb) %o% Lambda0$hazard
            }
            if ("cumHazard" %in% type || "survival" %in% type){
                cumHazard <- exp(Xb) %o% Lambda0$cumHazard
            }
            if ("survival" %in% type){
                survival <- exp(-cumHazard)
            }
        }else{
            hazard <- matrix(0, nrow = NROW(newdata), ncol = length(etimes))
            if ("cph" %in% class(object)){
                tmp <- model.frame(sterms,newdata)
                colnames(tmp) <- names(prodlim::parseSpecialNames(names(tmp),"strat"))
                tmp <- data.frame(lapply(1:NCOL(tmp),function(j){factor(paste0(names(tmp)[j],"=",tmp[,j,drop=TRUE]))}))
                newstrata <- apply(tmp,1,paste,collapse=".")
                ## newstrata <- prodlim::model.design(sterms,data=newdata,specialsFactor=TRUE)$strat[[1]]                
            }else{
                newstrata <- prodlim::model.design(sterms,data=newdata,xlev=stratalevels,specialsFactor=TRUE)$strata[[1]]
            }
            ## loop across strata
            for(S in unique(newstrata)){
                id.S <- Lambda0$strata==S
                newid.S <- newstrata==S
                hazard.S <- Lambda0$hazard[id.S]
                maxtime <- max(Lambda0$time[id.S])
                ## take care of early ending strata
                hazard[newid.S,etimes>maxtime] <- NA
                ## find event times within this strata
                time.index.S <- match(Lambda0$time[id.S],etimes,nomatch=0L)
                hazard[newid.S,time.index.S] <- exp(Xb[newid.S]) %o% Lambda0$hazard[id.S]
            }
            if ("cumHazard" %in% type || "survival" %in% type){
                cumHazard <- rowCumSum(hazard)
            }
            if ("survival" %in% type){
                survival <- exp(-cumHazard)
            }
        }
    }
    out <- list()
    if (missing(times)){
        if ("hazard" %in% type) out <- c(out,list(hazard=hazard))
        if ("cumHazard" %in% type) out <- c(out,list(cumHazard=cumHazard))
        if ("survival" %in% type) out <- c(out,list(survival=survival))
        if (keep.times==TRUE) out <- c(out,list(times=Lambda0$time))
    }else{
      if(any(times < 0)){
        stop("Time points at which to evaluate the predictions must be positive \n")
      }
      
        if ("hazard" %in% type) {
            hits <- etimes%in%times
            if (sum(hits)<length(times)){
                if (sum(hits==0)) {
                    out <- c(out,list(hazard=matrix(0,ncol=length(times),nrow=NROW(newdata))))
                }else{
                    hh <- matrix(0,ncol=length(times),nrow=NROW(newdata))
                    hh[,hits] <- hazard[times%in%etimes,,drop=FALSE]
                    out <- c(out,list(hazard=hh))
                }
            }else{
                out <- c(out,list(hazard=hazard[,hits,drop=FALSE]))
            }
        }
        if ("cumHazard" %in% type || "survival" %in% type){
            tindex <- prodlim::sindex(jump.times=etimes,eval.times=times)
        }
        if ("cumHazard" %in% type) out <- c(out,list(cumHazard=cbind(0,cumHazard)[,tindex+1]))
        if ("survival" %in% type) out <- c(out,list(survival=cbind(1,survival)[,tindex+1]))
        if (keep.times==TRUE) out <- c(out,list(times=times))
    }
    if (is.strata && keep.strata==TRUE) out <- c(out,list(strata=newstrata))
    out
}

