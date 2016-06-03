#' @title Fast prediction of survival, hazard and cumulative hazard from Cox regression model 
#'
#' Without strata, hazard, cumulative hazard and survival
#' probabilities are evaluated at the event times of the object.  With
#' strata the hazard, cumulative hazard and survival probabilities are
#' returned at unique event times across strata
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata A data frame or table containing the values of the
#'     predictor variables
#' @param type should the hazard or the cumulative hazard be returned
#' @param keep.strata If \code{TRUE} return an additional column
#'     containing the strata level.
#' @param keep.times  If \code{TRUE} return event times (including censored times) as an additional element 
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds
#'     tag@@biostat.ku.dk
#' @details Note: for Cox regression models with time varying
#'     covariates this function may not work, and it does not make
#'     sense to ask for survival or cumulative hazard predictions.
#' 
#' @return A data table or a matrix containing the predictions for each subject (in rows)
#'         and each time (in columns) and the strata (if requested).
#' 
#' @examples 
#' 
#' library(survival)
#' 
#' set.seed(10)
#' d <- SimSurv(1e2)
#' nd <- SimSurv(10)
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 + X2,data=d)
#' # table(duplicated(d$time))
#' ttt <- sample(x = unique(sort(d$time)), size = 10) 
#'
#' predictSurv(fit, newdata = d)
#' predictSurv(fit, newdata = nd)
#' predictSurv(fit, newdata = d)
#' predictSurv(fit, newdata = d, type = "cumHazard")
#' predictSurv(fit, newdata = d, type = "cumHazard")
#' 
#' # strata
#' fitS = coxph(Surv(time,status)~strata(X1)+X2,
#' data=d, method="breslow",y=TRUE)
#' 
#' predictSurv(fitS, newdata = d)
#' predictSurv(fitS, newdata = nd)
#' all.equal(predictSurv(fitS, newdata = nd,type="survival"),
#'           predictSurvProb(fitS, newdata = nd,times=fitS$y[,"time"]))
#'
#' d$U=sample(letters[1:5],replace=TRUE,size=NROW(d))
#' nd$U=sample(letters[1:5],replace=TRUE,size=NROW(nd))
#' fitSU = coxph(Surv(time,status)~strata(U)+X2+X1,data=d, method="breslow",y=TRUE)
#' cphSU = cph(Surv(time,status)~strat(U)+X2+X1,data=d, method="breslow",y=TRUE,surv=TRUE)
#' predictSurv(fitSU, newdata = nd)
#' predictSurv(cphSU, newdata = nd)
#'
#' d$V=factor(d$X1)
#' nd$V=factor(nd$X1)
#' fitSU2 = coxph(Surv(time,status)~strata(U)+X2+strata(V),data=d, method="breslow",y=TRUE)
#' cphSU2 = cph(Surv(time,status)~strat(U)+X2+strat(X1),data=d, method="breslow",y=TRUE,surv=TRUE)
#' predictSurv(fitSU2, newdata = nd)
#' predictSurv(cphSU2, newdata = nd)
#' 
#' predictSurv(fitS, newdata = d, type = "cumHazard")
#' predictSurv(fitS, newdata = d, type = "cumHazard")
#' 
#' @export
predictSurv <- function(object,
                        newdata,
                        type = "survival",
                        keep.strata = FALSE,
                        keep.times = FALSE){
    ## strata?
    if ("cph" %in% class(object))
        strataspecials <- attr(object$terms,"specials")$strat
    else
        strataspecials <- attr(object$terms,"specials")$strata 
    no.strata <- is.null(strataspecials)
    BaseVar <- attr(object$terms, "term.labels")
    if(no.strata){ 
        stratalabels <- NULL
        strataVar <- NULL
    } else {
        stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
        BaseVar <- BaseVar[ BaseVar %in% stratalabels == FALSE]
        if  ("cph" %in% class(object))
            strataVar <- names(prodlim::parseSpecialNames(stratalabels, special = "strat"))
        else
            strataVar <- names(prodlim::parseSpecialNames(stratalabels, special = "strata"))
    }
    ## linear predictor
    if  ("cph" %in% class(object)){
        if(length(BaseVar) > 0){
            Xb <- stats::predict(object, newdata, type = "lp")}
        else{ 
            Xb <- rep(0, nrow(newdata)) 
        }
    }
    else{
        if(length(BaseVar) == 0){ 
            Xb <- rep(0, nrow(newdata))
        } else if(no.strata){ 
            Xb <- stats::predict(object, newdata, type = "lp") 
        }else { 
            Xb <- rowSums(stats::predict(object, newdata = newdata, type = "terms")) 
        }}
    ## baseline hazard
    Lambda0 = baseHazRR(object, centered = TRUE)
    ## preparation
    type <- match.arg(type, choices = c("survival", "cumHazard", "hazard"), several.ok = TRUE)
    #### main
    out <- list()
    if(no.strata){ ## no strata
        if ("hazard" %in% type){
            hazard <- data.table(exp(Xb) %o% Lambda0[, hazard])
            out <- c(out,list(hazard=hazard))
        }
        if ("cumHazard" %in% type || "survival" %in% type){
            cumHazard <- data.table(exp(Xb) %o% Lambda0[, cumHazard])
            if ("cumHazard" %in% type){out <- c(out,list(cumHazard=cumHazard))}
        }
        if ("survival" %in% type){out <- c(out,list(survival=exp(-cumHazard)))}
    } else{ ## strata
        alltimes <- Lambda0[,sort(unique(time))]
        hazard <- matrix(0, nrow = NROW(newdata), ncol = length(alltimes))
        if(!is.null(attr(Xb,"strata"))){
            newStrata <- attr(Xb,"strata")
            newLevels <- levels(newStrata)
        } else {
            strataFormula <- as.formula(paste0("~ 0 + ",paste(strataVar,collapse = " + ")))
            ## strataFormula <- as.formula(paste0("~ 0 + ",paste(paste("strata(",strataVar,")"),collapse = " + ")))
            newStrata <- interaction(stats::model.frame(strataFormula,newdata), drop = TRUE, sep = ".")
            tempoLevels <- levels(newStrata)
            if (length(grep("=",levels(Lambda0$strata)[1]))>0)
                newLevels <- unlist(lapply(strsplit(tempoLevels, split = ".", fixed = TRUE), function(x){paste(paste(strataVar,x,sep = "="), collapse = ".")}))
            else
                newLevels <- tempoLevels
        }
        originLevels <- levels(Lambda0$strata)
        if(any(newLevels %in% originLevels == FALSE)){
            stop("predictSurv: unknown strata ",paste(newLevels[newLevels %in% originLevels == FALSE], collapse = " ")," \n",
                 "existing strata: \"",paste(originLevels, collapse = "\" \""),"\"\nTo fix this you could try to transform all strata variables into factors.")
        }
        Lambda0$strata <- as.numeric(factor(Lambda0$strata, levels = c(newLevels, setdiff(originLevels,newLevels)))) - 1
        newStrata <- as.numeric(newStrata) - 1
        ## compute hazard in strata
        maxtime.strata <- Lambda0[,.(maxtime=time[.N]),by=strata]
        for(S in 0:(length(newLevels)-1)){
            ## take care of early ending strata
            hazard[newStrata==S,alltimes>maxtime.strata[strata==S,maxtime]] <- NA
            ## find event times within this strata
            time.index.S <- match(Lambda0[strata==S,time],alltimes,nomatch=0L)
            hazard[newStrata == S,time.index.S] <- exp(Xb[newStrata == S]) %o% Lambda0[strata == S, hazard]
        }
        if ("hazard" %in% type){
            if(keep.strata == TRUE){hazard[,strata:=strata]}
            out <- c(out,list(hazard=data.table(hazard)))
        }
        if ("cumHazard" %in% type || "survival" %in% type){
            cumHazard <- rowCumSum(hazard)
        }        
        if ("survival" %in% type){
            survival=data.table(exp(-cumHazard))
            if(keep.strata == TRUE){survival[,strata:=strata]}
            out <- c(out,list(survival=survival))
        }
        if("cumHazard" %in% type){
            if(keep.strata == TRUE){cumHazard[,strata:=strata]}
            out <- c(out,list(cumHazard=data.table(cumHazard)))
        }
    }
    if (keep.times==TRUE) out <- c(out,list(times=Lambda0[,sort(unique(time))]))
    return(out)
}




