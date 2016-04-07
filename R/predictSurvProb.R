#' @title Predicting survival, hazard or cumulative hazard from Cox regression model 
#' 
#' @param object The fitted Cox regression model object either obtained with \code{coxph} (survival package) or \code{cph} (rms package).
#' @param newdata A data frame or table containing the values of the predictor variables 
#' @param times Vector of times at which to return the estimated hazards/survival
#' @param type should the hazard or the cumulative hazard be returned
#' @param method.baseHaz The implementation to be used for computing the baseline hazard: "dt" or "cpp"
#' @param keep.strata If \code{TRUE} return an additional column containing the strata level. 
#' 
#' @details 
#' Note: for Cox regression models with time varying covariates this function may not work, and
#'       it does not make sense to ask for survival or cumulative hazard predictions.
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
#' d$time <- round(d$time,1)
#' fit <- coxph(Surv(time,status)~X1 * X2,data=d)
#' # table(duplicated(d$time))
#' ttt <- sample(x = unique(sort(d$time)), size = 10) 
#'
#' res1 <- predictSurv(fit, newdata = d, times = 10)
#' res2 <- predictSurv(fit, newdata = d, times = ttt)
#' res3 <- predictSurv(fit, newdata = d, times = 10, type = "cumHazard")
#' res4 <- predictSurv(fit, newdata = d, times = ttt, type = "cumHazard")
#' 
#' # strata
#' fitS <- coxph(Surv(time,status)~strata(X1)+X2,data=d, ties="breslow")
#' 
#' res1S <- predictSurv(fitS, newdata = d, times = 10)
#' res2S <- predictSurv(fitS, newdata = d, times = ttt)
#' res3S <- predictSurv(fitS, newdata = d, times = 10, type = "cumHazard")
#' res4S <- predictSurv(fitS, newdata = d, times = ttt, type = "cumHazard")
#' @export
predictSurv <- function(object,
                        newdata,
                        times,
                        type = "Survival",
                        method.baseHaz = "cpp",
                        keep.strata = FALSE) {
    
    ## strata?
    if ("cph" %in% class(object))
        strataspecials <- attr(object$terms,"specials")$strat
    else
        strataspecials <- attr(object$terms,"specials")$strata 
    no.strata <- is.null(strataspecials)

    ## linear predictor
    BaseVar <- attr(object$terms, "term.labels")
    
    if(no.strata){ 
        
        stratalabels <- NULL
        strataVar <- NULL
        
    } else {
        
        stratalabels = attr(object$terms,"term.labels")[strataspecials-1]
        BaseVar <- BaseVar[ BaseVar %in% stratalabels == FALSE]
        
        names.newdata <- names(newdata)
        if  ("cph" %in% class(object))
            strataVar <- names.newdata[which(paste0("strat(", names.newdata,")") %in% stratalabels)]
        else
            strataVar <- names.newdata[which(paste0("strata(", names.newdata,")") %in% stratalabels)]
    }
    ## main
    if  ("cph" %in% class(object)){
        if(length(BaseVar) > 0){
            Xb=predict(object, newdata, type = "lp")}
        else{ 
            Xb=rep(0, nrow(newdata)) 
        }
    }
    else{
        if(length(BaseVar) == 0){ 
            Xb=rep(0, nrow(newdata))
        } else if(no.strata){ 
            Xb=survival:::predict.coxph(object, newdata, type = "lp") 
        }else { 
            Xb=rowSums(survival:::predict.coxph(object, newdata = newdata, type = "terms")) 
        }}
    ## baseline hazard
    Lambda0 = baseHazRR(object,
                        method = method.baseHaz,
                        centered = TRUE,
                        maxtime =  max(times))
    ## preparation
    match.arg(type, choices = c("Survival", "cumHazard", "hazard"), several.ok = TRUE)
    n.type <- length(type)
    resPred <- vector(length(type),mode="list")
    names(resPred) <- type
    ## predicted hazard 
    ## time.index <- match(Lambda0$time,times,nomatch=NA)
    if(no.strata){ ## no strata
        if ("hazard" %in% type) 
            resPred[["hazard"]] <- data.table(exp(Xb) %o% Lambda0[, hazard])
        cumHazard <- data.table(exp(Xb) %o% Lambda0[, cumHazard])
        if ("cumHazard" %in% type) 
            resPred[["cumHazard"]] <- cumHazard
        if ("Survival" %in% type)
            resPred[["Survival"]] <- exp(-cumHazard)
    } else{ ## strata
        nPatient <- NROW(newdata)  
        XXstrata <- attr(Xb,"strata")
        levelsStrata <- levels(XXstrata)
        resStrata <- character(nPatient)
    
        ## compute hazard
        for(iterType in 1:n.type){
      
      resPred[[iterType]] <- matrix(nrow = nPatient, ncol = length(times))
      colnames(resPred[[iterType]]) <- paste0("t",times)
      
      for(iterS in levelsStrata){
        
         indexNew <- which(XXstrata == iterS)
         indexHaz <- Lambda0[,.I[strata == iterS]]
        if(length(indexHaz) == 0){
          stop("predictSurv_internal: unknown strata ",iterS," \n",
               "existing strata: \"",paste(levels(Lambda0$strata), collapse = "\" \""),"\"\n")
        }
        
        if(type[iterType] == "Survival"){ # step function interpolation
          time.index <- prodlim::sindex(jump.times = Lambda0[indexHaz,time], eval.times=times)
          res.tempo <- Lambda0[indexHaz[time.index], exp(- exp(Xb[indexNew]) %o% .SD[[1]]), .SDcols = "cumHazard"]
        }else if(type[iterType] == "cumHazard"){ # step function interpolation
          time.index <- prodlim::sindex(jump.times = Lambda0[indexHaz,time], eval.times=times)
          res.tempo <- Lambda0[indexHaz[time.index], exp(Xb[indexNew]) %o% .SD[[1]], .SDcols = "cumHazard"]
        }else if(type[iterType] == "hazard"){ # diract function interpolation
          res.tempo <- matrix(0, nrow = length(indexNew), ncol = length(times))
          time.index <- which(Lambda0[indexHaz,time] %in% times)
          if(length(time.index)>0){
            res.tempo[,which(times %in% Lambda0[indexHaz,time])] <- Lambda0[indexHaz[time.index], exp(Xb[indexNew]) %o% .SD[[1]], .SDcols ="hazard"]
          }
        }
         
          resPred[[iterType]][indexNew,] <- res.tempo
          if(keep.strata == TRUE){
              resStrata[indexNew] <- iterS
          }
      }
      
        if(keep.strata == TRUE){
            resPred[[iterType]][,strata := resStrata]
        }
      
    }
    
    }
    ## export
    names(resPred) <- type
    switch(as.character(n.type),
           "1" = return(resPred[[1]]),
           resPred
           )
  
}
  
  


