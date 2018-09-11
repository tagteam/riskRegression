## ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2018 (17:47) 
## Version: 
## Last-Updated: Sep 11 2018 (12:22) 
##           By: Thomas Alexander Gerds
##     Update #: 1053
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ateRobust (documentation)
#' @title Average Treatment Effects (ATE) for survival outcome (with competing risks) using doubly robust estimating equations 
#' @description Compute the average treatment effect using different estimators:
#' G-formula based on Cox regression, inverse probability of treatment weighting (IPTW)
#' combined with inverse probability of censoring weighting (IPCW), augmented inverse probability weighting (AIPTW, AIPCW).
#' @name ateRobust
#'
#' @param data [data.frame or data.table] Data set in which to evaluate the ATE.
#' @param formula.event [formula] Cox model for the event of interest (outcome model).
#' Typically \code{Surv(time,event)~treatment}.
#' @param formula.censor [formula] Cox model for the censoring (censoring model).
#' Typically \code{Surv(time,event==0)~treatment}.
#' @param formula.treatment [formula] Logistic regression for the treatment (propensity score model).
#' Typically \code{treatment~1}.
#' @param times [numeric] Time point at which to evaluate average treatment effects.
#' @param fitter [character] Routine to fit the Cox regression models.
#' If \code{coxph} use \code{survival::coxph} else use \code{rms::cph}.
#' @param product.limit [logical] If \code{TRUE} the survival is computed using the product limit estimator.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
#' @param cause [numeric/character] The cause of interest. Defaults to the first cause.
#' @param type [character] When set to \code{"survival"} uses a cox model for modeling the survival,
#' otherwise when set to \code{"competing.risks"} uses a Cause Specific Cox model for modeling the absolute risk of the event.
#' @param augment.cens [logical] If \code{TRUE} add an censoring model augmentation term to the estimating equation 
#' @param na.rm [logical] If \code{TRUE} ignore observations whose influence function is NA.
#'
#' @rdname ateRobust
#' @examples
#' # survival outcome, binary treatment X1 
#' library(survival)
#' library(lava)
#' library(data.table)
#' ds=sampleData(101,outcome="survival")
#' ateRobust(data = ds, type = "survival",
#'           formula.event = Surv(time, event) ~ X1+X6,
#'          formula.censor = Surv(time, event==0) ~ X6,
#'          formula.treatment = X1 ~ X6+X2+X7, times = 1,
#'                      product.limit = FALSE)
#' # competing risk outcome, binary treatment X1 
#' library(survival)
#' library(lava)
#' library(data.table)
#' dc=sampleData(101,outcome="competing.risks")
#' ateRobust(data = dc, type = "competing.risks",
#'           formula.event = list(Hist(time, event) ~ X1+X6,Hist(time, event) ~ X6),
#'          formula.censor = Surv(time, event==0) ~ X6,
#'          formula.treatment = X1 ~ X6+X2+X7, times = 1,cause=1,
#'                      product.limit = FALSE)
#'
#'
## * ateRobust (code)
#' @rdname ateRobust
#' @export
ateRobust <- function(data, times, cause, type,
                      formula.event, formula.censor, formula.treatment, 
                      fitter = "coxph", product.limit = FALSE,
                      augment.cens=TRUE,
                      na.rm = FALSE){

    ## ** normalize arguments
    if(!is.data.table(data)){
        data <- data.table::as.data.table(data)
    }else{
        data <- data.table::copy(data)
    }
    n.obs <- NROW(data)
    if(length(times)!=1){
        stop("Argument \'times\' must have length 1 \n")
    }
    treatment <- all.vars(formula.treatment)[1]
    if(treatment %in% names(data) == FALSE){
        stop("Variable \"",treatment,"\" not found in data \n")
    }
    if(!is.factor(data[[treatment]])){
        data[, c(treatment) := as.factor(.SD[[treatment]])]
    }
    level.treatment <- levels(data[[treatment]])
    if(length(level.treatment)!=2){
        stop("only implemented for binary treatment variables \n")
    }
    type <- match.arg(type, c("survival","competing.risks"))
    if(product.limit){
        predictor.cox <- "predictCoxPL"
    }else{
        predictor.cox <- "predictCox"        
    }
    
    ## times since it is an argument that will be used in the data table
    reserved.name <- c("times","time.tau","status.tau","censoring.tau","treatment.bin",
                       "prob.event","prob.event0","prob.event1",
                       "prob.treatment",
                       "weights","prob.censoring","prob.indiv.censoring",
                       "Lterm")
    if(any(reserved.name %in% names(data))){
        txt <- reserved.name[reserved.name %in% names(data)]
        stop("Argument \'data\' should not contain column(s) named \"",paste0(txt, collapse = "\" \""),"\"\n")
    }
    if(type == "survival"){
        varTempo <- SurvResponseVar(formula.event)$status
        if(length(unique(data[[varTempo]]))>2){
            stop("type=\"survival\" can handle at most 2 types of events")
        }
    }
    
    ## ** fit models
    ## event
    if(type == "survival"){
        model.event <- do.call(fitter, args = list(formula = formula.event, data = data, x = TRUE, y = TRUE))
        coxMF <- coxModelFrame(model.event)
        ##  start       stop status strata ...
    }else{
        model.event <- CSC(formula.event, data = data, fitter = fitter, cause = cause, surv.type = "hazard")
        coxMF <- as.data.frame(unclass(model.event$response))
        names(coxMF)[names(coxMF)=="time"] <- "stop"
        ## stop status event ...
    }
    n.censor <- sum( (coxMF$status==0) * (coxMF$stop <= times) )

    
    ## treatment
    model.treatment <- do.call(glm, args = list(formula = formula.treatment, data = data, family = stats::binomial(link = "logit")))       
    data[, c("prob.treatment") := predict(model.treatment, newdata = data, type = "response", se = FALSE)]

    ## censoring
    if(n.censor>0){
        model.censor <- do.call(fitter, args = list(formula = formula.censor, data = data, x = TRUE, y = TRUE)) ## , no.opt = TRUE
    }
    
    ## ** prepare dataset
    ## convert to binary
    data[, c("treatment.bin") := as.numeric(as.factor(.SD[[treatment]]))-1]

    ## counterfactual
    data0 <- copy(data)
    data0[,c(treatment) := factor(level.treatment[1], level.treatment)]

    data1 <- copy(data)
    data1[,c(treatment) := factor(level.treatment[2], level.treatment)]

    ## random variables stopped at times
    if(type == "survival"){
        ## coxMF contains three columns: start, stop, status (0 censored, 1 event), + strata?
        data[,c("time.tau") := pmin(coxMF$stop,times)]
        data[,c("status.tau") := (coxMF$stop<=times)*(coxMF$status==1)]
        ## 0 survival or censored, 1 event
        data[,c("Ncensoring.tau") := as.numeric(((coxMF$stop>=times) + (coxMF$status!=0)) > 0)]
        ## 0 censored, 1 survival or event
    }else if(type=="competing.risks"){
        ## coxMF contains three columns: stop, status (0 censored, 1 any event), event (1 event of interest, 2 competing event, 3 censoring)
        data[,c("time.tau") := pmin(coxMF$stop,times)]
        data[,c("status.tau") := (coxMF$stop<=times)*(coxMF$event==cause)]
        ## 0 survival competing event or censored, 1 event of interest
        data[,c("Ncensoring.tau") := as.numeric(((coxMF$stop>=times) + (coxMF$status!=0)) > 0)]
        ## 0 censored, 1 survival or event (of interest or comepting)
    }
    
    ## ** outcome model: conditional expectation
    ## if(nuisance.iid){
        ## Computation of the influence function (Gformula, AIPW)
        ## additional term: eq:IF-AIPWfull (green terms)
        ## this is sent to predictCox to multiply each individual IF before averaging

        nuisance.iid0 <- TRUE
        attr(nuisance.iid0, "factor") <- list("Gformula" = matrix(1, nrow = n.obs, ncol = 1)
                                              ## "AIPW" = cbind(data[, 1 - (1-.SD$treatment.bin) / (1-.SD$prob.treatment)])
                                              )
        nuisance.iid1 <- TRUE
        attr(nuisance.iid1, "factor") <- list("Gformula" = matrix(1, nrow = n.obs, ncol = 1)
                                              ## "AIPW" = cbind(data[, 1 - .SD$treatment.bin / .SD$prob.treatment])
                                              )
        if(type == "competing.risks"){
            attr(nuisance.iid0, "factor") <- do.call(cbind,attr(nuisance.iid0, "factor"))
            attr(nuisance.iid1, "factor") <- do.call(cbind,attr(nuisance.iid1, "factor"))
        }
    ## }else{
        ## nuisance.iid0 <- FALSE
        ## nuisance.iid1 <- FALSE
    ## }
    
    if(type == "survival"){
        ## Estimation of the survival + IF
        prediction.event <- do.call(predictor.cox, args = list(model.event, newdata = data, times = times, type = "survival"))
        prediction.event0 <- do.call(predictor.cox, args = list(model.event, newdata = data0, times = times, type = "survival", average.iid = nuisance.iid0))
        prediction.event1 <- do.call(predictor.cox, args = list(model.event, newdata = data1, times = times, type = "survival", average.iid = nuisance.iid1))

        ## store results
        data[, c("prob.event") := 1 - prediction.event$survival[,1]]
        data[, c("prob.event0") := 1 - prediction.event0$survival[,1]]
        data[, c("prob.event1") := 1 - prediction.event1$survival[,1]]

        ## if(nuisance.iid){
            iidG.event0 <- -prediction.event0$survival.average.iid[[1]][,1]
            ## iidAIPW.event0 <- -prediction.event0$survival.average.iid[[2]][,1]
            iidG.event1 <- -prediction.event1$survival.average.iid[[1]][,1]
            ## iidAIPW.event1 <- -prediction.event1$survival.average.iid[[2]][,1]
        ## }

    }else if(type=="competing.risks"){

        ## Estimation of the survival + IF
        prediction.event <- predict(model.event, newdata = data, times = times, cause = cause, product.limit = product.limit)
        ## prediction.event0 <- predict(model.event, newdata = data0, times = times, cause = cause, product.limit = product.limit, iid = nuisance.iid)
        prediction.event0 <- predict(model.event, newdata = data0, times = times, cause = cause, product.limit = product.limit, average.iid = nuisance.iid0)
        ## prediction.event1 <- predict(model.event, newdata = data1, times = times, cause = cause, product.limit = product.limit, iid = nuisance.iid)
        prediction.event1 <- predict(model.event, newdata = data1, times = times, cause = cause, product.limit = product.limit, average.iid = nuisance.iid1)

        data[, c("prob.event") := prediction.event$absRisk[,1]]
        data[, c("prob.event0") := prediction.event0$absRisk[,1]]
        data[, c("prob.event1") := prediction.event1$absRisk[,1]]

        ## store results
        ## if(nuisance.iid){
            ## weight0 <- attr(nuisance.iid0,"factor")[,2]
            iidG.event0 <- prediction.event0$absRisk.average.iid[[1]]
            ## range(iidG.event0 - prediction.event0$absRisk.average.iid[[1]])
            ## iidG.event0 <- colMeans(prediction.event0$absRisk.iid[,1,])
            ## iidAIPW.event0 <- prediction.event0$absRisk.average.iid[[2]]
            ## range(iidAIPW.event0 - prediction.event0$absRisk.average.iid[[2]])
            ## iidAIPW.event0 <- apply(prediction.event0$absRisk.iid[,1,],2,function(iCol){mean(weight0*iCol)})

            ## weight1 <- attr(nuisance.iid1, "factor")[,2]
            iidG.event1 <- prediction.event1$absRisk.average.iid[[1]]
            ## range(iidG.event1 - prediction.event1$absRisk.average.iid[[1]])
            ## iidG.event1 <- colMeans(prediction.event1$absRisk.iid[,1,])
            ## iidAIPW.event1 <- prediction.event1$absRisk.average.iid[[2]]
            ## range(iidAIPW.event1 - prediction.event1$absRisk.average.iid[[2]])
            ## iidAIPW.event1 <- apply(prediction.event1$absRisk.iid[,1,],2,function(iCol){mean(weight1*iCol)})
        ## }

    }
    
    ## ** Censoring model: weights
    ## fit model
    if(n.censor==0){
        
        data[,c("prob.censoring") := 1]
        data[,c("prob.indiv.censoring") := 1]
        data[,c("weights") := 1]
        
        }else{
        
            ## survival = P[C>min(T,tau)] = P[Delta(min(T,tau))==1] - ok

        ## stopped at tau
        predTau.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = times, type = "survival"))
        ## predTau.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = times, type = "survival", iid = nuisance.iid))
        data[,c("prob.censoring") := predTau.censor$survival[,1]]
        
        ## at each time
        predIndiv.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), type = "survival", diag = TRUE))
        ## predIndiv.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), type = "survival", diag = TRUE, iid = nuisance.iid))
        data[,c("prob.indiv.censoring") := predIndiv.censor$survival[,1]]
        
        ## store
        data[,c("weights") := as.numeric(NA)]
        data[coxMF$stop<=times, c("weights") := .SD$Ncensoring.tau / .SD$prob.indiv.censoring]
        data[coxMF$stop>times, c("weights") := .SD$Ncensoring.tau / .SD$prob.censoring]

        }
    ## ** Propensity score model: weights
    ## needs to be after censoring to get the weights
    ## if(nuisance.iid){

        factor <- cbind("IPW0" = data[,  .SD$weights  * .SD$status.tau * (1-.SD$treatment.bin) / (1-.SD$prob.treatment)^2],
                        "IPW1" = data[, .SD$weights * .SD$status.tau * .SD$treatment.bin / (.SD$prob.treatment)^2])

        average.iid <- TRUE
        attr(average.iid, "factor") <- factor
        
        prediction.treatment.iid <- attr(.predictGLM(model.treatment, newdata = data, average.iid = average.iid), "iid")
        iidIPW.treatment0 <-  prediction.treatment.iid[,1]
        iidIPW.treatment1 <- -prediction.treatment.iid[,2]
        ## iidAIPW.treatment0 <- -prediction.treatment.iid[,3]
        ## iidAIPW.treatment1 <- prediction.treatment.iid[,4]

    ## ** correction for efficiency
    ## data is updated within the function .calcLterm
    if(augment.cens){
        if(n.censor==0){
            data[,c("Lterm") := 0]
        }else{
            data[,c("Lterm") := .calcLterm(data = data,
                                     n.obs = n.obs,
                                     times = times,
                                     model.censor = model.censor,
                                     model.event = model.event,
                                     type = type,
                                     predictor.cox = predictor.cox,
                                     product.limit = product.limit,
                                     cause = cause)]
        }
    }

    ## ** Compute parameter of interest
    IF <- list()

    ## *** Gformula
    ## eq:IF-Gformula (blue term)
    IF$Gformula <- data[,cbind(
        .SD$prob.event0,
        .SD$prob.event1
    )] / n.obs
    ## if (nuisance.iid){
    nuisanceEvent.Gformula <- cbind(iidG.event0, iidG.event1)
    IF$Gformula <- IF$Gformula + nuisanceEvent.Gformula
    ## }

    ## *** IPW
    ## eq:IF-IPWc (blue term)
    IF$IPTW.IPCW <- data[,cbind(
        .SD$weights * .SD$status.tau * (1-.SD$treatment.bin) / (1-.SD$prob.treatment),
        .SD$weights * .SD$status.tau * (.SD$treatment.bin) / (.SD$prob.treatment)
    )]/ n.obs
    ## colSums(IF$IPTW.IPCW)
    ## colSums(IF$Gformula)
    
    if(augment.cens){
        ## additional term: eq:IF-IPWeff (red term)
        augment.cens.IPW <- data[,cbind(
            .SD$Lterm * (1-.SD$treatment.bin) / (1-.SD$prob.treatment), ## .SD$prob.event
            .SD$Lterm * (.SD$treatment.bin) / (.SD$prob.treatment) ## .SD$prob.event
        )]/ n.obs
        IF$IPTW.AIPCW <- IF$IPTW.IPCW + augment.cens.IPW
        ## colSums(augment.cens.IPW)
        ## sum(data$Lterm)/n.obs
        ## data$prob.treatment
    }
    
    ## additional term: eq:IF-IPWfull (green term)
    nuisanceTreatment.IPW <- cbind(iidIPW.treatment0, iidIPW.treatment1)
    ## nuisanceCensoring.IPW <- cbind(iidIPW.censoring0, iidIPW.censoring1)
    ## colMeans(IPWadd)
    ## apply(IF$IPTW.IPCW + nuisanceTreatment.IPW, 2, sd)
    ## apply(IF$IPTW.IPCW + nuisanceTreatment.IPW + nuisanceCensoring.IPW, 2, sd)

    ## mean(IF$IPTW.IPCW)
    ## mean(nuisanceTreatment.IPW)
    ## mean(nuisanceCensoring.IPW)
    ## IF$IPTW.IPCW <- IF$IPTW.IPCW + nuisanceTreatment.IPW
    ## if(augment.cens){
        ## IF$IPTW.AIPCW.mod <- IF$IPTW.AIPCW + nuisanceTreatment.IPW
    ## }

    ## *** AIPW
    ## additional term:
    AIPWadd <- data[,cbind(
        .SD$prob.event0 * (1-(1-.SD$treatment.bin)/(1-.SD$prob.treatment)),
        .SD$prob.event1 * (1-.SD$treatment.bin/(.SD$prob.treatment))
    )]/ n.obs
    IF$AIPTW.IPCW <- IF$IPTW.IPCW + AIPWadd
    if(augment.cens){
        IF$AIPTW.AIPCW <- IF$IPTW.AIPCW + AIPWadd
    }
    

    ## ** export
    out <- list()

    ## value
    n.estimator <- length(IF)
    name.risk <- paste0("risk.",level.treatment)
    out$ate.value <- matrix(NA, nrow = 3, ncol = n.estimator,
                            dimnames = list(c(name.risk,"ate.diff"),
                                              names(IF)))
    for(iL in 1:n.estimator){ ## iL <- 1
        if(na.rm){
            IF[[iL]] <- IF[[iL]][which(rowSums(is.na(IF[[iL]]))==0),,drop=FALSE]
        }
        out$ate.value[name.risk,iL] <- colSums(IF[[iL]])
        IF[[iL]] <- rowCenter_cpp(IF[[iL]], center = out$ate.value[paste0("risk.",level.treatment),iL]/n.obs)
    }
    out$ate.value["ate.diff",] <- out$ate.value[name.risk[2],] - out$ate.value[name.risk[1],]
    ##    out$ate.value["ate.ratio",] <- out$ate.value[name.risk[2],] / out$ate.value[name.risk[1],]

    ## standard error
    out$ate.se <- do.call(cbind,lapply(IF, function(iIF){ ## iIF <- IF[[1]]
        sqrt(c(colSums(iIF^2), sum((iIF[,2]-iIF[,1])^2)))
    }))
    rownames(out$ate.se) <- rownames(out$ate.value)

    class(out) <- "ateRobust"
    out$augment.cens <- augment.cens
    return(out)
    

}

## * .calcLterm
.calcLterm <- function(data, n.obs, times,
                       model.censor,
                       model.event, type, predictor.cox, product.limit, cause){

    info.censor <- SurvResponseVar(coxFormula(model.censor))
    timeVar.censor <- info.censor$time
    statusVar.censor <- info.censor$status
    
    X.censor <- coxModelFrame(model.censor)
    new.time <- X.censor$stop
    new.status <- X.censor$status
    
    jump.time <- sort(X.censor$stop[X.censor$status == 1]) ##  only select jumps
    jump.time <- jump.time[jump.time <= times] ## before time horizon
    njump <- length(jump.time)

    ## ** compute conditional risk
    riskTau <- matrix(data$prob.event, nrow = n.obs, ncol = njump, byrow = FALSE)
    if(type == "survival"){
        riskTime <- 1 - predictCox(model.event, newdata = data, times = jump.time, type = "survival")$survival
        riskConditional <- (riskTau - riskTime)/(1-riskTime)
    }else if(type == "competing.risks"){
        riskTime <- predict(model.event, newdata = data, times = jump.time,
                            cause = cause, product.limit = product.limit)$absRisk

        survTime <- predictSurv(model.event, newdata = data, times = jump.time, product.limit = product.limit)
        
        riskConditional <- (riskTau - riskTime)/(survTime)

        ## check        
        ## index.test <- 5
        ## GS <- predict(model.event, newdata = data, times = times, landmark = jump.time[index.test],
        ## cause = cause, product.limit = product.limit)$absRisk
        ## range(GS-riskConditional[,index.test])
        
    }
    
    ## ** at risk indicator and counting process for the censoring 
   atRisk <- matrix(NA, nrow = n.obs, ncol = njump)
    dN <- matrix(NA, nrow = n.obs, ncol = njump)
    for(iN in 1:n.obs){
        atRisk[iN,] <- as.numeric(jump.time <= new.time[iN])
        dN[iN,] <- (jump.time == new.time[iN]) * new.status[iN]
    }
    ## if (any(is.na(atRisk)))browser()
    ## ** compensator for the censoring
    dLambda <- predictCox(model.censor, newdata = data, times = jump.time, type = "hazard")$hazard

    ## ** survival for censoring at t-
    if(njump>1){
        survCensoring <- do.call(predictor.cox,
                                 args = list(model.censor, newdata = data, times = jump.time, type = "survival")
                                 )$survival
        survCensoring <- cbind(1,survCensoring[,1:(njump-1)])
    }else{
        survCensoring <- matrix(1, nrow = n.obs, ncol = 1)
    }
    ## ** integral
    
    out <- rowSums(atRisk * riskConditional/survCensoring * (dN-dLambda))
    ## ** export
    return(out)
}
