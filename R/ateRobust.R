## ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne, Thomas A. Gerds
## Created: jun 27 2018 (17:47) 
## Version: 
## Last-Updated: sep  6 2019 (12:08) 
##           By: Brice Ozenne
##     Update #: 1252
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
#' @description Compute the average treatment effect using different methods:
#' G-formula based on (cause-specific) Cox regression, inverse probability of treatment weighting (IPTW)
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
#' @param se [logical] If \code{TRUE} compute and add the standard errors relative to the G-formula and IPTW method to the output.
#' @param times [numeric] Time point at which to evaluate average treatment effects.
#' @param fitter [character] Routine to fit the Cox regression models.
#' If \code{coxph} use \code{survival::coxph} else use \code{rms::cph}.
#' @param product.limit [logical] If \code{TRUE} the survival is computed using the product limit method.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
#' @param cause [numeric/character] The cause of interest. Defaults to the first cause.
#' @param type [character] When set to \code{"survival"} uses a cox model for modeling the survival,
#' otherwise when set to \code{"competing.risks"} uses a Cause Specific Cox model for modeling the absolute risk of the event.
#' @param augment.cens [logical] If \code{TRUE} add an censoring model augmentation term to the estimating equation 
#' @param na.rm [logical] If \code{TRUE} ignore observations whose influence function is NA.
#'
#' @details The standard errors/confindence intervals/p-values output by ateRobust
#' do not account for the uncertainty related to the estimation of the parameters of the censoring model (only relevant for IPCW/AIPCW estimators).
#' Note that for the AIPTW, this uncertainty is neglectable (i.e. o_p(n^{-1/2})) in correctly specified models.
#'
#' 
#' @rdname ateRobust
#' @examples
#' library(survival)
#' library(lava)
#' library(data.table)
#' library(prodlim)
#' set.seed(10)
#' # survival outcome, binary treatment X1 
#'
#' ds <- sampleData(101,outcome="survival")
#' out <- ateRobust(data = ds, type = "survival",
#'          formula.event = Surv(time, event) ~ X1+X6,
#'          formula.censor = Surv(time, event==0) ~ X6,
#'          formula.treatment = X1 ~ X6+X2+X7, times = 1)
#' out
#' dt.out=as.data.table(out)
#' dt.out
#' 
#' # competing risk outcome, binary treatment X1 
#' dc=sampleData(101,outcome="competing.risks")
#' x=ateRobust(data = dc, type = "competing.risks",
#'           formula.event = list(Hist(time, event) ~ X1+X6,Hist(time, event) ~ X6),
#'          formula.censor = Surv(time, event==0) ~ X6,
#'          formula.treatment = X1 ~ X6+X2+X7, times = 1,cause=1,
#'                      product.limit = FALSE)
#' @seealso \code{\link{ate}} when considering more than 2 treatments or more than 1 timepoint.
## * ateRobust (code)
#' @rdname ateRobust
#' @export
ateRobust <- function(data, times, cause, type,
                      formula.event, formula.censor, formula.treatment, 
                      fitter = "coxph", product.limit = NULL, se = TRUE,
                      augment.cens = TRUE, 
                      na.rm = FALSE){

    ## ** normalize arguments
    data <- data.table::copy(data)
    setDT(data)
    n.obs <- NROW(data)

    if(length(times)!=1){
        stop("Argument \'times\' must have length 1 \n")
    }
    treatment <- all.vars(formula.treatment)[1]
    if(treatment %in% names(data) == FALSE){
        stop("Variable \"",treatment,"\" not found in data \n")
    }
    if(!is.factor(data[[treatment]])){
        data[, c(treatment) := as.factor(.SD[[1]]),.SDcols=treatment]
    }
    level.treatment <- levels(data[[treatment]])
    if(length(level.treatment)!=2){
        stop("only implemented for binary treatment variables \n")
    }
    type <- match.arg(type, c("survival","competing.risks"))
    if(is.null(product.limit)){
        product.limit <- switch(type,
                                "survival" = FALSE,
                                "competing.risks" = TRUE)
    }
    if(product.limit){
        predictor.cox <- "predictCoxPL"
    }else{
        predictor.cox <- "predictCox"        
    }
    
    ## times since it is an argument that will be used in the data table
    reserved.name <- c("times","time.tau","status.tau","censoring.tau","treatment.bin0","treatment.bin1",
                       "prob.event","prob.event0","prob.event1",
                       "prob.treatment0","prob.treatment1", ## "prob.stabilizing0","prob.stabilizing1",
                       "weights","prob.censoring","prob.indiv.censoring",
                       "Lterm")
    if(any(reserved.name %in% names(data))){
        txt <- reserved.name[reserved.name %in% names(data)]
        stop("Argument \'data\' should not contain column(s) named \"",paste0(txt, collapse = "\" \""),"\"\n")
    }

    ## Extract the types of events: check only 2 for survival and output warning if NA
    if(type == "survival"){
        if(!inherits(formula.event, "formula")){
            stop("Argument \'formula.event\' must be a formula when argument \'type\' is set to \"survival\"")
        }
        Ustatus <- unique(data[[SurvResponseVar(formula.event)$status]])        
        if(length(Ustatus)>2){
            stop("The status variable in \'formula.event\' contains ",length(Ustatus)," distinct types of events: ",
                 "\"",paste(Ustatus, collapse = "\" \""),"\"\n",
                 "but ateRobust can handle at most 2 types of events when the argument \'type\' is set to \"survival\" \n")
        }
        if(any(is.na(Ustatus))){
            warning("The status variable in \'formula.event\' contains NA \n")
        }
    }else{
        if(inherits(formula.event, "formula")){
            ls.Ustatus <- unique(data[[SurvResponseVar(formula.event)$status]])
        }else{
            ls.Ustatus <- lapply(formula.event, function(iE){unique(data[[SurvResponseVar(iE)$status]])})
        }
        if(any(is.na(unlist(ls.Ustatus)))){
            warning("One of the status variables in \'formula.event\' contains NA \n")
        }
    }

    
    ## ** fit models (need to be before section "prepare dataset")
    ## *** outcome model
    if(type == "survival"){
        model.event <- do.call(fitter, args = list(formula = formula.event, data = data, x = TRUE, y = TRUE))
        coxMF <- coxModelFrame(model.event)
    }else if(type == "competing.risks"){
        model.event <- CSC(formula.event, data = data, fitter = fitter, cause = cause, surv.type = "hazard")
        coxMF <- as.data.frame(unclass(model.event$response))
        names(coxMF)[names(coxMF)=="time"] <- "stop"
    }
    
    ## *** treatment model
    model.treatment <- do.call(glm, args = list(formula = formula.treatment, data = data, family = stats::binomial(link = "logit")))    

    ## *** censoring model
    n.censor <- sum( (coxMF$status==0) * (coxMF$stop <= times) )
    if(n.censor>0){
        model.censor <- do.call(fitter, args = list(formula = formula.censor, data = data, x = TRUE, y = TRUE)) ## , no.opt = TRUE
    }
    
    ## ** prepare dataset
    ## convert to binary
    data[, c("treatment.bin1") := as.numeric(as.factor(.SD[[1]]))-1,.SDcols=treatment]
    data[, c("treatment.bin0") := 1 - .SD[[1]],.SDcols="treatment.bin1"]

    ## counterfactual
    data0 <- copy(data)
    data0[,c(treatment) := factor(level.treatment[1], level.treatment)]

    data1 <- copy(data)
    data1[,c(treatment) := factor(level.treatment[2], level.treatment)]

    ## random variables stopped at times
    data[,c("time.tau") := pmin(coxMF$stop,times)]
    if(type == "survival"){
        ## coxMF contains three columns: start, stop, status (0 censored, 1 event)
        data[,c("status.tau") := (coxMF$stop<=times)*(coxMF$status==1)]
        ## 0 survival or censored, 1 event
    }else if(type=="competing.risks"){
        ## coxMF contains three columns: stop, status (0 censored, 1 any event), event (1 event of interest, 2 competing event, 3 censoring)
        data[,c("status.tau") := (coxMF$stop<=times)*(coxMF$event==cause)]
        ## 0 survival competing event or censored, 1 event of interest
    }
    
    data[,c("Ncensoring.tau") := (coxMF$stop<times)*(coxMF$status!=0)] ## 1(Tobs<tau,delta!=0)
    ## previously 1(C>min(Tobs,tau)): (coxMF$stop>=times) + (coxMF$stop<times)*(coxMF$status!=0)]

    ## print(data$Ncensoring.tau)

    ## ** Evaluate fitted distributions  
    ## *** outcome model: conditional expectations
    if(type == "survival"){
        ## Estimation of the survival
        prediction.event <- do.call(predictor.cox, args = list(model.event, newdata = data, times = times, type = "survival"))
        prediction.event0 <- do.call(predictor.cox, args = list(model.event, newdata = data0, times = times, type = "survival"))
        prediction.event1 <- do.call(predictor.cox, args = list(model.event, newdata = data1, times = times, type = "survival"))

        ## store results
        data[, c("prob.event") := 1 - prediction.event$survival[,1]]
        data[, c("prob.event0") := 1 - prediction.event0$survival[,1]]
        data[, c("prob.event1") := 1 - prediction.event1$survival[,1]]

    }else if(type=="competing.risks"){
        
        ## Estimation of the cumulative incidence function
        prediction.event <- predict(model.event, newdata = data, times = times, cause = cause, product.limit = product.limit)
        prediction.event0 <- predict(model.event, newdata = data0, times = times, cause = cause, product.limit = product.limit)
        prediction.event1 <- predict(model.event, newdata = data1, times = times, cause = cause, product.limit = product.limit)

        ## store results
        data[, c("prob.event") := prediction.event$absRisk[,1]]
        data[, c("prob.event0") := prediction.event0$absRisk[,1]]
        data[, c("prob.event1") := prediction.event1$absRisk[,1]]
        
    }

    ## *** treatment model: proba
    data[, c("prob.treatment1") := predict(model.treatment, newdata = data, type = "response", se = FALSE)]
    data[, c("prob.treatment0") := 1-.SD[[1]],.SDcols="prob.treatment1"]
    ## if(stabilizing.weight){
    ## model.treatment0 <- do.call(glm, args = list(formula = update(formula.treatment,".~1"), data = data, family = stats::binomial(link = "logit")))
    ## data[, c("stabilizing.treatment0") := predict(model.treatment0, newdata = data, type = "response", se = FALSE)]
    ## data[, c("stabilizing.treatment1") := 1 - .SD$stabilizing.treatment0]
    ## }else{
    ## data[, c("stabilizing.treatment0") := 1]
    ## data[, c("stabilizing.treatment1") := 1]
    ## }

    ## *** censoring model: weights
    ## fit model
    if(n.censor==0){
        
        data[,c("weights") := 1]
        
    }else{
        ## survival = P[C>min(T,tau)] = P[Delta(min(T,tau))==1] - ok
        
        ## stopped at tau
        ## predTau.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = times, type = "survival"))
        ## data[,c("prob.censoring") := predTau.censor$survival[,1]]
        
        ## at each time
        predIndiv.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), type = "survival", diag = TRUE))
        data[,c("weights") := .SD$Ncensoring.tau / predIndiv.censor$survival[,1]]

        ## data[,c("prob.indiv.censoring") := predIndiv.censor$survival[,1]]
        ## store
        ## data[,c("weights") := as.numeric(NA)]
        ## data[coxMF$stop<=times, c("weights") := .SD$Ncensoring.tau / .SD$prob.indiv.censoring]
        ## data[coxMF$stop>times, c("weights") := .SD$Ncensoring.tau / .SD$prob.censoring]
    }

    ## ** Augmentation term for the censoring 
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
    
    ## ** Influence function for the nuisance parameters
    ## needs to be after censoring to get the weights
    if(se>0){
        ## *** outcome model
        nuisance.iid0 <- TRUE
        attr(nuisance.iid0, "factor") <- list("Gformula" = matrix(1, nrow = n.obs, ncol = 1),
                                              "AIPW" = cbind(data[, 1 - .SD$treatment.bin0 / .SD$prob.treatment0])
                                              )
        
        nuisance.iid1 <- TRUE
        attr(nuisance.iid1, "factor") <- list("Gformula" = matrix(1, nrow = n.obs, ncol = 1),
                                              "AIPW" = cbind(data[, 1 - .SD$treatment.bin1 / .SD$prob.treatment1])
                                              )

        if(type == "survival"){
            prediction.event0 <- do.call(predictor.cox, args = list(model.event, newdata = data0, times = times, type = "survival", average.iid = nuisance.iid0))
            prediction.event1 <- do.call(predictor.cox, args = list(model.event, newdata = data1, times = times, type = "survival", average.iid = nuisance.iid1))

            ## there is a minus because predictor.cox return iid(survival) = -iid(risk)
            iidG.event0 <- -prediction.event0$survival.average.iid[[1]][,1]
            iidAIPW.event0 <- -prediction.event0$survival.average.iid[[2]][,1]

            iidG.event1 <- -prediction.event1$survival.average.iid[[1]][,1]
            iidAIPW.event1 <- -prediction.event1$survival.average.iid[[2]][,1]
            
        }else if(type == "competing.risks"){
            attr(nuisance.iid0, "factor") <- do.call(cbind,attr(nuisance.iid0, "factor"))
            attr(nuisance.iid1, "factor") <- do.call(cbind,attr(nuisance.iid1, "factor"))

            prediction.event0 <- predict(model.event, newdata = data0, times = times, cause = cause, product.limit = product.limit, average.iid = nuisance.iid0)
            prediction.event1 <- predict(model.event, newdata = data1, times = times, cause = cause, product.limit = product.limit, average.iid = nuisance.iid1)

            ## weight0 <- attr(nuisance.iid0,"factor")[,2]
            iidG.event0 <- prediction.event0$absRisk.average.iid[[1]]
            ## range(iidG.event0 - prediction.event0$absRisk.average.iid[[1]])
            ## iidG.event0 <- colMeans(prediction.event0$absRisk.iid[,1,])
            iidAIPW.event0 <- prediction.event0$absRisk.average.iid[[2]]
            ## range(iidAIPW.event0 - prediction.event0$absRisk.average.iid[[2]])
            ## iidAIPW.event0 <- apply(prediction.event0$absRisk.iid[,1,],2,function(iCol){mean(weight0*iCol)})

            ## weight1 <- attr(nuisance.iid1, "factor")[,2]
            iidG.event1 <- prediction.event1$absRisk.average.iid[[1]]
            ## range(iidG.event1 - prediction.event1$absRisk.average.iid[[1]])
            ## iidG.event1 <- colMeans(prediction.event1$absRisk.iid[,1,])
            iidAIPW.event1 <- prediction.event1$absRisk.average.iid[[2]]
            ## range(iidAIPW.event1 - prediction.event1$absRisk.average.iid[[2]])
            ## iidAIPW.event1 <- apply(prediction.event1$absRisk.iid[,1,],2,function(iCol){mean(weight1*iCol)})
        }

        ## *** treatment model
        factor <- cbind("IPW0" = data[, .SD$weights * .SD$treatment.bin0 * (.SD$status.tau==1) / .SD$prob.treatment0^2],
                        "IPW1" = data[, .SD$weights * .SD$treatment.bin1 * (.SD$status.tau==1) / .SD$prob.treatment1^2],
                        "AIPW0" = data[, .SD$treatment.bin0 * .SD$prob.event0 / .SD$prob.treatment0^2],
                        "AIPW1" = data[, .SD$treatment.bin1 * .SD$prob.event1 / .SD$prob.treatment1^2])
        
        average.iid <- TRUE
        attr(average.iid, "factor") <- factor
        
        prediction.treatment.iid <- attr(predictGLM(model.treatment, newdata = data, average.iid = average.iid), "iid")
        iidIPW.treatment0 <-  prediction.treatment.iid[,1]
        iidIPW.treatment1 <- -prediction.treatment.iid[,2]
        iidAIPW.treatment0 <- -prediction.treatment.iid[,3]
        iidAIPW.treatment1 <- prediction.treatment.iid[,4]
        
        ## *** censoring model
        ## difficult because we need the iid not for all combinaisons time x obs but only the "diagonal"
        ## nuisance.iid <- TRUE
        ## attr(nuisance.iid, "factor") <- list("IPCW" = data[,.SD$weights^2 * (.SD$status.tau==1) * (.SD$treatment.bin1/ .SD$prob.treatment1 - .SD$treatment.bin0/ .SD$prob.treatment0)])
        ## predIndiv.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), type = "survival", diag = TRUE, average.iid = nuisance.iid1))g
        ## Error in predictCox(list(coefficients = c(X6 = -0.00036957593518861),  : 
        ##   Arguments 'se', 'band', and 'average.iid' must be FALSE when 'diag' is TRUE 
    }
    
    ## ** Compute the (non-centered) influence function
    IF <- list()

    ## *** Gformula
    IF$Gformula <- data[,cbind(
        .SD$prob.event0,
        .SD$prob.event1
    )] / n.obs
    
    ## *** IPW
    IF$IPTW.IPCW <- data[,cbind(
        .SD$weights * .SD$status.tau * .SD$treatment.bin0 / .SD$prob.treatment0,
        .SD$weights * .SD$status.tau * .SD$treatment.bin1 / .SD$prob.treatment1
    )]/ n.obs

    ## *** AIPW
    AIPWadd <- data[,cbind(
        .SD$prob.event0 * (1 - .SD$treatment.bin0  /.SD$prob.treatment0),
        .SD$prob.event1 * (1 - .SD$treatment.bin1  /.SD$prob.treatment1)
    )]/ n.obs
    
    IF$AIPTW.IPCW_knownNuisance <- IF$IPTW.IPCW + AIPWadd

    ## *** iid for the nuisance parameters
    if (se>0){
        nuisanceEvent.Gformula <- cbind(iidG.event0, iidG.event1)
        IF$Gformula <- IF$Gformula + nuisanceEvent.Gformula

        nuisanceTreatment.IPW <- cbind(iidIPW.treatment0, iidIPW.treatment1)
        IF$IPTW.IPCW <- IF$IPTW.IPCW + nuisanceTreatment.IPW

        nuisanceEvent.AIPW <- cbind(iidAIPW.event0, iidAIPW.event1) ## no outcome here so no censoring
        nuisanceTreatment.AIPW <- cbind(iidAIPW.treatment0, iidAIPW.treatment1) ## no outcome here so no censoring

        IF$AIPTW.IPCW_estimatedNuisance <- IF$AIPTW.IPCW_knownNuisance + (nuisanceEvent.AIPW + nuisanceTreatment.IPW + nuisanceTreatment.AIPW)
        ## (nuisanceTreatment.IPW + nuisanceEvent.AIPW + nuisanceTreatment.AIPW)
    }

    ## *** augmentation term for censoring
    if(augment.cens){       
        AUGMENTadd <- data[,cbind(
            .SD$Lterm * .SD$treatment.bin0 / .SD$prob.treatment0, ## .SD$prob.event
            .SD$Lterm * .SD$treatment.bin1 / .SD$prob.treatment1 ## .SD$prob.event
        )]/ n.obs
        
        IF$IPTW.AIPCW <- IF$IPTW.IPCW + AUGMENTadd        
        IF$AIPTW.AIPCW_knownNuisance <- IF$AIPTW.IPCW_knownNuisance + AUGMENTadd
        if(se>0){
            IF$AIPTW.AIPCW_estimatedNuisance <- IF$AIPTW.IPCW_estimatedNuisance + AUGMENTadd
        }
    }
    
    ## ** Estimate ATE with standard error
    out <- list()

    ## *** value
    n.method <- length(IF)
    name.risk <- paste0("risk.",c(0,1))
    out$ate.value <- matrix(NA, nrow = 3, ncol = n.method,
                            dimnames = list(c(name.risk,"ate.diff"),
                                            names(IF)))
    for(iL in 1:n.method){ ## iL <- 1
        if(na.rm){
            IF[[iL]] <- IF[[iL]][which(rowSums(is.na(IF[[iL]]))==0),,drop=FALSE]
        }
        out$ate.value[name.risk,iL] <- colSums(IF[[iL]])
        IF[[iL]] <- rowCenter_cpp(IF[[iL]], center = out$ate.value[paste0("risk.",level.treatment),iL]/n.obs)
    }

    out$ate.value["ate.diff",] <- out$ate.value[name.risk[2],] - out$ate.value[name.risk[1],]
    ##    out$ate.value["ate.ratio",] <- out$ate.value[name.risk[2],] / out$ate.value[name.risk[1],]

    ## *** standard error
    if (se>0){
        out$ate.se <- do.call(cbind,lapply(IF, function(iIF){ ## iIF <- IF[[1]]
            sqrt(c(colSums(iIF^2), sum((iIF[,2]-iIF[,1])^2)))
        }))
        rownames(out$ate.se) <- rownames(out$ate.value)
    }

    ## ** export
    class(out) <- "ateRobust"

    out$se <- se>0
    out$level.treatment <- level.treatment
    out$augment.cens <- augment.cens
    out$product.limit <- product.limit
    out$augment.cens <- augment.cens
    out$weight.treatment <- data[,.SD, .SDcols = c("prob.treatment0","prob.treatment1")]
    if (se>0){out <- confint(out)}

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
        riskTime <- 1 - do.call(predictor.cox,
                                args = list(model.event, newdata = data, times = jump.time, type = "survival"))$survival
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
    
    if(any(riskTau == riskTime)){
        riskConditional[riskTau == riskTime] <- 0
    }

    ## ** at risk indicator and counting process for the censoring 
   atRisk <- matrix(NA, nrow = n.obs, ncol = njump)
    dN <- matrix(NA, nrow = n.obs, ncol = njump)
    for(iN in 1:n.obs){
        atRisk[iN,] <- as.numeric(jump.time <= new.time[iN])
        dN[iN,] <- (jump.time == new.time[iN]) * new.status[iN]
    }
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
