#' Risk Regression
#' Fits a regression model for the risk of an event -- allowing for competing
#' risks.
#' 
#' This is the twin-sister of the function \code{comp.risk} from the timereg package.
#' 
#' @aliases riskRegression ARR LRR
#' @param formula Formula where the left hand side specifies the event
#' history event.history and the right hand side the linear predictor.  See
#' examples.
#' @param data The data for fitting the model in which includes all
#' the variables included in formula.
#' @param times Vector of times. For each time point in \code{times}
#' estimate the baseline risk and the timevarying coefficients.
#' @param link \code{"relative"} for the absolute risk regression
#' model.  \code{"logistic"} for the logistic risk regression model.
#' \code{"prop"} for the Fine-Gray regression model.
#' @param cause The cause of interest.
#' @param confint If \code{TRUE} return the iid decomposition, that
#' can be used to construct confidence bands for predictions.
#' @param cens.model Specified the model for the (conditional)
#' censoring distribution used for deriving weights (ICPW). Defaults
#' to "KM" (the Kaplan-Meier method ignoring covariates) alternatively
#' it may be "Cox" (Cox regression).
#' @param cens.formula Right hand side of the formula used for fitting
#' the censoring model.  If not specified the right hand side of
#' \code{formula} is used.
#' @param numSimu Number of simulations in resampling.
#' @param maxiter Maximal number of iterations.
#' @param silent Set this to \code{0} to see some system messages
#' during fitting.
#' @param convLevel Integer between 1 and 10, specifying the
#' convergence criterion for the fitting algorithm as
#' \code{10^{-convLevel}}.
#' @param conservative If \code{TRUE} use variance formula that ignores the contribution by the estimate of the inverse of the probability of censoring weights 
#' @param ... Not used.
#' @author Thomas H. Scheike \email{ts@@biostat.ku.dk}
#' Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @references Gerds, TA and Scheike, T and Andersen, PK (2011) Absolute risk
#' regression for competing risks: interpretation, link functions and
#' prediction Research report 11/8. Department of Biostatistics, University of
#' Copenhagen
#' 
#' Scheike, Zhang and Gerds (2008), Predicting cumulativeincidence probability
#' by direct binomial regression,Biometrika, 95, 205-220.
#' 
#' Scheike and Zhang (2007), Flexible competing risks regression modelling and
#' goodness of fit, LIDA, 14, 464-483.
#' 
#' Martinussen and Scheike (2006), Dynamic regression models for survival data,
#' Springer.
##' @examples
##' 
##' 
##' data(Melanoma,package="riskRegression")
##' ## tumor thickness on the log-scale
##' Melanoma$logthick <- log(Melanoma$thick)
##' 
##' # Single binary factor
##' 
##' ## absolute risk regression
##' library(prodlim)
##' fit.arr <- ARR(Hist(time,status)~sex,data=Melanoma,cause=1)
##' print(fit.arr)
##' # show predicted cumulative incidences
##' plot(fit.arr,col=3:4,newdata=data.frame(sex=c("Female","Male")))
##' 
##' ## compare with non-parametric Aalen-Johansen estimate
##' library(prodlim)
##' fit.aj <- prodlim(Hist(time,status)~sex,data=Melanoma)
##' plot(fit.aj,confint=FALSE)
##' plot(fit.arr,add=TRUE,col=3:4,newdata=data.frame(sex=c("Female","Male")))
##' 
##' ## with time-dependent effect
##' fit.tarr <- ARR(Hist(time,status)~strata(sex),data=Melanoma,cause=1)
##' plot(fit.tarr,newdata=data.frame(sex=c("Female","Male")))
##' 
##' ## logistic risk regression
##' fit.lrr <- LRR(Hist(time,status)~sex,data=Melanoma,cause=1)
##' summary(fit.lrr)
##' 
##' 
##' # Single continuous factor
##' 
##' ## tumor thickness on the log-scale
##' Melanoma$logthick <- log(Melanoma$thick)
##' 
##' ## absolute risk regression 
##' fit2.arr <- ARR(Hist(time,status)~logthick,data=Melanoma,cause=1)
##' print(fit2.arr)
##' # show predicted cumulative incidences
##' plot(fit2.arr,col=1:5,newdata=data.frame(logthick=quantile(Melanoma$logthick)))
##' 
##' ## comparison with nearest neighbor non-parametric Aalen-Johansen estimate
##' library(prodlim)
##' fit2.aj <- prodlim(Hist(time,status)~logthick,data=Melanoma)
##' plot(fit2.aj,confint=FALSE,newdata=data.frame(logthick=quantile(Melanoma$logthick)))
##' plot(fit2.arr,add=TRUE,col=1:5,lty=3,newdata=data.frame(logthick=quantile(Melanoma$logthick)))
##' 
##' ## logistic risk regression
##' fit2.lrr <- LRR(Hist(time,status)~logthick,data=Melanoma,cause=1)
##' summary(fit2.lrr)
##' 
##' ## change model for censoring weights
##' library(rms)
##' fit2a.lrr <- LRR(Hist(time,status)~logthick,
##'                  data=Melanoma,
##'                  cause=1,
##'                  cens.model="cox",
##'                  cens.formula=~sex+epicel+ulcer+age+logthick)
##' summary(fit2a.lrr)
##' 
##' ##  compare prediction errors
##' \dontrun{
##'     library(pec)
##'     plot(pec(list(ARR=fit2.arr,AJ=fit2.aj,LRR=fit2.lrr),data=Melanoma,maxtime=3000))
##' }
##' 
##' 
##' # multiple regression
##' library(pec)
##' library(riskRegression)
##' library(prodlim)
##' # absolute risk model
##' multi.arr <- ARR(Hist(time,status)~logthick+sex+age+ulcer,data=Melanoma,cause=1)
##' 
##' # stratified model allowing different baseline risk for the two gender
##' multi.arr <- ARR(Hist(time,status)~thick+strata(sex)+age+ulcer,data=Melanoma,cause=1)
##' 
##' # stratify by a continuous variable: strata(age)
##' multi.arr <- ARR(Hist(time,status)~tp(thick,power=0)+strata(age)+sex+ulcer,
##'                  data=Melanoma,
##'                  cause=1)
##' 
##' fit.arr2a <- ARR(Hist(time,status)~tp(thick,power=1),data=Melanoma,cause=1)
##' summary(fit.arr2a)
##' fit.arr2b <- ARR(Hist(time,status)~timevar(thick),data=Melanoma,cause=1)
##' summary(fit.arr2b)
##' 
##' ## logistic risk model
##' fit.lrr <- LRR(Hist(time,status)~thick,data=Melanoma,cause=1)
##' summary(fit.lrr)
##' 
##' 
##' 
##' 
##' 
##' ## nearest neighbor non-parametric Aalen-Johansen estimate
##' library(prodlim)
##' fit.aj <- prodlim(Hist(time,status)~thick,data=Melanoma)
##' plot(fit.aj,confint=FALSE)
##' 
##' \dontrun{
##'     # prediction error
##'     library(pec)
##'     x <- pec(list(fit.arr2a,fit.arr2b,fit.lrr),
##'              data=Melanoma,
##'              formula=Hist(time,status)~1,
##'              cause=1,
##'              B=10,
##'              splitMethod="none")
##' }
##'
##' 
#' @keywords survival
#' @export
riskRegression <- function(formula,
                           data,
                           times,
                           link="relative",
                           cause,
                           confint=TRUE,
                           cens.model,
                           cens.formula,
                           numSimu=0,
                           maxiter=50,
                           silent=1,
                           convLevel=6,
                           conservative=TRUE,
                           ...){
    # {{{ preliminaries
    weighted=0
    detail=0
    stopifnot(is.numeric(maxiter)&&maxiter>0&&(round(maxiter)==maxiter))
    stopifnot(silent %in% c(0,1))
    stopifnot(convLevel %in% 1:10)
    conv <- 10^{-convLevel}
    # trans=1 P_1=1-exp(- ( x' b(b)+ z' gam t^time.pow) ), 
    # trans=2 P_1=1-exp(-exp(x a(t)+ z` b )
    # trans=not done P_1=1-exp(-x a(t) exp(z` b )) is not good numerically
    # trans=3 P_1=exp(-exp(x a(t)+ z` b )
    ##   trans <- switch(link,"log"=
    trans <- switch(link,
                    "additive"=1, # 
                    "prop"=2,     # Proportional hazards (Cox, FG)
                    "logistic"=3, # Logistic absolute risks 
                    "relative"=4) # Relative absolute risks
    if (numSimu==0) sim <- 0 else sim <- 1
    # }}}
    # {{{ read the data and the design
    call <- match.call()
    EHF <- prodlim::EventHistory.frame(formula,
                                       data,
                                       specials=c("timevar","strata","prop","const","tp"),
                                       stripSpecials=c("timevar","prop"),
                                       stripArguments=list("prop"=list("power"=0),"timevar"=list("test"=0)),
                                       stripAlias=list("timevar"=c("strata"),"prop"=c("tp","const")),
                                       stripUnspecials="prop",
                                       specialsDesign=TRUE,
                                       dropIntercept=TRUE)
    ## Terms <- attr(EHF,"Terms")
    Terms <- terms(formula)
    event.history <- EHF$event.history
    n <- NROW(event.history)
    # }}}
    # {{{ time-constant design matrix Z
    ## Specials <- attr(Terms,"stripped.specials")
    ## strippedArguments <- attr(Terms,"stripped.arguments")
    Z <- EHF$prop
    if (!is.null(Z)){
        power.args <- attr(EHF$prop,"arguments.terms")$power
        stopifnot(all(match(colnames(Z),names(power.args),nomatch=0)))
        timeconst.power <- as.numeric(power.args[colnames(Z)])
        factorLevelsZ <- attr(EHF$prop,"levels")
        refLevelsZ <- lapply(factorLevelsZ,function(x)x[1])
        colnamesZ <- colnames(Z)
        dimZ <- NCOL(Z)
        fixed <- 1
        names(timeconst.power) <- colnames(Z)
        stopifnot(length(timeconst.power)==dimZ)
        if (!(all(timeconst.power %in% 0:2)))
            stop("Only powers of time in 0,1,2 can be multipled to constant covariates")
    } else{
        Z <- matrix(0,n,1)
        dimZ <- 1
        colnamesZ <- NULL
        fixed <- 0
        factorLevelsZ <- NULL
        refLevelsZ <- NULL
        timeconst.power <- NULL
    }
    # }}}
    # {{{ time-varying design matrix X
    X <- EHF$timevar
    if (!is.null(X)){
        test.args <- attr(EHF$timevar,"arguments.terms")$test
        stopifnot(all(match(colnames(X),names(test.args),nomatch=0)))
        timevar.test <- as.numeric(test.args[colnames(X)])
        names(timevar.test) <- colnames(X)
        factorLevelsX <- attr(EHF$timevar,"levels")
        refLevelsX <- lapply(factorLevelsX,function(x)x[1])
        ## intercept 
        X <- cbind("Intercept"=rep(1,n),X)
        timevar.test <- c(0,as.numeric(timevar.test))
        dimX <- NCOL(X)
    } else{
        ## X <- matrix(0,n,1)
        dimX <- 1
        colnamesX <- NULL
        factorLevelsX <- NULL
        refLevelsX <- NULL
        timevar.test <- 0
        ## intercept 
        X <- cbind("Intercept"=rep(1,n))
    }
    stopifnot(length(timevar.test)==dimX)
    colnamesX <- colnames(X)
    if (!(all(timevar.test %in% 0:2)))
        stop("Time power tests only available for powers 0,1,2")
    theData <- data.frame(cbind(unclass(event.history),do.call("cbind",EHF[-1])))
    # }}}
    # {{{ event.history and order the data
    delayed <- !(is.null(attr(event.history,"entry.type"))) && !(attr(EHF$event.history,"entry.type")=="")
    if (delayed){
        stop("Delayed entry is not (not yet) supported.")
    }
    cens.code <- attr(event.history,"cens.code")
    model.type <- attr(event.history,"model")
    states <- prodlim::getStates(event.history)
    stopifnot(model.type %in% c("survival","competing.risks"))
    cens.type <- attr(event.history,"cens.type")
    stopifnot(cens.type %in% c("rightCensored","uncensored"))
    neworder <- order(event.history[,"time"],-event.history[,"status"])
    event.history <- event.history[neworder,,drop=FALSE]
    Z <- Z[neworder,,drop=FALSE]
    X <- X[neworder,,drop=FALSE]
    theData <- theData[neworder,]
    if (model.type!="survival" && !("event" %in% colnames(event.history)))
        warning("Only one cause of failure found in data.")
    eventtime <- as.vector(event.history[,"time"])
    time  <- numeric(length(eventtime))
    if (model.type %in% c("competing.risks","survival")){
        if (cens.type %in% c("rightCensored","uncensored")){
            delta <- as.vector(event.history[,"status"])
            if (model.type=="competing.risks"){
                if (missing(cause)) {
                    cause <- states[1]
                    warning(paste("Argument cause is missing, analysing cause 1: ",cause,". Other causes are:",paste(states[-1],collapse=","),sep=""))
                }else {
                    if (!(cause %in% states)) stop(paste("Cause",cause," is not among the causes in data; these are:",paste(states,collapse=",")))
                    cause <- match(cause,states,nomatch=0)
                }
                ## event is 1 if the event of interest occured and 0 otherwise
                event <-  event.history[,"event"] == cause
                if (sum(event)==0) stop(paste("No events of type:", cause, "in data."))
            } else{
                event <- delta
            }
        }else{
            stop("Works only for right-censored data")
        }
    }else{stop("Response is neither competing risks nor survival.")}
    # }}}
    # {{{ cluster variable
    clusters <- EHF$cluster
    if(is.null(clusters)){
        clusters  <-  0:(NROW(X) - 1)
        antclust  <-  NROW(X)
    } else {
        clusters  <-  as.integer(factor(clusters))-1
        antclust  <-  length(unique(clusters))
    }
    # }}}
    # {{{ time points for timevarametric components
    if (missing(times) || length(times)==0) {
        times <- sort(unique(eventtime[event]))
        ## times <- times[-c(1:5)] 
    }
    else{
        times <- sort(unique(times))
    }
    ntimes <- length(times)
    if (ntimes>1) silent <- c(silent,rep(0,ntimes-1))
    # }}}
    # {{{ estimate ipcw
    if (missing(cens.formula)){
        cterms <- prodlim::strip.terms(terms(update(formula,NULL~.),specials=c("tp","timevar","strata")),
                                       specials=c("tp","timevar","strata"),
                                       arguments=list("tp"="power","timevar"="test","strata"="test"))
        cens.formula <- formula(cterms)
    }
    else
        cens.formula <- update(cens.formula,NULL~.)
    if (missing(cens.model)) cens.model <- "KM"
    imodel <- switch(tolower(cens.model),"km"="marginal","cox"="cox","forest"="forest","aalen"="aalen","uncensored"="none")
    if (imodel == "cox")
        iFormula <- update(cens.formula,"survival::Surv(time,status)~.")
    else
        iFormula <- update(cens.formula,"Surv(time,status)~.")
    if (imodel=="marginal")
        iData <- data.frame(event.history)
    else{
        iData <- cbind(event.history,get_all_vars(cens.formula,
                                                  data)[neworder,,drop=FALSE])
    }
    ## browser()
    stopifnot(NROW(iData)==NROW(event.history))
    Gcx <- subjectWeights(formula=iFormula,
                          data=iData,
                          method=cens.model,
                          lag=1)$weights
    ## }
    # }}}
    # {{{ prepare fitting
    if (confint == TRUE){
        biid  <-  double(ntimes* antclust * dimX)
        gamiid <-  double(antclust *dimZ)
    } else {
        gamiid  <-  biid  <-  NULL
    }
    # }}}
    # {{{ C does the hard work

    line <- ifelse(trans==1,1,0)
    # if line=1 then test "b(t) = gamma t"
    # if line=0 then test "b(t) = gamma "
    ordertime <- order(eventtime,-delta)
    out <- .C("itfit",
              times=as.double(times),
              ntimes=as.integer(ntimes),
              eventtime=as.double(eventtime),
              cens.code=as.integer(cens.code),
              event=as.integer(event),
              Gcx=as.double(Gcx),
              X=as.double(X),
              n=as.integer(n),
              dimX=as.integer(dimX),
              maxiter=as.integer(maxiter),
              betas=double(dimX),
              score=double(ntimes*(dimX+1)),
              double(dimX*dimX),
              est=double(ntimes*(dimX+1)),
              var=double(ntimes*(dimX+1)),
              sim=as.integer(sim),
              numSimu=as.integer(numSimu),
              test=double(numSimu*3*dimX),
              testOBS=double(3*dimX),
              Ut=double(ntimes*(dimX+1)),
              simUt=double(ntimes*50*dimX),
              weighted=as.integer(weighted),
              gamma=double(dimZ),
              var.gamma=double(dimZ*dimZ),
              fixed=as.integer(fixed),
              Z=as.double(Z),
              dimZ=as.integer(dimZ),
              trans=as.integer(trans),
              gamma2=double(dimX),
              cause=as.integer(1),
              line=as.integer(line),
              detail=as.integer(detail),
              biid=as.double(biid),
              gamiid=as.double(gamiid),
              resample.iid=as.integer(as.numeric(confint)),
              timeconst.power=as.double(timeconst.power),
              clusters=as.integer(clusters),
              antclust=as.integer(antclust),
              timevar.test=as.double(timevar.test),
              silent=as.integer(silent),
              conv=as.double(conv),
              # new since 10 Sep 2014 (07:07)
              weights=as.double(rep(1,n)),
              entry=as.double(rep(0,length(eventtime))),
              trunkp=as.double(rep(1,n)),
              estimator=as.integer(1),
              fixgamma=as.integer(0),
              stratum=as.integer(0),
              ordertime=as.integer(ordertime-1),
              ssf=double(0),
              conservative=as.integer(conservative),
              PACKAGE="riskRegression")

    # }}}
    # {{{ prepare the output

    if (fixed==1){
        timeConstantCoef <- out$gamma
        names(timeConstantCoef) <- colnamesZ
        timeConstantVar <- matrix(out$var.gamma,dimZ,dimZ,dimnames=list(colnamesZ,colnamesZ))
    }
    else{
        timeConstantCoef <- NULL
        timeConstantVar <- NULL
    }
    ## FIXME out$est should not include time
    timeVaryingCoef <- matrix(out$est,ntimes,dimX+1,dimnames=list(NULL,c("time",colnamesX)))
    timeVaryingVar <- matrix(out$var,ntimes,dimX+1,dimnames=list(NULL,c("time",colnamesX)))
    score <- matrix(out$score,ntimes,dimX+1,dimnames=list(NULL,c("time",colnamesX)))
    if (is.na(sum(score[,-1])))
        score <- NA
    else 
        if (sum(score[,-1])<0.00001)
            score <- sum(score[,-1])
    ## time power test
    testedSlope <- out$gamma2
    names(testedSlope) <- colnamesX
    if (confint==1)  {
        biid <- matrix(out$biid,ntimes,antclust*dimX)
        if (fixed==1) gamiid <- matrix(out$gamiid,antclust,dimZ) else gamiid <- NULL
        B.iid <- list()
        for (i in (0:(antclust-1))*dimX) {
            B.iid[[i/dimX+1]] <- matrix(biid[,i+(1:dimX)],ncol=dimX)
            colnames(B.iid[[i/dimX+1]]) <- colnamesX
        }
        if (fixed==1) colnames(gamiid) <- colnamesZ
    } else B.iid <- gamiid <- NULL
    if (sim==1) {
        simUt <- matrix(out$simUt,ntimes,50*dimX)
        UIt <- list()
        for (i in (0:49)*dimX) UIt[[i/dimX+1]] <- as.matrix(simUt[,i+(1:dimX)])
        Ut <- matrix(out$Ut,ntimes,dimX+1)
        colnames(Ut) <-  c("time",colnamesX)
        test <- matrix(out$test,numSimu,3*dimX)
        testOBS <- out$testOBS
        supUtOBS <- apply(abs(Ut[,-1,drop=FALSE]),2,max)
        # {{{ confidence bands
    
        percen<-function(x,per){
            n<-length(x)
            tag<-round(n*per)+1
            out<-sort(x)[tag]
            return(out)
        }
        unifCI <- do.call("cbind",lapply(1:dimX,function(i)percen(test[,i],0.95)))
        colnames(unifCI) <- colnamesX
        # }}}
        # {{{ Significance test
        posSig <- 1:dimX
        sTestSig <- testOBS[posSig]
        names(sTestSig) <- colnamesX
        pval<-function(simt,Otest)
            {
                simt<-sort(simt)
                p<-sum(Otest<simt)/length(simt)
                return(p)
            }
    
        pTestSig <- sapply(posSig,function(i){pval(test[,i],testOBS[i])})
        names(pTestSig) <- colnamesX
        timeVarSignifTest <- list(Z=sTestSig,pValue=pTestSig)
        # }}}
        # {{{ Kolmogoroff-Smirnoff test
        posKS <- (dimX+1):(2*dimX)
        sTestKS <- testOBS[posKS]
        names(sTestKS) <- colnamesX
        pTestKS <- sapply(posKS,function(i){pval(test[,i],testOBS[i])})
        names(pTestKS) <- colnamesX
        timeVarKolmSmirTest <- list(Z=sTestKS,pValue=pTestKS)
        # }}}
        # {{{ Kramer-von-Mises test
        posKvM <- (2*dimX+1):(3*dimX)
        sTestKvM <- testOBS[posKvM]
        names(sTestKvM) <- colnamesX
        pTestKvM <- sapply(posKvM,function(i){pval(test[,i],testOBS[i])})
        names(pTestKvM) <- colnamesX
        timeVarKramvMisTest <- list(Z=sTestKvM,pValue=pTestKvM)
        # }}}
    }
    # }}}
    # {{{ return results
    timeConstantEffects <- list(coef=timeConstantCoef,var=timeConstantVar)
    class(timeConstantEffects) <- "timeConstantEffects"
    timeVaryingEffects <- list(coef=timeVaryingCoef,
                               var=timeVaryingVar)
    class(timeVaryingEffects) <- "timeVaryingEffects"
    out <- list(call=call,
                response=event.history,
                design=list(Terms=Terms,
                    const=colnamesZ,
                    timevar=colnamesX,
                    timepower=timeconst.power),
                link=link,
                time=times,
                timeConstantEffects=timeConstantEffects,
                timeconst.power=timeconst.power,
                timeVaryingEffects=timeVaryingEffects,
                score=score,
                censModel= cens.model,
                factorLevels=c(factorLevelsX,factorLevelsZ),
                refLevels=c(refLevelsX,refLevelsZ),
                "na.action"=attr(EHF,"na.action"))
    if (confint && sim==1)
        out <- c(out,list(resampleResults=list(conf.band=unifCI,
                              B.iid=B.iid,
                              gamma.iid=gamiid,
                              test.procBeqC=Ut,
                              sim.test.procBeqC=UIt)))
    if (sim==1)
        out <- c(out,list(timeVarSigTest=timeVarSignifTest,
                          timeVarKolmSmirTest=timeVarKolmSmirTest,
                          timeVarKramvMisTest=timeVarKramvMisTest))
    if (is.null(out$call$cause))
        out$call$cause <- cause
    class(out) <- "riskRegression"
    return(out)
    # }}}
}




