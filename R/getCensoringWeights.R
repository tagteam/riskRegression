getCensoringWeights <- function(formula,
                                data,
                                response,
                                times,
                                cens.model,
                                response.type,
                                influence.curve=FALSE){
    if((cens.model != "KaplanMeier")){
        if (length(attr(terms(formula),"factors"))==0){
            cens.model <- "marginal"
        }
    }
    else{
        cens.model <- "marginal"
    }
    switch(cens.model,
           "marginal"={
               sFormula <- update(formula,"Surv(time,status)~1")
               fit <- prodlim::prodlim(sFormula,data=data,reverse=TRUE,bandwidth="smooth")
               IPCW.times <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
               IPCW.subject.times <- prodlim::predictSurvIndividual(fit,lag=1)
               out <- list(IPCW.times=IPCW.times,
                           IPCW.subject.times=IPCW.subject.times,
                           method=cens.model,IC.data=NULL)
               if (influence.curve==TRUE){
                  #Should not have to calculate IF for Nelson-AAlen
                  out <- c(out,
                           list(IC=NULL))
                  
                   # out <- c(out,
                   #          list(IC=IC_Nelson_Aalen_cens_time(time=data$time,status=data$status)))
                   ## list(IC=data[,getInfluenceCurve.NelsonAalen(time=time,status=status)]))
               }
               out
           },"cox"={
               event = event2 = NULL
               sFormula <- update(formula,"Surv(time,status == 0)~.")
               tFormula <- update(formula,"Surv(time,status != 0)~.")
               wdata <- copy(data)
               if (length(unique(wdata[["event"]])) > 2){
                   wdata[,event2 := as.numeric(event)]
                   wdata[status == 0,event2 := 0]
                   FormulaCSC <- update(formula,"Hist(time,event2)~.")
                   suppressWarnings(fitCSC <- CSC(formula = FormulaCSC, data = wdata))
                   wdata[,event2 := NULL]
               }
               else {
                   fitCSC <- NULL
               }
               # wdata[,status:=1-status]
               Y <- data[["time"]]
               status <- data[["status"]]
               ## fit Cox model for censoring times
               args <- list(x=TRUE,y=TRUE,eps=0.000001)
               args$surv <- TRUE
               fit <- do.call(rms::cph,c(list(sFormula,data=wdata),args))
               # fit.time <- do.call(rms::cph,c(list(tFormula,data=wdata),args))
               IC.data <- NULL #list(wdata=wdata,fit.time=fit.time,fit.cens=fit,fitCSC=fitCSC)

               ## need G(Ti-|Xi) only for i where status=1 && Ti < max(times)
               if (length(times)==1){
                   IPCW.times <- matrix(rms::survest(fit,newdata=wdata,times=times,se.fit=FALSE)$surv,ncol=1)
               } else{
                   IPCW.times <- rms::survest(fit,newdata=wdata,times=times,se.fit=FALSE)$surv
               }
               ## only one per subject, so must be a flat vector
               ## FIXME: really need subject.times only where events occur before times ???
               times.data.minus <- c(0,Y[-length(Y)])
               IPCW.subject.times <- as.numeric(rms::survest(fit,times=times.data.minus,what='parallel'))
               subject.times <- Y[(((Y<=max(times))*status)==1)]
               ## TRY TO CORRECT LARGE WEIGHTS
               warn.weights <- FALSE
               f <- function(x) { ### HAVE TO NOT DIVIDE BY ZERO
                 if (x==0){
                   1
                 }
                 else {
                   1/x
                 }
               }
               ### FOR THE SUBJECT WEIGHTS
               for (t.ind in 1:length(times)){
                 number.of.censorings <- 0
                 upperTie <- 0
                 for (i in 1:length(Y)){
                   newcensorings <- 0
                   if (i > upperTie){
                     tieIter <- i
                     while(Y[tieIter] == Y[i] && tieIter < length(Y)){
                       if (status[tieIter] == 0 && Y[i] < times[t.ind]){
                         newcensorings <- newcensorings + 1
                       }
                       tieIter <- tieIter + 1
                     }
                     upperTie <- tieIter - 1
                   }
                   if (Y[i] <= times[t.ind] && status[i] != 0){
                     if (IPCW.subject.times[i] < f(number.of.censorings)){ ### THIS IS TYPICALLY TRUE FOR THE FIRST OBSERVATION!
                       warn.weights <- TRUE
                       IPCW.subject.times[i] <- f(number.of.censorings)
                     }
                   }
                   number.of.censorings <- number.of.censorings + newcensorings
                 }
                 ## FOR THE WEIGHTS GTAU
                 for (i in 1:length(Y)){
                   if (length(times)==1){
                     if (IPCW.times[i] < f(number.of.censorings)){ ### THIS IS TYPICALLY TRUE FOR THE FIRST OBSERVATION!
                       warn.weights <- TRUE
                       IPCW.times[i] <- f(number.of.censorings)
                     }
                   } else{
                     if (IPCW.times[i,t.ind] < f(number.of.censorings)){ ### THIS IS TYPICALLY TRUE FOR THE FIRST OBSERVATION!
                       warn.weights <- TRUE
                       IPCW.times[i,t.ind] <- f(number.of.censorings)
                     }
                   }
                 }
               }
               if (warn.weights){
                 warning("Some inverse censoring weights were found to be large.")
               }
               
               out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model,IC.data=IC.data)
               if (influence.curve==TRUE){
                   ## IC is an array with dimension (nlearn, times, newdata)
                   ##                           IC_G(t,z;x_k)
                   IC <- list(IC.subject=predictCox(fit, iid = TRUE,
                                                    newdata = wdata,
                                                    times = subject.times,
                                                    type = "survival")$survival.iid*NROW(data), ## we want the influence function   survival.iid*n
                              IC.times=predictCox(fit, iid = TRUE,
                                                  newdata = wdata,
                                                  times = times,
                                                  type = "survival")$survival.iid*NROW(data))  ## we want the influence function = survival.iid*n
                   ## IC <- predictCox(fit, iid = TRUE,newdata = wdata,times = c(subject.times,times),type = "survival")$survival.iid
                   IC.weights <- list()
                   for (t.ind in 1:length(times)){
                     N <- length(Y)
                     ic.weights <- matrix(0,N,N)
                     k=0 ## counts subject-times with event before t
                     for (i in 1:N){
                       if (Y[i]<=times[t.ind] && status[i] != 0){ ## min(T,C)<=t, note that (residuals==0) => (status==0)
                         k=k+1
                         ic.weights[i,] <- IC$IC.subject[i,k,]/IPCW.subject.times[i]
                       }else if (Y[i] > times[t.ind]){## min(T,C)>t
                         if (length(times) == 1){
                           ic.weights[i,] <- IC$IC.times[i,t.ind,]/IPCW.times[i]
                         }
                         else {
                           ic.weights[i,] <- IC$IC.times[i,t.ind,]/IPCW.times[i,t.ind]
                         }
                       }
                     }
                     IC.weights[[t.ind]] <- ic.weights
                   }
                   out <- c(out,list(IC=IC.weights))
               }
           }, "discrete"={
             stop("Do not use discrete option.")
             # vv <- all.vars(formula(delete.response(terms(formula))))
             # mod.frame<-model.frame(paste0("~",paste0(paste(vv,collapse = "+"))),data=data)
             # u.levels <- unique(mod.frame)
             # n <- length(data$time)
             # prob.risk <- rep(NA,n)
             # have.same.covariate <- list()
             # for (i in 1:nrow(u.levels)){
             #   wheres <- rep(NA,n)
             #   for (j in 1:n){
             #     if (all(mod.frame[j,]==u.levels[i,])){
             #       wheres[j] <- TRUE
             #       
             #     }
             #     else {
             #       wheres[j] <- FALSE
             #     }
             #   }
             #   prob.risk[wheres] <- sum(wheres)/n
             #   for (j in 1:n){
             #     if (wheres[j]){
             #       have.same.covariate[[j]] <- wheres
             #     }
             #   }
             # }
             # new.formula<-as.formula(paste0("Surv(time,status)",paste0("~",paste0(paste(vv,collapse = "+")))))
             # wdata <- copy(data)
             # fit.cens <- prodlim::prodlim(new.formula, data=wdata,reverse=TRUE)
             # fit.time <- prodlim::prodlim(new.formula, data=wdata)
             # IC.data <- list(Stimes=prodlim::predictSurvIndividual(fit.time,lag=0),Gtimes=prodlim::predictSurvIndividual(fit.cens,lag=0), Stau = c(1-predictRisk(fit.time,wdata,times,1)), prob.risk=prob.risk,have.same.covariate=have.same.covariate)
             # IPCW.subject.times <- prodlim::predictSurvIndividual(fit.cens,lag=1) #computational problem with predictRisk
             # IPCW.times <- predict(fit.cens,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
             # out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model,IC.data=IC.data)
           },
           {
               warning("Using other models (than Cox) for getting the censoring weights is under construction.  ")
               vv <- all.vars(formula(delete.response(terms(formula))))
               new.formula<-as.formula(paste0("Surv(time,status)",paste0("~",paste0(paste(vv,collapse = "+")))))
               wdata <- copy(data)
               wdata[,status:=1-status]
               input <- list(formula=new.formula,data=wdata)
               message("Fitting censoring model to data ...", appendLF = FALSE)
               fit <- do.call(cens.model,input)
               message("done!")
               if (influence.curve){
                 if (is.null(data[["event"]])){
                   new.formula<-as.formula(paste0("Surv(time,status==1)",paste0("~",paste0(paste(vv,collapse = "+")))))
                 }
                 else {
                   new.formula<-as.formula(paste0("Surv(time,event==1)",paste0("~",paste0(paste(vv,collapse = "+")))))
                 }
                 input <- list(formula=new.formula,data=wdata)
                 message("Fitting time model to data ...", appendLF = FALSE)
                 fit.time <- do.call(cens.model,input)
                 message("done!")
               }
               else {
                 fit.time <- NULL
               }
               
               # IC.data <- list(Stimes = diag(1-predictRisk(fit.time,wdata,wdata$time,1)), Gtimes=diag(1-predictRisk(fit,wdata,wdata$time,1)))
               IC.data <- list(fit.time=fit.time,fit.cens=fit,wdata=wdata)

               # fit<-Hal9001(new.formula,wdata)
               times.data.minus <- c(0,wdata$time[-length(wdata$time)]) #have to compute the weights for T_i minus not just Ti
               IPCW.subject.times <- diag(1-predictRisk(fit,wdata,times.data.minus,1)) #computational problem with predictRisk
               IPCW.times <- 1-predictRisk(fit,wdata,times,1)
               out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model,IC.data=IC.data)
           })
    out$dim <- ifelse(cens.model=="cox",1,0)
    out
}
