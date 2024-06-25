getCensoringWeights <- function(formula,
                                data,
                                times,
                                cens.model,
                                influence.curve=FALSE,
                                censoring.save.memory,
                                verbose){
    data <- copy(data)
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
               sFormula <- update(formula,"Surv(riskRegression_time,riskRegression_status)~1")
               fit <- prodlim::prodlim(sFormula,data=data,reverse=TRUE,bandwidth="smooth")
               IPCW.times <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
               IPCW.subject.times <- prodlim::predictSurvIndividual(fit,lag=1)
               out <- list(IPCW.times=IPCW.times,
                           IPCW.subject.times=IPCW.subject.times,
                           method=cens.model,IC.data=NULL)
               if (influence.curve==TRUE){
                   # Do not yet calculate IF for Nelson-AAlen
                   # this happens later
                  out <- c(out,list(IC=NULL))
                   # out <- c(out,
                   #          list(IC=IC_Nelson_Aalen_cens_time(time=data$riskRegression_time,status=data$riskRegression_status)))
                   ## list(IC=data[,getInfluenceCurve.NelsonAalen(time=riskRegression_time,riskRegression_status=riskRegression_status)]))
               }
               out
           },"cox"={
               event = event2 = NULL
               sFormula <- update(formula,"Surv(riskRegression_time,riskRegression_status == 0)~.")
               tFormula <- update(formula,"Surv(riskRegression_time,riskRegression_status != 0)~.")
               wdata <- copy(data)
               if (length(unique(wdata[["event"]])) > 2){
                   wdata[,event2 := as.numeric(event)]
                   wdata[riskRegression_status == 0,event2 := 0]
                   wdata[,event2 := NULL]
               }
               # else {
               #     fitCSC <- NULL
               # }
               # wdata[,riskRegression_status:=1-riskRegression_status]
               Y <- data[["riskRegression_time"]]
               riskRegression_status <- data[["riskRegression_status"]]
               ## fit Cox model for censoring times
               args <- list(x=TRUE,y=TRUE,eps=0.000001,linear.predictors=TRUE)
               args$surv <- TRUE
               fit <- do.call(rms::cph,c(list(sFormula,data=wdata),args))
               if (length(times)==1){
                   IPCW.times <- matrix(rms::survest(fit,newdata=wdata,times=times,se.fit=FALSE)$surv,ncol=1)
               } else{
                   IPCW.times <- rms::survest(fit,newdata=wdata,times=times,se.fit=FALSE)$surv
               }
               ## only one per subject, so must be a flat vector
               ## FIXME: really need subject.times only where events occur before times ???
               # times.data.minus <- c(0,Y[-length(Y)])
               incrementsOfTime <- diff(sort(Y))
               minimumIncrement <- min(incrementsOfTime[incrementsOfTime > 0]) #need > 0 in case of ties
               IPCW.subject.times <- as.numeric(rms::survest(fit,newdata=wdata,times=Y-minimumIncrement/2,what='parallel',se.fit=FALSE)) #Y-minimumIncrement/2
               # subject.times <- Y[(((Y<=max(times))*riskRegression_status)==1)]
               out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model)
               if (influence.curve==TRUE){
                   TiMinus = Y-minimumIncrement/2
                   if (!censoring.save.memory){
                     ## IC is an array with dimension (nlearn, times, newdata)
                     ##                           IC_G(t,z;x_k)
                     IC.subject=predictCox(fit, iid = TRUE,
                                           newdata = wdata,
                                           times = TiMinus,
                                           diag=TRUE,
                                           type = "cumhazard")$cumhazard.iid*NROW(data) ## we want the influence function   survival.iid*n
                     IC.times=predictCox(fit, iid = TRUE,
                                         newdata = wdata,
                                         times = times,
                                         type = "cumhazard")$cumhazard.iid*NROW(data)
                     IC.weights <- list()
                     for (t.ind in 1:length(times)){
                       N <- length(Y)
                       ic.weights <- matrix(0,N,N)  ## (i,j)'th entry is f_j(tilde{T}_i-;X_i)/G(tilde{T}_i|X_i) when delta_i != 0 and time_i <= tau; otherwise it is f_j(tau-;X_i)/G(tau | X_i)
                       # k=0 ## counts subject-times with event before t
                       for (i in 1:N){
                         if (Y[i]<=times[t.ind] && riskRegression_status[i] != 0){ ## min(T,C)<=t, note that (residuals==0) => (status==0)
                           # k=k+1
                           ic.weights[i,] <- IC.subject[,1,i]
                           
                         }else if (Y[i] > times[t.ind]){## min(T,C)>t
                             ic.weights[i,] <- IC.times[,t.ind,i]
                         }
                       }
                       IC.weights[[t.ind]] <- ic.weights
                     }
                       IC <- list(censoring.save.memory = censoring.save.memory, IC.weights = IC.weights)
                   }
                   else {
                       IC <- list(censoring.save.memory = censoring.save.memory,
                                  fit=fit,
                                  wdata=wdata,
                                  TiMinus = TiMinus)
                   }
                   out <- c(out,list(IC=IC))
               }
           }, "discrete"={
               stop("Use stratification with Cox instead. ")
           },
           {
               if (verbose>0)warning("Using other models (than Cox and Kaplan-Meier) for getting the censoring weights is under construction.  ")
               vv <- all.vars(formula(delete.response(terms(formula))))
               new.formula<-as.formula(paste0("Surv(riskRegression_time,riskRegression_status)",paste0("~",paste0(paste(vv,collapse = "+")))))
               wdata <- copy(data)
               wdata[,riskRegression_status:=1-riskRegression_status]
               input <- list(formula=new.formula,data=wdata)
               if (verbose>1)message("Fitting censoring model to data ...", appendLF = FALSE)
               fit <- do.call(cens.model,input)
               if (verbose>1)message("done!")
               if (influence.curve){
                   if (is.null(data[["event"]])){
                       new.formula<-as.formula(paste0("Surv(riskRegression_time,riskRegression_status==1)",paste0("~",paste0(paste(vv,collapse = "+")))))
                   }
                   else {
                       new.formula<-as.formula(paste0("Surv(riskRegression_time,riskRegression_event==1)",paste0("~",paste0(paste(vv,collapse = "+")))))
                   }
                   input <- list(formula=new.formula,data=wdata)
                   if (verbose>0)message("Fitting time model to data ...", appendLF = FALSE)
                   fit.time <- do.call(cens.model,input)
                   if (verbose>0)message("done!")
               }
               else {
                   fit.time <- NULL
               }
               
               # IC.data <- list(Stimes = diag(1-predictRisk(fit.time,wdata,wdata$time,1)), Gtimes=diag(1-predictRisk(fit,wdata,wdata$time,1)))
               IC.data <- list(fit.time=fit.time,fit.cens=fit,wdata=wdata)

               # fit<-Hal9001(new.formula,wdata)
               # times.data.minus <- c(0,wdata$time[-length(wdata$time)]) #have to compute the weights for T_i minus not just Ti
               IPCW.subject.times <- diag(1-predictRisk(fit,wdata,wdata$time,1)) #computational problem with predictRisk
               IPCW.times <- 1-predictRisk(fit,wdata,times,1)
               out <- list(
                   IPCW.times=IPCW.times,
                   IPCW.subject.times=IPCW.subject.times,
                   method=cens.model,
                   IC.data=IC.data)
           })
    out$dim <- ifelse(cens.model=="KaplanMeier" || cens.model == "marginal",0,1)
    out
}
