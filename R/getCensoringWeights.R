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
                   out <- c(out,
                            list(IC=IC_Nelson_Aalen_cens_time(time=data$time,status=data$status)))
                   ## list(IC=data[,getInfluenceCurve.NelsonAalen(time=time,status=status)]))
               }
               out
           },"cox"={
               sFormula <- update(formula,"Surv(time,status)~.")
               tFormula <- update(formula,"Surv(time,event==1)~.")
               wdata <- copy(data)
               wdata[,status:=1-status]
               Y <- data[["time"]]
               status <- data[["status"]]
               ## fit Cox model for censoring times
               args <- list(x=TRUE,y=TRUE,eps=0.000001)
               args$surv <- TRUE
               fit <- do.call(rms::cph,c(list(sFormula,data=wdata),args))
               fit.time <- do.call(rms::cph,c(list(tFormula,data=wdata),args))
               IC.data <- list(Stimes = as.numeric(rms::survest(fit.time,times=Y,what='parallel')), Gtimes=as.numeric(rms::survest(fit,times=Y,what='parallel')))

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
                   out <- c(out,list(IC=IC))
               }
           },
           {
               vv <- all.vars(formula(delete.response(terms(formula))))
               new.formula<-as.formula(paste0("Surv(time,status)",paste0("~",paste0(paste(vv,collapse = "+")))))
               wdata <- copy(data)
               wdata[,status:=1-status]
               input <- list(formula=new.formula,data=wdata)
               fit <- do.call(cens.model,input)
               new.formula<-as.formula(paste0("Surv(time,event==1)",paste0("~",paste0(paste(vv,collapse = "+")))))
               input <- list(formula=new.formula,data=wdata)
               fit.time <- do.call(cens.model,input)
               IC.data <- list(Stimes = diag(1-predictRisk(fit.time,wdata,wdata$time,1)), Gtimes=diag(1-predictRisk(fit,wdata,wdata$time,1)))
               # IC.data <- list(fit.time=fit.time,fit.cens=fit,wdata=wdata)

               # fit<-Hal9001(new.formula,wdata)
               times.data.minus <- c(0,wdata$time[-length(wdata$time)]) #have to compute the weights for T_i minus not just Ti
               IPCW.subject.times <- diag(1-predictRisk(fit,wdata,times.data.minus,1)) #computational problem with predictRisk
               IPCW.times <- 1-predictRisk(fit,wdata,times,1)
               out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model,IC.data=IC.data)
           })
    out$dim <- ifelse(cens.model=="cox",1,0)
    out
}
