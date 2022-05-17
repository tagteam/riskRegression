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
                           method=cens.model)
               if (influence.curve==TRUE){
                   out <- c(out,
                            list(IC=IC_Nelson_Aalen_cens_time(time=data$time,status=data$status)))
                   ## list(IC=data[,getInfluenceCurve.NelsonAalen(time=time,status=status)]))
               }
               out
           },"cox"={
               sFormula <- update(formula,"Surv(time,status)~.")
               wdata <- copy(data)
               wdata[,status:=1-status]
               Y <- data[["time"]]
               status <- data[["status"]]
               ## fit Cox model for censoring times
               args <- list(x=TRUE,y=TRUE,eps=0.000001)
               args$surv <- TRUE
               fit <- do.call(rms::cph,c(list(sFormula,data=wdata),args))
               ## need G(Ti-|Xi) only for i where status=1 && Ti < max(times)
               subject.times <- Y[(((Y<=max(times))*status)==1)]
               if (length(times)==1){
                   IPCW.times <- matrix(rms::survest(fit,newdata=wdata,times=times,se.fit=FALSE)$surv,ncol=1)
               } else{
                   IPCW.times <- rms::survest(fit,newdata=wdata,times=times,se.fit=FALSE)$surv
               }
               ## only one per subject, so must be a flat vector
               ## FIXME: really need subject.times only where events occur before times
               IPCW.subject.times <- as.numeric(rms::survest(fit,times=Y-min(diff(c(0,unique(Y))))/2,what='parallel'))
               out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model)
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
           "hal9001"={
               vv <- all.vars(formula(delete.response(terms(formula))))
               formula.covariates<-as.formula(paste0("~",paste0(paste(vv,collapse = "+"),"-1")))
               wdata <- copy(data)
               wdata[,status:=1-status]
               covariates.data <- wdata[,..vv]
               x.surv <- model.matrix(formula.covariates,data=wdata)
               y.surv <- Surv(wdata$time,wdata$status)
               ## fit Cox model for censoring times
               fit<-hal9001::fit_hal(X=x.surv,Y=y.surv,family="cox",yolo=FALSE)
               surv_info <- list()
               surv_info$start <- 0
               surv_info$stop <- wdata[["time"]]
               surv_info$status <- wdata[["status"]]
               fit$surv_info <- surv_info
               ## need G(Ti-|Xi) only for i where status=1 && Ti < max(times)
               IPCW.subject.times <- diag(1-predictRisk(fit,x.surv,wdata[["time"]],1))
               IPCW.times <- 1-predictRisk(fit,x.surv,times,1)
               out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model)
           },
           {
               stop("IPCW works only for nuisance model obtained with Kaplan-Meier (marginal) or Cox regression (cox).")
           })
    out$dim <- ifelse(cens.model=="cox",1,0)
    out
}
