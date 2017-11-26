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
               fit <- prodlim::prodlim(sFormula,data=data,reverse=TRUE)
               IPCW.times <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
               IPCW.subject.times <- prodlim::predictSurvIndividual(fit,lag=1)
               out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model)
               if (influence.curve==TRUE){
                   out <- c(out,list(IC=data[,getInfluenceCurve.KM(time=time,status=status)]))
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
               subject.position <- (((Y<=max(times))*status)==1)
               subject.times <- Y[subject.position]
               IPCW.times <- rms::survest(fit,newdata=wdata,times=times,se.fit=FALSE)$surv
               IPCW.subject.times <- rms::survest(fit,times=subject.times-min(diff(c(0,unique(subject.times))))/2,what='parallel')
               out <- list(IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,method=cens.model)
               ## array: nlearn x times x newdata
               if (influence.curve==TRUE){
                   IC <- predictCox(fit, iid = TRUE,
                                    newdata = wdata,
                                    times = c(subject.times,times),
                                    log.transform = FALSE,
                                    type = "survival")$survival.iid
                   out <- c(out,list(IC=IC))
               }
           },{
               stop("IPCW works only for nuisance model obtained with Kaplan-Meier (marginal) or Cox regression (cox).")
           })
    out$dim <- ifelse(cens.model=="cox",0,1)
    out
}
