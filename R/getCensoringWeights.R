getCensoringWeights <- function(formula,
                                data,
                                response,
                                times,
                                cens.model,
                                responseType){
    getIPCW.marginal.competing.risk.AUC <- NULL
    getIPCW.marginal.competing.risk.cindex <- NULL
    if((cens.model != "KaplanMeier")){
        if (length(attr(terms(formula),"factors"))==0){
            cens.model <- "marginal"
        }
    }
    else{
        cens.model <- "marginal"
    }
    if (responseType=="competing.risks"){
        iFormula <- as.formula(paste("Surv(itime,istatus)","~",as.character(formula)[[3]]))
        iData <- data
        iData$itime <- response[["time"]]
        iData$istatus <- response[["status"]]
        weights <- ipcw(formula=iFormula,
                        data=iData,
                        method=cens.model,
                        times=times,
                        subject.times=iData$itime,
                        lag=1,keep=NULL)
    }
    else{
        weights <- ipcw(formula=formula,
                        data=data,
                        method=cens.model,
                        times=times,
                        subject.times=data[["time"]],
                        lag=1,keep=NULL)
    }
    weights$dim <- if (cens.model %in% c("marginal","none")) 0 else 1
    weights
}
