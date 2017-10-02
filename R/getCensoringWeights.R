getCensoringWeights <- function(formula,
                                data,
                                response,
                                times,
                                censModel,
                                responseType){
    getIPCW.marginal.competing.risk.AUC <- NULL
    getIPCW.marginal.competing.risk.cindex <- NULL
    if((censModel != "KaplanMeier")){
        if (length(attr(terms(formula),"factors"))==0){
            censModel <- "marginal"
        }
    }
    else{
        censModel <- "marginal"
    }
    if (responseType=="competing.risks"){
        iFormula <- as.formula(paste("Surv(itime,istatus)","~",as.character(formula)[[3]]))
        iData <- data
        iData$itime <- response[["time"]]
        iData$istatus <- response[["status"]]
        weights <- ipcw(formula=iFormula,
                        data=iData,
                        method=censModel,
                        times=times,
                        subject.times=iData$itime,
                        lag=1,keep=NULL)
    }
    else{
        weights <- ipcw(formula=formula,
                        data=data,
                        method=censModel,
                        times=times,
                        subject.times=data[["time"]],
                        lag=1,keep=NULL)
    }
    weights$dim <- if (censModel %in% c("marginal","none")) 0 else 1
    weights
}
