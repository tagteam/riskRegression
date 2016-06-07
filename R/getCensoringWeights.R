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
        iData$itime <- response[,"time",with=FALSE]
        iData$istatus <- response[,"status",with=FALSE]
        weights <- ipcw(formula=iFormula,data=iData,method=censModel,times=times,subjectTimes=iData$itime,subjectTimesLag=1)
    }
    else{
        weights <- ipcw(formula=formula,
                        data=data,
                        method=censModel,
                        times=times,
                        subjectTimes=data[["time"]],
                        subjectTimesLag=1)
    }
    weights$dim <- if (censModel %in% c("marginal","none")) 0 else 1
    weights
}
