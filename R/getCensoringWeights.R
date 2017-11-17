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
    sFormula <- update(formula,"Surv(time,status)~.")
    if (influence.curve==TRUE && cens.model=="cox")
        keep <- "fit"
    else keep <- NULL
    weights <- ipcw(formula=sFormula,
                    data=data,
                    method=cens.model,
                    times=times,
                    subject.times=data[["time"]],
                    lag=1,
                    keep=keep)
    ## browser(skipCalls=1L)
    if (influence.curve==TRUE){
        switch(cens.model,
               "marginal"={
                   IC <- data[,getInfluenceCurve.KM(time=time,
                                                    status=status)]
               },
               "cox"={
                   ## array: nlearn x times x newdata
                   IC <- predictCox(weights$fit,
                                    iid=TRUE,
                                    newdata=data,
                                    times=data[["time"]],
                                    log.transform=FALSE,
                                    type="survival")$survival.iid
               })
        weights <- c(weights,list(influence.curve=IC))
    }
    weights$dim <- if (cens.model %in% c("marginal","none")) 0 else 1
    weights
}
