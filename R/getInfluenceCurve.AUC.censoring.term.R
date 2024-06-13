### getInfluenceCurve.AUC.censoring.term.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2024 (11:47) 
## Version: 
## Last-Updated: Jun 13 2024 (17:28) 
##           By: Thomas Alexander Gerds
##     Update #: 9
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getInfluenceFunction.AUC.censoring.term <- function(time,
                                                    event,
                                                    t,
                                                    IFcalculationList,
                                                    MC,
                                                    cens.model,
                                                    Wt,
                                                    auc,
                                                    nth.times){
    if (cens.model[[1]] == "KaplanMeier"){
        ind.controls<-rep(NA,length(time))
        controls.index1 <- IFcalculationList[["controls1"]]
        controls.index2 <- IFcalculationList[["controls2"]]
        if ((sum(controls.index1)+sum(controls.index2)) == 0){
            return(numeric(length(time)))
        }
        ind.controls[controls.index1] <- 1
        ind.controls[controls.index2] <- 0
        start.controls1 <- sindex(ind.controls[controls.index1 | controls.index2],0)
        getInfluenceFunctionAUCKMCensoringTerm(time,
                                               event,
                                               t,
                                               IFcalculationList[["ic0Case"]],
                                               IFcalculationList[["ic0Control"]],
                                               IFcalculationList[["weights"]],
                                               IFcalculationList[["firsthit"]],
                                               IFcalculationList[["muCase"]],
                                               IFcalculationList[["muControls"]], 
                                               IFcalculationList[["nu"]],
                                               Wt[1],
                                               auc,
                                               start.controls1)
    }
    else if (cens.model[[1]] == "cox"){
        n <- length(time)
        cases.index <- IFcalculationList[["cases"]]
        controls.index <- IFcalculationList[["controls"]]
        if ((sum(controls.index)+sum(controls.index)) == 0){
            return(numeric(length(time)))
        }
        ic0Case <- IFcalculationList[["ic0Case"]]
        ic0Control <- IFcalculationList[["ic0Control"]]
        Phi <- IFcalculationList[["muCase"]] * IFcalculationList[["muControls"]] / (n*n)
        weights <- IFcalculationList[["weights"]]
        muCase <- IFcalculationList[["muCase"]]
        muControls <- IFcalculationList[["muControls"]]
        aucLPO <- auc
        w.cases <- weights[cases.index]
        w.controls <- weights[controls.index]
        if (!MC$censoring.save.memory){
            ic.weights <- MC[[2]][[nth.times]] ## load IF from Censoring weights
            icPart <- as.numeric(rowSumsCrossprod(as.matrix(1/(Phi*n^2)*ic0Case-(aucLPO/Phi)*(1/n^2)*muControls*w.cases), ic.weights[cases.index,], 0)) +
                as.numeric(rowSumsCrossprod(as.matrix(1/(Phi*n^2)*ic0Control-(aucLPO/Phi)*(1/n^2)*muCase*w.controls),ic.weights[controls.index,],0))
        }
        else {
            wdata <- MC[[3]]
            fit <- MC[[2]]
            TiMinus <- MC[[4]]
            ic0CaseOld <- rep(0,n)
            ic0CaseOld[cases.index] <- ic0Case
            ic0ControlOld <- rep(0,n)
            ic0ControlOld[controls.index] <- ic0Control
            controls.index1 <- IFcalculationList[["controls1"]]
            controls.index2 <- IFcalculationList[["controls2"]]
            Wbeforet <- (1/(Phi*n^2))*(ic0CaseOld*cases.index+ic0ControlOld*controls.index2)-
                (1/n)*(aucLPO/Phi)*((cases.index)*weights*(1/n)*muControls + (controls.index2)*weights*(1/n)*muCase)
            Waftert <- (1/(Phi*n^2))*ic0ControlOld*controls.index1- 
                (1/n)*(aucLPO/Phi)*(controls.index1)*weights*(1/n)*muCase
            ## First term gives for i'th entry: 1/n \sum_j weights[j] * \hat{f}_i(\tilde{T}_j-,X_j); 
            ## Next one does: 1/n \sum_j weights[j] * \hat{f}_i(tau,X_j) for Cox
            icPart <- predictCoxWeights(fit, diag=TRUE,newdata = wdata, times = TiMinus,weights=Wbeforet, isBeforeTau = TRUE, tau = t)+
                predictCoxWeights(fit, diag=FALSE,newdata = wdata,times = t,weights=Waftert)
        }
        icPart
    }
    else {
        warning("Censoring model not yet implemented. Reverting to conservative = TRUE for AUC. ")
        0
    }
}
######################################################################
### getInfluenceCurve.AUC.censoring.term.R ends here
