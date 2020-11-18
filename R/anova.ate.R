### anova.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug 19 2020 (09:18) 
## Version: 
## Last-Updated: okt  6 2020 (16:05) 
##           By: Brice Ozenne
##     Update #: 71
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * anova.ate (documentation)
#' @title Risk Comparison Over Time
#' @description Comparison of risk differences or risk ratios over all timepoints.
#'
#' @param object A \code{ate} object, i.e. output of the \code{ate} function.
#' @param allContrast [matrix] contrast for which the risks should be compared.
#' Matrix with two rows, the first being the sequence of reference treatments and the second the sequence of alternative treatments. 
#' @param type [character vector] the functionnal used to compare the risks: \code{"diffRisk"} or \code{"ratioRisk"}.
#' @param estimator [character] The type of estimator relative to which the comparison should be performed. 
#' @param test [character] The type of statistic used to compare the risks over times:
#' \code{"KM"} (extremum risk), \code{"CvM"} (sum of squares of the risk), or \code{"sum"} (sum of the risks).
#' @param transform [character] Should a transformation be used, e.g. the test is performed after log-transformation of the estimate, standard error, and influence function.
#' @param alternative [character] a character string specifying the alternative hypothesis, must be one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
#' @param n.sim [integer, >0] the number of simulations used to compute the p-values.
#' @param print [logical] should the results be displayed?
#' @param ... Not used.
#'
#' @details Experimental!!!

## * confint.ate (examples)
##' @examples
##' library(survival)
##' library(data.table)
##' 
##' \dontrun{
##' ## simulate data
##' set.seed(12)
##' n <- 200
##' dtS <- sampleData(n,outcome="survival")
##' dtS$X12 <- LETTERS[as.numeric(as.factor(paste0(dtS$X1,dtS$X2)))]
##' dtS <- dtS[dtS$X12!="D"]
##'
##' ## model fit
##' fit <- cph(formula = Surv(time,event)~ X1+X6,data=dtS,y=TRUE,x=TRUE)
##' seqTime <- 1:10
##' ateFit <- ate(fit, data = dtS, treatment = "X1", contrasts = NULL,
##'               times = seqTime, B = 0, iid = TRUE, se = TRUE, verbose = TRUE, band = TRUE)
##'
##' ## display
##' autoplot(ateFit)
##' 
##' ## inference (two sided)
##' statistic <- ateFit$diffRisk$estimate/ateFit$diffRisk$se
##' confint(ateFit, p.value = TRUE, method.band = "bonferroni")$diffRisk
##' confint(ateFit, p.value = TRUE, method.band = "maxT-simulation")$diffRisk
##'
##' anova(ateFit, test = "KS")
##' anova(ateFit, test = "CvM")
##' anova(ateFit, test = "sum")
##'
##' ## manual calculation (one sided)
##' n.sim <- 1e4
##' statistic <- ateFit$diffRisk[, estimate/se]
##' iid.norm <- scale(ateFit$iid$GFORMULA[["1"]]-ateFit$iid$GFORMULA[["0"]],
##'                   scale = ateFit$diffRisk$se)
##'
##' ls.out <- lapply(1:n.sim, function(iSim){
##' iG <- rnorm(NROW(iid.norm))
##' iCurve <- t(iid.norm) %*% iG
##' data.table(max = max(iCurve), L2 = sum(iCurve^2), sum = sum(iCurve),
##' maxC = max(iCurve) - max(statistic),
##' L2C = sum(iCurve^2) - sum(statistic^2),
##' sumC = sum(iCurve) - sum(statistic),
##' sim = iSim)
##' })
##'
##' dt.out <- do.call(rbind,ls.out)
##' dt.out[,.(max = mean(.SD$maxC>=0),
##'           L2 = mean(.SD$L2C>=0),
##'           sum = mean(.SD$sumC>=0))]
##' 
##' ## permutation
##' n.sim <- 250
##' stats.perm <- vector(mode = "list", length = n.sim)
##' pb <- txtProgressBar(max = n.sim, style=3)
##' treatVar <- ateFit$variables["treatment"]
##' 
##' for(iSim in 1:n.sim){ ## iSim <- 1
##' iData <- copy(dtS)
##' iIndex <- sample.int(NROW(iData), replace = FALSE)
##' iData[, c(treatVar) := .SD[[treatVar]][iIndex]]
##'
##' iFit <- update(fit, data = iData)
##' iAteSim <- ate(iFit, data = iData, treatment = treatVar,
##'                times = seqTime, verbose = FALSE)
##' iStatistic <- iAteSim$diffRisk[,estimate/se]
##' stats.perm[[iSim]] <- cbind(iAteSim$diffRisk[,.(max = max(iStatistic),
##'                                                 L2 = sum(iStatistic^2),
##'                                                 sum = sum(iStatistic))],
##'                             sim = iSim)
##' stats.perm[[iSim]]$maxC <- stats.perm[[iSim]]$max - max(statistic)
##' stats.perm[[iSim]]$L2C <- stats.perm[[iSim]]$L2 - sum(statistic^2)
##' stats.perm[[iSim]]$sumC <- stats.perm[[iSim]]$sum - sum(statistic)
##' setTxtProgressBar(pb, iSim)
##' }
##'
##' dtstats.perm <- do.call(rbind,stats.perm)
##' dtstats.perm[,.(max = mean(.SD$maxC>=0),
##'                 L2 = mean(.SD$L2C>=0),
##'                 sum = mean(.SD$sumC>=0))]
##' }

## * anova.ate (code)
#' @rdname anova.ate
#' @method anova ate
#' @export
anova.ate <- function(object,
                      allContrast = NULL, type = "diff", estimator = object$estimator[1],
                      test = "CvM", transform = NULL,  alternative = "two.sided", n.sim = 1e4,
                      print = TRUE, ...){

    ## ** initialize and check arguments
    if(is.null(allContrast)){
        contrasts <- object$contrasts
        allContrasts <- utils::combn(contrasts, m = 2)
    }
    n.allContrasts <- NCOL(allContrasts)
    type <- match.arg(type, choices = c("diff", "ratio"))
    null <- switch(type, "diff" = 0, "ratio" = 1)
    estimator <- match.arg(estimator, choices = object$estimator)
    if(object$inference$se==FALSE){
        stop("Cannot perform the test without the standard errors \n",
             "Set argument \'se\' to TRUE when calling the ate function \n")
    }
    if(is.null(object$iid)){
        stop("Cannot perform the test without the standard errors \n",
             "Set argument \'iid\' to TRUE when calling the ate function \n")
    }
    if(is.null(transform)){
        if(!is.null(object$transform[[paste0(type,"Risk")]])){
            transform <- "none"
        }else{
            transform <- object$transform[[paste0(type,"Risk")]]
        }
    }
    n.sample <- NROW(object$iid[[1]][[1]])
    n.time <- NCOL(object$iid[[1]][[1]])
    test <- match.arg(test, c("KS","CvM","sum"))
    alternative <- match.arg(alternative, c("two.sided","greater","less"))
    
    ## ** prepare arguments for cpp routine
    beta <- matrix(NA, nrow = n.time, ncol = n.allContrasts)
    iid2cpp <- array(NA, dim = c(n.time, n.sample, n.allContrasts))
    statistic2cpp <- matrix(NA, nrow = n.time, ncol = n.allContrasts)

    typeRisk <- object[[paste0(type,"Risk")]]
    
    for(iC in 1:n.allContrasts){## iC <- 1
        ## *** gather information
        iC.A <- allContrasts[1,iC]
        iC.B <- allContrasts[2,iC]

        iRowIndex <- which((typeRisk$estimator==estimator)*(typeRisk$A==iC.A)*(typeRisk$B==iC.B)==1)
        beta[,iC] <- typeRisk[iRowIndex, .SD$estimate]
        beta.se <- typeRisk[iRowIndex, .SD$se]

        iid.A <- object$iid[[estimator]][[iC.A]]
        iid.B <- object$iid[[estimator]][[iC.B]]
        
        if(type=="diff"){
            iid.AB <- iid.B-iid.A
        }else if(type=="ratio"){
            iRowIndex <- which((object$meanRisk$estimator==estimator)*(object$meanRisk$treatment==iC.A)==1)
            risk.A <- object$meanRisk[iRowIndex, .SD$estimate]
            iRowIndex <- which((object$meanRisk$estimator==estimator)*(object$meanRisk$treatment==iC.B)==1)
            risk.B <- object$meanRisk[iRowIndex, .SD$estimate]
            
            iid.AB <- rowScale_cpp(iid.B, scale = risk.A)-rowScale_cpp(iid.A, scale = risk.B/risk.A^2)
        }
        ## *** transformation
        beta.se <- transformSE(estimate = beta[,iC], se = beta.se, type = transform)
        statistic2cpp[,iC] <- transformT(estimate = beta[,iC], se = beta.se, null = null, type = transform)
        iid2cpp[,,iC] <- t(rowScale_cpp(transformIID(estimate = beta[,iC], iid = iid.AB, type = transform), scale = beta.se))
    }

    ## ** run simulation
    resCpp <- sampleMaxProcess_cpp(nSample = n.sample,
                                   nContrast = n.allContrasts,
                                   nSim = n.sim,
                                   value = statistic2cpp,
                                   iid =  iid2cpp,
                                   global = FALSE,
                                   alternative = switch(alternative,
                                                        "two.sided" = 3,
                                                        "greater" = 2,
                                                        "less" = 1),
                                   type = switch(test, "KS"=1, "CvM"=2, "sum"=3)
                                   )

    ## ** process results
    out <- data.frame(matrix(NA, nrow = n.allContrasts, ncol = 4,
                             dimnames = list(NULL,c("treatment.A","treatment.B","statistic","p.value"))))
    out$treatment.A <- allContrasts[1,]
    out$treatment.B <- allContrasts[2,]
    if(test == "KS"){
        if(alternative == "less"){
            out$statistic <- apply(statistic2cpp,2,min)
        }else if(alternative == "greater"){
            out$statistic <- apply(statistic2cpp,2,max)
        }else if(alternative == "two.sided"){
            out$statistic <- apply(abs(statistic2cpp),2,max)
        }
    }else if(test == "CvM"){
        out$statistic <- colSums(statistic2cpp^2)
    }else if(test == "sum"){
        out$statistic <- colSums(statistic2cpp)
    }
    out$p.value <- colMeans(resCpp>0)

    ## ** display
    if(print){
        txt.alternative <- switch(alternative,
                                  "two.sided" = "unequal mean risk between treatment groups",
                                  "greater" = "greater mean risk with treatment B compared to treatment A",
                                  "less" = "lower mean risk with treatment B compared to treatment A")
        txt.test <- switch(test,
                           "KS" = "maximum",
                           "CvM" = "sum of squares",
                           "sum" = "sum")
        txt.type <- switch(type,
                           "diff" = "differences",
                           "ratio" = "ratios")

        cat("    Comparison of the Average Treatment Effect over all times \n\n")
        cat(" - test statistic        : ",txt.test," of the risk ",txt.type," \n",sep="")
        cat(" - alternative hypothesis: ",txt.alternative," \n",sep="")
        print(out)
    }
    
    ## ** export
    return(invisible(out))
        
}

######################################################################
### anova.ate.R ends here
