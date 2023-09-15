### AUC.survival.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 11 2022 (17:06)
## Version:
## Last-Updated: Jun 30 2023 (16:09)
##           By: Thomas Alexander Gerds
##     Update #: 13
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:

AUC.survival <- function(DT,
                         breaks = NULL,
                         MC,
                         se.fit,
                         conservative,
                         cens.model,
                         keep.vcov = FALSE,
                         keep.iid = FALSE,
                         multi.split.test,
                         alpha,
                         N,
                         NT,
                         NF,
                         dolist,
                         ROC,
                         IC.data,
                         cutpoints,
                         ...) {
  ID <- model <- times <- risk <- Cases <- time <- status <- Controls <- TPR <- FPR <- WTi <- Wt <- ipcwControls <- ipcwCases <- IF.AUC <- lower <- se <- upper <- AUC <- nth.times <- NULL
  cause <- 1
  aucDT <- DT[model > 0]
  ## remove null model comparisons
  dolist <- dolist[sapply(dolist, function(do) {
    match("0", do, nomatch = 0L)
  }) == 0]
  ## assign Weights before ordering
  aucDT[, ipcwControls := 1 / (Wt * N)]
  aucDT[, ipcwCases := 1 / (WTi * N)]
  ## order data
  data.table::setorder(aucDT, model, times, -risk)
  ## identify cases and controls
  aucDT[, Cases := (time <= times & status == cause)]
  aucDT[, Controls := (time > times)]
  ## prepare Weights
  aucDT[Cases == 0, ipcwCases := 0]
  aucDT[Controls == 0, ipcwControls := 0]
  ## compute denominator
  aucDT[, TPR := cumsum(ipcwCases) / sum(ipcwCases), by = list(model, times)] # technically sum_i I(M_i >= M_j) not M_i > M_j
  aucDT[, FPR := (cumsum(ipcwControls)) / (sum(ipcwControls)), by = list(model, times)]
  nodups <- aucDT[, c(!duplicated(risk)[-1], TRUE), by = list(model, times)]$V1
  AireTrap <- function(FP, TP) {
    N <- length(FP)
    sum((FP - c(0, FP[-N])) * ((c(0, TP[-N]) + TP) / 2))
  }
  score <- aucDT[nodups, list(AUC = AireTrap(FPR, TPR)), by = list(model, times)]
  if (!is.null(cutpoints)) {
    breaks <- sort(cutpoints, decreasing = TRUE)
    aucDT[, nth.times := as.numeric(factor(times))]
    cutpoint.helper.fun <- function(FPR, TPR, risk, ipcwCases, ipcwControls, N, time, times, event, cens.model, nth.times, conservative, IC.G, cutpoints, se.fit) {
      den_TPR <- sum(ipcwCases) ## estimate the cumulative incidence via IPCW
      den_FPR <- sum(ipcwControls)
      indeces <- sindex(risk, cutpoints, comp = "greater", TRUE)
      res <- list()
      ordered <- order(time) ## can probably move this outside to improve computation time, for now keep it
      SE.TPR <- SE.FPR <- SE.PPV <- SE.NPV <- NA
      for (i in seq_len(length(cutpoints))) {
        den_PPV <- sum(ipcwCases[risk > cutpoints[i]] + ipcwControls[risk > cutpoints[i]])
        den_NPV <- 1 - den_PPV
        if (indeces[i] != 0) {
          TPRi <- TPR[indeces[i]]
          FPRi <- FPR[indeces[i]]
          if (se.fit) {
            IC0.TPR <- ipcwCases * N * ((risk > cutpoints[i]) - TPRi) / den_TPR
            SE.TPR <- sd(getInfluenceCurve.Brier(times, time[ordered], IC0.TPR[ordered], IC0.TPR[ordered], IC.G, cens.model, nth.times, conservative, event[ordered])) / sqrt(N)
            IC0.FPR <- (ipcwControls) * N * ((risk > cutpoints[i]) - FPRi) / (1 - den_TPR)
            SE.FPR <- sd(getInfluenceCurve.Brier(times, time[ordered], IC0.FPR[ordered], IC0.FPR[ordered], IC.G, cens.model, nth.times, conservative, event[ordered])) / sqrt(N)
          }
        } else {
          TPRi <- FPRi <- 0
        }
        if (den_PPV > 1e-10) {
          PPV <- (TPRi * den_TPR) / den_PPV
          if (se.fit) {
            IC0.PPV <- (risk > cutpoints[i]) / den_PPV * (((ipcwCases) * N) * (1 * (event == 1) - 1 * (event != 0) * PPV) - ipcwControls * N * PPV) # OBS, check other causes, paul's implementation
            SE.PPV <- sd(getInfluenceCurve.Brier(times, time[ordered], IC0.PPV[ordered], IC0.PPV[ordered], IC.G, cens.model, nth.times, conservative, event[ordered])) / sqrt(N)
          }

          ## Alternative implementation
          # IC0.PPV <- (risk > cutpoints[i])/den_PPV*(ipcwCases*N - PPV)
          # weights.PPV <- ((risk > cutpoints[i])/(den_PPV))*((1-PPV)*ipcwCases*N - PPV*(N*ipcwControls1+N*ipcwControls2)) #include censoring from denominator
          # weights.PPV <- ((risk > cutpoints[i])/(den_PPV))*(ipcwCases*N) #exclude censoring from denominator
          # SE.PPV <- sd(getInfluenceCurve.Brier(times,time[ordered],IC0.PPV[ordered],weights.PPV[ordered],IC.G,cens.model,nth.times,conservative,event[ordered]))/sqrt(N)
        } else {
          PPV <- NA
        }
        if (den_NPV > 1e-10) {
          NPV <- ((1 - FPRi) * den_FPR) / den_NPV
          if (se.fit) {
            IC0.NPV <- (risk <= cutpoints[i]) / den_NPV * (((ipcwCases) * N) * (1 * (event != 1 & event != 0) - 1 * (event != 0) * NPV) + ipcwControls * N * (1 - NPV)) # OBS, check other causes, paul's implementation
            SE.NPV <- sd(getInfluenceCurve.Brier(times, time[ordered], IC0.NPV[ordered], IC0.NPV[ordered], IC.G, cens.model, nth.times, conservative, event[ordered])) / sqrt(N)
          }

          ## Alternative implementation
          # IC0.NPV <- (risk <= cutpoints[i])/den_NPV*((ipcwControls1+ipcwControls2)*N - NPV)
          # weights.NPV <- ((risk <= cutpoints[i])/(den_NPV))*((ipcwControls1+ipcwControls2)*N)
          # SE.NPV <- sd(getInfluenceCurve.Brier(times,time[ordered],IC0.NPV[ordered],weights.NPV[ordered],IC.G,cens.model,nth.times,conservative,event[ordered]))/sqrt(N)
        } else {
          NPV <- NA
        }
        res[[i]] <- data.table(risk = cutpoints[i], TPR = TPRi, SE.TPR = SE.TPR, FPR = FPRi, SE.FPR = SE.FPR, PPV = PPV, SE.PPV = SE.PPV, NPV = NPV, SE.NPV = SE.NPV)
      }
      do.call("rbind", res)
    }
    output <- list(res.cut = aucDT[, cutpoint.helper.fun(FPR, TPR, risk, ipcwCases, ipcwControls, N, time, times[1], status, cens.model, nth.times[1], conservative, MC, cutpoints, se.fit), by = list(model, times)])
  } else if (ROC == TRUE) {
    if (is.null(breaks)) {
      output <- list(ROC = aucDT[nodups, c("model", "times", "risk", "TPR", "FPR"), with = FALSE])
    } else {
      breaks <- sort(breaks, decreasing = TRUE)
      helper.fun <- function(FPR, TPR, risk, breaks) {
        indeces <- sindex(risk, breaks, comp = "greater", FALSE)
        data.table(risk = breaks, TPR = c(rep(0, sum(indeces == 0)), TPR[indeces[indeces != 0]]), FPR = c(rep(0, sum(indeces == 0)), FPR[indeces[indeces != 0]]))
      }
      output <- list(ROC = aucDT[, helper.fun(FPR, TPR, risk, breaks = breaks), by = list(model, times)])
    }
  } else {
    output <- NULL
  }

  data.table::setkey(score, model, times)
  aucDT <- merge(score, aucDT, all = TRUE)
  if (se.fit[[1]] == 1L || multi.split.test[[1]] == TRUE) {
    aucDT[, nth.times := as.numeric(factor(times))]
    ## compute influence function
    ## data.table::setorder(aucDT,model,times,time,-status)
    data.table::setorder(aucDT, model, times, ID)
    aucDT[, IF.AUC := getInfluenceCurve.AUC(
      times[1],
      time,
      status,
      WTi,
      Wt,
      risk,
      MC,
      AUC[1],
      nth.times[1],
      conservative[[1]],
      cens.model
    ), by = list(model, times)]
    se.score <- aucDT[, list(se = sd(IF.AUC) / sqrt(N)), by = list(model, times)]

    data.table::setkey(se.score, model, times)
    score <- score[se.score]
    if (se.fit == 1L) {
      score[, lower := pmax(0, AUC - qnorm(1 - alpha / 2) * se)]
      score[, upper := pmin(1, AUC + qnorm(1 - alpha / 2) * se)]
    } else {
      score[, se := NULL]
    }
    data.table::setkey(aucDT, model, times)
    aucDT <- aucDT[score]
    if (keep.vcov[[1]] == TRUE) {
      output <- c(output, list(vcov = getVcov(aucDT, "IF.AUC", times = TRUE)))
    }
    if (keep.iid[[1]] == TRUE && se.fit[[1]] == TRUE) {
      output <- c(
        output,
        list(iid.decomp = aucDT[, data.table::data.table(ID, model, times, IF.AUC)])
      )
    }
  }
  ## add score to object
  output <- c(list(score = score), output)
  if (length(dolist) > 0) {
    if (se.fit[[1]] == TRUE || multi.split.test[[1]] == TRUE) {
      contrasts.AUC <- aucDT[, getComparisons(data.table(x = AUC, IF = IF.AUC, model = model),
        NF = NF,
        N = N,
        alpha = alpha,
        dolist = dolist, multi.split.test = multi.split.test, se.fit = se.fit
      ), by = list(times)]
    } else {
      contrasts.AUC <- score[, getComparisons(data.table(x = AUC, model = model),
        NF = NF,
        N = N,
        alpha = alpha,
        dolist = dolist,
        multi.split.test = FALSE,
        se.fit = FALSE
      ), by = list(times)]
    }
    setnames(contrasts.AUC, "delta", "delta.AUC")
    output <- c(list(score = score, contrasts = contrasts.AUC), output)
  }
  output
}

#----------------------------------------------------------------------
### AUC.survival.R ends here
