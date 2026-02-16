### IPA.Compute.R ---
#----------------------------------------------------------------------
## Author: Asbjørn Risom
## Created: aug 14 2025
## Version:
## Last-Updated: feb  6 2026 (12:13) 
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

IPACompute <- function(brier.results, 
                       response.type,
                       se.fit,
                       alpha) {
    Brier = model = IF.IPA = Brier.Null = IF.Brier = IF.Brier.Null = . = lower = se = upper = times =  NULL
    score <- brier.results$score
    # FIXME: this is silly but necessary
    if ("IBS" %in% names(score)){
        score[["IBS"]] <- NULL
    }
    # Calculate IPA point estimate
    if (response.type == "binary") {
        byvars <- "model"
        bystat <- NULL
        score[, IPA := 1 - Brier / Brier[model == "Null model"]]
    } else {
        byvars <- c("model","times")
        bystat <- "times"
        score[, IPA := 1 - Brier / Brier[model == "Null model"], by = times]
    }
    if (se.fit == 1L) {
        # Ensure Influence curve is available and has necessary columns
        if (!is.null(brier.results$iid.decomp) && "IF.Brier" %in% names(brier.results$iid.decomp)) {
            # Get brier and IF.Brier for all models and the corresponding null model. Merge to null model to compute the IF for IPA
            model.IF.data <- brier.results$score[,.SD,.SDcols = c(byvars,"Brier")][brier.results$iid.decomp[,.SD,.SDcols = c("riskRegression_ID",byvars,"IF.Brier")],on = byvars]
            null.IF.data <- model.IF.data[model == "Null model"][,-"model"]
            merged.IF.data <- merge(model.IF.data, null.IF.data, by = c("riskRegression_ID", if (response.type != "binary") "times"),suffixes = c("",".Null"))
            # Calculate IF for IPA and corresponding standard errors
            merged.IF.data[, IF.IPA := -(1/Brier.Null) * IF.Brier + (Brier / Brier.Null^2) * IF.Brier.Null]
            ipa.se.dt <- merged.IF.data[, .(se = sd(IF.IPA) / sqrt(.N)), by = c("model", if (response.type != "binary") "times")]
            #Merge and compute confidence intervals for IPA
            score <- score[,-c("se","lower","upper")]
            score <- merge(score, ipa.se.dt, by = c("model", if (response.type != "binary") "times"))
            score[, lower := IPA - qnorm(1 - alpha/2) * se]
            score[, upper := IPA + qnorm(1 - alpha/2) * se]
        } else {
            warning("Brier influence curves not available for IPA confidence interval computation.")
        }
    }
    score <- score[,-"Brier"]
    return(list(score = score))
}
