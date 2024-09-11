### tidy.ate.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 24 2023 (18:04) 
## Version: 
## Last-Updated: okt 24 2023 (18:21) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## ##' @description EXPERIMENTAL (inspired from broom:::tidy.glm)
## tidy.ate <- function(x, conf.int = FALSE, conf.level = 0.95, ...){
##     broom:::warn_on_subclass(x, "tidy")
##     x.confint <- stats::confint(x, level = conf.level)

##     x.meanRisk <- data.frame(estimate = x.confint$meanRisk$estimate,
##                              std.error = x.confint$meanRisk$se,
##                              statistic = NA,
##                              p.value = NA,
##                              lower = x.confint$meanRisk$lower,
##                              upper = x.confint$meanRisk$upper,
##                              check.names = FALSE)
##     rownames(x.meanRisk) <- paste0("meanRisk:",x.confint$meanRisk$treatment,"(t=",x.confint$meanRisk$time,")")
##     colnames(x.meanRisk)[5:6] <- c(paste0(100*(1-conf.level)/2," %"),paste0(100*(1-(1-conf.level)/2)," %"))
##     x.diffRisk <- data.frame(estimate = x.confint$diffRisk$estimate,
##                              std.error = x.confint$diffRisk$se,
##                              statistic = x.confint$diffRisk$estimate/x.confint$diffRisk$se,
##                              p.value = x.confint$diffRisk$p.value,
##                              lower = x.confint$diffRisk$lower,
##                              upper = x.confint$diffRisk$upper,
##                              check.names = FALSE)
##     rownames(x.diffRisk) <- paste0("diffRisk:",x.confint$diffRisk$B,"-",x.confint$diffRisk$A,"(t=",x.confint$diffRisk$time,")")
##     colnames(x.diffRisk)[5:6] <- c(paste0(100*(1-conf.level)/2," %"),paste0(100*(1-(1-conf.level)/2)," %"))
##     x.ratioRisk <- data.frame(estimate = x.confint$ratioRisk$estimate,
##                               std.error = x.confint$ratioRisk$se,
##                               statistic = x.confint$ratioRisk$estimate/x.confint$ratioRisk$se,
##                               p.value = x.confint$ratioRisk$p.value,
##                               lower = x.confint$ratioRisk$lower,
##                               upper = x.confint$ratioRisk$upper,
##                               check.names = FALSE)
##     rownames(x.ratioRisk) <- paste0("ratioRisk:",x.confint$ratioRisk$B,"/",x.confint$ratioRisk$A,"(t=",x.confint$ratioRisk$time,")")
##     colnames(x.ratioRisk)[5:6] <- c(paste0(100*(1-conf.level)/2," %"),paste0(100*(1-(1-conf.level)/2)," %"))
##     x.allRisk <- rbind(x.meanRisk, x.diffRisk, x.ratioRisk)

##     if (conf.int) {
##         ret <- tibble:::as_tibble(x.allRisk, rownames = "term")
##     }else{
##         ret <- tibble:::as_tibble(x.allRisk[,1:4], rownames = "term")
##     }
##     ret
## }

##----------------------------------------------------------------------
### tidy.ate.R ends here
