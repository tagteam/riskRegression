#' Paquid sample
#'
#' PAQUID is a prospective cohort study initiated in 1988 in South Western
#' France to explore functional and cerebral ageing. This sample includes
#' n=2561 subjects. Data contains a time-to-event, a type of event and
#' two cognitive scores measured at baseline.
#'
#'
#' @name Paquid
#' @docType data
#' @format A data frame with 2561 observations on the following 4 variables.
#' \describe{
#' \item{\code{time}}{the time-to-event (in years).}
#' \item{\code{status}}{the type of event \code{0} = censored, \code{1} = dementia onset and \code{2} = death without dementia.}
#' \item{\code{DSST}}{score at the Digit Symbol Substitution Score Test. This test explores attention and psychomotor speed.}
#' \item{\code{MMSE}}{score at the Mini Mental State Examination. This test is often used as an index of global cognitive performance.}
#' }
#' @references
#'
#' Dartigues, J., Gagnon, M., Barberger-Gateau, P., Letenneur, L., Commenges, D.,
#' Sauvel, C., Michel, P., and Salamon, R. (1992). The paquid epidemiological program
#' on brain ageing. Neuroepidemiology, 11(1):14--18.
#'
#'
#' Blanche, P., Dartigues, J. F., & Jacqmin-Gadda, H. (2013). Estimating and
#' comparing time-dependent areas under receiver operating characteristic curves
#' for censored event times with competing risks. Statistics in Medicine, 32(30),
#' 5381-5397.
#'
#' @source
#' The data have been first made publicly available via the package timeROC.
#' @keywords datasets
#' @examples
#' data(Paquid)
NULL

#' Malignant melanoma data
#'
#' In the period 1962-77, 205 patients with malignant melanoma (cancer of the
#' skin) had a radical operation performed at Odense University Hospital,
#' Denmark. All patients were followed until the end of 1977 by which time 134
#' were still alive while 71 had died (of out whom 57 had died from cancer and
#' 14 from other causes).
#'
#' The object of the study was to assess the effect of risk factors on
#' survival. Among such risk factors were the sex and age of the patients and
#' the histological variables tumor thickness and ulceration (absent vs.
#' present).
#
#'
#' @name Melanoma
#' @docType data
#' @format A data frame with 205 observations on the following 12 variables.
#' \describe{
#' \item{time}{ time in days from operation}
#' \item{status}{a numeric with values \code{0=censored} \code{1=death.malignant.melanoma} \code{2=death.other.causes}}
#' \item{event}{a factor with levels \code{censored} \code{death.malignant.melanoma} \code{death.other.causes}}
#' \item{invasion}{a factor with levels \code{level.0}, \code{level.1}, \code{level.2}}
#' \item{ici}{inflammatory cell infiltration (IFI): 0, 1, 2 or 3}
#' \item{epicel}{a factor with levels \code{not present} \code{present}}
#' \item{ulcer}{a factor with levels \code{not present} \code{present}}
#' \item{thick}{tumour thickness (in 1/100 mm)}
#' \item{sex}{a factor with levels \code{Female} \code{Male}}
#' \item{age}{age at operation (years)}
#' \item{logthick}{tumour thickness on log-scale}
#' }
#' @references Regression with linear predictors (2010)
#'
#' Andersen, P.K. and Skovgaard, L.T.
#'
#' Springer Verlag
#'
#' @keywords datasets
##' @examples
##'
##' data(Melanoma)
NULL

#' @docType package
#' @name riskRegression
#' @useDynLib riskRegression, .registration=TRUE
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom grDevices rainbow
#' @importFrom rms strat cph
#' @importFrom mets phreg
#' @importFrom ggplot2 autoplot aes aes_string element_blank element_line element_rect geom_errorbar geom_line geom_point geom_ribbon ggplot labs guide_legend guides scale_colour_manual scale_color_continuous scale_fill_manual scale_linetype_manual scale_y_continuous theme theme_bw "%+replace%"  unit xlab  ylab
#' @importFrom survival Surv strata coxph survreg
#' @importFrom lava sim iid information score transform<- exogenous endogenous regression<-
#' @importFrom data.table data.table set dcast setkeyv setDT as.data.table copy data.table is.data.table melt rbindlist setnames setorder setcolorder setkey ":=" ".N" ".SD"
#' @importFrom prodlim Hist dimColor prodlim
#' @importFrom foreach "%dopar%" foreach "%do%"
#' @importFrom cmprsk predict.crr
#' @importFrom timereg comp.risk Event
#' @importFrom prodlim Hist jackknife prodlim sindex
#' @importFrom grDevices col2rgb gray
#' @importFrom graphics bxp abline axis box legend lines mtext par plot points segments text title polygon par boxplot
#' @importFrom utils capture.output find head select.list setTxtProgressBar tail txtProgressBar
#' @importFrom stats confint cov as.formula coef delete.response drop.terms family formula get_all_vars lm glm median model.frame model.matrix model.response na.fail na.omit nobs optim pnorm predict qnorm quantile rbinom reformulate rexp runif sd setNames smooth terms terms.formula time uniroot update update.formula var wilcox.test
NULL

