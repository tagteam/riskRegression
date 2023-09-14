#' Backward variable selection in the Cox regression model
#' 
#' This is a wrapper function which first selects variables in the Cox
#' regression model using \code{fastbw} from the \code{rms} package and then
#' returns a fitted Cox regression model with the selected variables.
#' 
#' This function first calls \code{cph} then \code{fastbw} and finally
#' \code{cph} again.
#' 
#' @param formula A formula object with a \code{Surv} object on the left-hand
#' side and all the variables on the right-hand side.
#' @param data Name of an data frame containing all needed variables.
#' @param rule The method for selecting variables. See \code{\link{fastbw}} for
#' details.
#' @references Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23. URL
#' http://www.jstatsoft.org/v50/i11/.
#' @keywords survival
#' @examples
#' library(survival)
#' set.seed(74)
#' d <- sampleData(89,outcome="survival")
#' f <- selectCox(Surv(time,event)~X1+X2+X3+X4+X6+X7+X8+X9, data=d)
#' 
#' @export 
selectCox <- function(formula,data,rule="aic"){
    fit <- rms::cph(formula, data, surv=TRUE,x=TRUE,y=TRUE)
    bwfit <- rms::fastbw(fit,rule=rule)
    ## print(bwfit$names.kept)
    if (length(bwfit$names.kept)==0){
        newform <- update(formula,".~1")
        newfit <- prodlim::prodlim(newform,data=data)
    }
    else{
        newform <- update(formula,paste(".~",paste(bwfit$names.kept,collapse="+")))
        ## reformulate(bwfit$names.kept, formula[[2]])
        newfit <- rms::cph(newform,data, surv=TRUE,x=TRUE,y=TRUE)
    }
    out <- list(fit=newfit,In=bwfit$names.kept)
    out$call <- match.call()
    class(out) <- "selectCox"
    out
}


