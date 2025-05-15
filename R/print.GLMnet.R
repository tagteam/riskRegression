### print.GLMnet.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May  9 2025 (06:14) 
## Version: 
## Last-Updated: May 15 2025 (06:57) 
##           By: Thomas Alexander Gerds
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Print of a glmnet regression model
#' 
#' Print of a penalized regression model which was fitted with \link{GLMnet}.
#' @param x Object obtained with GLMnet
#' @param ... Passed to print
#' 
#' @method print GLMnet
#' @export
print.GLMnet <- function(x,...){
    cat(paste0(switch(x$fit$call$family, "binomial" = "Penalized logistic",
                      "cox" = "Penalized Cox","Penalized linear")," regression: "),
        "glmnet object fitted via formula interface GLMnet.\n",
        "The fitted glmnet object is storted at x$fit.\n",
        if(x$cv){
            paste0("The penalty parameter lambda was selected via ",x$fit$call$nfolds," cross-validation with loss function ",x$fit$call$type.measure,": lambda=")
        }else{
            if(x$selector == "undersmooth"){
                paste0("The largest possible penalty parameter lambda was selected such that the model converged: lambda=")
            }else{
                paste0("A prespecified penalty parameter lambda value was selected: lambda=")
            }
        },x$selected.lambda,"\n\nThe regression coefficients at the selected lambda value:\n\n",sep = "")
    print(x$selected.beta)
}


######################################################################
### print.GLMnet.R ends here
