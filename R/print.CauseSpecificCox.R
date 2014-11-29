#' @S3method print CauseSpecificCox
print.CauseSpecificCox <- function(x,...){
    print(x$call)
    print(x$response)
    if (x$survtype=="hazard"){
        nix <- lapply(1:length(x$causes),function(c){
            cat("\n\n----------> Cause: ",x$causes[c],"\n\n")
            xc <- x$models[[c]]
            xc$call$data <- NULL
            print(summary(xc),...)
        })
    }
    else{ # survtype=="survival"
        cat("\n\n----------> Cause: ",x$theCause,"\n\n")
        x1 <- x$models[[1]]
        x1$call$data <- NULL
        print(summary(x1),...)
        cat("\n\n----------> Event-free survival:\n\n")
        x2 <- x$models[[2]]
        x2$call$data <- NULL
        print(summary(x2),...)
    }
}
