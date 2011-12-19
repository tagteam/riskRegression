print.CauseSpecificCox <- function(x,...){
  print(x$call)
  print(x$response)
  if (x$survtype=="hazard"){
    nix <- lapply(1:length(x$causes),function(c){
      cat("\n\n----------> Cause: ",x$cause[c],"\n\n")
      print(summary(x$models[[c]]),...)
    })
  }
  else{ # survtype=="survival"
    cat("\n\n----------> Cause: ",x$theCause,"\n\n")
    print(summary(x$models[[1]]),...)
    cat("\n\n----------> Event-free survival:\n\n")
    print(summary(x$models[[2]]),...)
  }
}
