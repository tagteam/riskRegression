print.subjectWeights <- function(x,digits=3){
  cat("\nEstimated inverse of the probability of censoring weights (subjectWeights)\n\n")
  method=switch(x$method,
    "cox"="Cox regression",
    "marginal"="Kaplan-Meier",
    "km"="Kaplan-Meier",
    "nonpar"="Stratified Kaplan-Meier",
    "aalen"="Additive Aalen regression",
    "none"="No weighting",
    "Unknown")
  cat("Method for estimation: ", method,"\n")
  cat("Handler function: ",paste(as.character(x$fit$call[1]),"()",sep=""),"\n\n")
  if(x$method %in% c("km","marginal","none")){
    if (x$lag==1)
    gstring <- "G(T_i-)"
    else
    gstring <- "G(T_i)"      
  }else{
    if (x$lag==1)
      gstring <- "G(T_i-|X_i)"
    else
      gstring <- "G(T_i|X_i)"
  }
  cat("Summary of the weights ",gstring,":\n\n")
  print(summary(x$weights),digits=digits,quote=FALSE)
  invisible(x)
}
