print.cumincCox <- function(x,...){
  print(x$response)
  cat("\n\nFitted Cox models:\n\n")
  print(lapply(x$models,function(m){
    m
  }))
}
