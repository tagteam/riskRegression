summary.cumIncCox <- function(object, cause){

  n.event <- object$n.event 

  if (missing(cause)) N <- n.event
  else N <- cause

  show=lapply(1:N, function(x){
    summary(object$cause.haz[[x]])
  })


  
}
  
