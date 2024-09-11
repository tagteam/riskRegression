##'
##' @title Export a \code{synth} object.
##' @description Function for exporting objects of type \code{synth} via \code{dput}.
##' @param object An object of type \code{synth}.
##' @param file File to print to.
##' @return Either a string or nothing if printed to a file.
##' @seealso lvm
##' @examples
##' # See the documentation for the \code{synthesize} function.
##'
##' @export
saveSynth <- function(object, file = ""){
  ## should probably write try around some of the eval(parse(text = ...)) calls at some point
  zzz <- NULL
  if (missing(object)) stop("No object specified")
  if (!inherits(object,"synth")) stop("Object is not of class 'synth'")
  ## handle various function objects in distribution attribute
  ## we replace the corresponding parts of the function with the actual values as default values
  for (l in names(object$lava.object$attributes$distribution)){
    if ((at<-attributes(object$lava.object$attributes$distribution[[l]])$family$family) %in% c("binomial", "gaussian", "weibull")){
      b_attr <- attributes(object$lava.object$attributes$distribution[[l]])
      if (at == "binomial"){
        b_fun <- function (n, mu, var, ...) 
        {
          rbinom(n, 1, stats::binomial()$linkinv(mu))
        }
      } else if (at == "gaussian"){
        b_fun <- function (n, mu, var, ...) 
        {
          stats::rnorm(n, mu, sqrt(var))
        }
      } else if (at == "weibull"){
        b_fun <- object$lava.object$attributes$distribution[[l]]
        attributes(b_fun) <- NULL
        dep_fun <- deparse(b_fun)
        dep_fun <- gsub("\\b(scale)\\b", b_attr$family$par[2], dep_fun)
        dep_fun <- gsub("\\b(shape)\\b", b_attr$family$par[1], dep_fun)
        b_fun <- eval(parse(text=dep_fun))
      }
      b_attr$srcref <- NULL
      object$lava.object$attributes$distribution[[l]] <- b_fun
      attributes(object$lava.object$attributes$distribution[[l]]) <- b_attr
    } else {
      stop("Distribution not implemented")
    }
  }
  
  ## case: categorical
  bb_fun <- function (y, p, idx = zzz, ...) 
  {
    if (length(p[idx])){
      theta <- p[idx]
      v <- theta[1]
      breaks <- c(-Inf, cumsum(c(v,exp(theta[seq(length(theta)-1L)+1L]))), Inf) 
      as.numeric(cut(y, breaks = breaks)) - 1
    } else {
      stop("error")
    }
  }
  deparsed_bb_fun <- deparse(bb_fun)
  
  for (l in names(object$lava.object$constrainY)){
    ll <- object$lava.object$attributes$ordinalparname[[l]]
    dep_fun <- gsub("\\b(zzz)\\b", paste0(deparse(ll),collapse=""), deparsed_bb_fun)
    object$lava.object$constrainY[[l]]$fun <- eval(parse(text=dep_fun))
  }
  for (l in names(object$lava.object$attributes$transform)){
    curr <- object$lava.object$attributes$transform[[l]]
    bbb_fun <- curr$fun
    dep_fun <- deparse(bbb_fun)
    dep_fun <- gsub("\\b(l)\\b", paste0("\"", gsub(curr$x,"",l,fixed=TRUE),"\""), dep_fun)
    object$lava.object$attributes$transform[[l]]$fun <- eval(parse(text=dep_fun)) 
  }
  
  ## survival/comp.risk
  if (length(object$lava.object$attributes$multitransform) > 0) {
    bbbb_fun <- function(z, events = zzz, ...) {
      idx <- apply(z,1,which.min)
      cbind(z[cbind(seq(NROW(z)),idx)],events[idx])
    }
    dep_fun <- deparse(bbbb_fun)
    dep_fun <- gsub("\\b(zzz)\\b", deparse(object$lava.object$attributes$eventHistory$time$events), dep_fun)
    object$lava.object$attributes$multitransform[[1]]$fun <- eval(parse(text=dep_fun)) 
  }
  dput(object, file = file,
       control = c("keepNA", "keepInteger", "niceNames","showAttributes", "useSource","quoteExpressions"))
}
