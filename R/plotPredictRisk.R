#' Plotting predicted risks curves.
#' 
#' Time-dependent event risk predictions.
#' 
#' Arguments for the invoked functions \code{legend} and \code{axis} can be
#' specified as \code{legend.lty=2}. The specification is not case sensitive,
#' thus \code{Legend.lty=2} or \code{LEGEND.lty=2} will have the same effect.
#' The function \code{axis} is called twice, and arguments of the form
#' \code{axis1.labels}, \code{axis1.at} are used for the time axis whereas
#' \code{axis2.pos}, \code{axis1.labels}, etc., are used for the y-axis.
#' 
#' These arguments are processed via \code{\dots{}} of
#' \code{plotPredictRisk} and inside by using the function
#' \code{SmartControl}.
#' 
#' @param x Object specifying an event risk prediction model.
#' @param newdata A data frame with the same variable names as those that were
#' used to fit the model \code{x}.
#' @param times Vector of times at which to return the estimated probabilities.
#' @param cause Show predicted risk of events of this cause
#' @param xlim Plotting range on the x-axis.
#' @param ylim Plotting range on the y-axis.
#' @param xlab Label given to the x-axis.
#' @param ylab Label given to the y-axis.
#' @param axes Logical. If \code{FALSE} no axes are drawn.
#' @param col Vector of colors given to the survival curve.
#' @param density Densitiy of the color -- useful for showing many
#' (overlapping) curves.
#' @param lty Vector of lty's given to the survival curve.
#' @param lwd Vector of lwd's given to the survival curve.
#' @param add Logical. If \code{TRUE} only lines are added to an existing
#' device
#' @param legend Logical. If TRUE a legend is plotted by calling the function
#' legend.  Optional arguments of the function \code{legend} can be given in
#' the form \code{legend.x=val} where x is the name of the argument and val the
#' desired value. See also Details.
#' @param percent Logical. If \code{TRUE} the y-axis is labeled in percent.
#' @param \dots Parameters that are filtered by \code{\link[prodlim]{SmartControl}} and
#' then passed to the functions: \code{\link{plot}}, \code{\link{axis}},
#' \code{\link{legend}}.
#' @return The (invisible) object.
#' @author Ulla B. Mogensen and Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{plotRisk}}
#' @references Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23. URL
#' http://www.jstatsoft.org/v50/i11/.
#' @keywords survival
#' @examples
#' library(survival)
#' # generate survival data
#' # no effect
#' set.seed(8)
#' d <- sampleData(80,outcome="survival",formula = ~f(X6, 0) + f(X7, 0))
#' d[,table(event)]
#' f <- coxph(Surv(time,event)~X6+X7,data=d,x=1)
#' plotPredictRisk(f)
#'
#' # large effect
#' set.seed(8)
#' d <- sampleData(80,outcome="survival",formula = ~f(X6, 0.1) + f(X7, -0.1))
#' d[,table(event)]
#' f <- coxph(Surv(time,event)~X6+X7,data=d,x=1)
#' plotPredictRisk(f)
#' 
#' # generate competing risk data
#' # small effect
#' set.seed(8)
#' d <- sampleData(40,formula = ~f(X6, 0.01) + f(X7, -0.01))
#' d[,table(event)]
#' f <- CSC(Hist(time,event)~X5+X6,data=d)
#' plotPredictRisk(f)
#'
#' # large effect
#' set.seed(8)
#' d <- sampleData(40,formula = ~f(X6, 0.1) + f(X7, -0.1))
#' d[,table(event)]
#' f <- CSC(Hist(time,event)~X5+X6,data=d)
#' plotPredictRisk(f)
#' @export
plotPredictRisk <- function(x,
                            newdata,
                            times,
                            cause=1,
                            xlim,
                            ylim,
                            xlab,
                            ylab,
                            axes=TRUE,
                            col,
                            density,
                            lty,
                            lwd,
                            add=FALSE,
                            legend=TRUE,
                            percent=FALSE,
                            ...){
    # {{{ call argument

    allArgs <- match.call()
    # }}}
    # {{{ find times
    
    if(missing(times)){
        # formula
        formula <- eval(x$call$formula)
        if (!inherits(x=formula,what="formula"))
            stop("Argument formula is missing.")
        # find data
        data <- eval(x$call$data)
        # extract response
        m <- model.frame(formula,data,na.action=na.fail)
        response <- model.response(m)
        # ordering time
        if (!(match("time",colnames(response),nomatch=FALSE))){
            stop("Cannot find time variable in response. Perhaps this is not a time-to-event model?")
        }
        neworder <- order(response[,"time"],-response[,"status"])
        response <- response[neworder,,drop=FALSE]
        times <- response[,"time"]
        # unique event times
        times <- unique(times)
    }

    # }}}
    # {{{ newdata
    if(missing(newdata)){
        newdata <- eval(x$call$data)
    }
    ## stop("newdata argument is missing")

    # }}}
    # {{{ xlim, ylim

    if (missing(xlim)) xlim <- c(0, max(times))
    at <- times <= xlim[2]
    orig.X <- times[at]
    X <- times[at]
  
    # }}}  
    # {{{ predict newdata at times
    y <- predictRisk(object=x,newdata=newdata,times=orig.X,cause=cause)
    # }}}
    # {{{  plot arguments

    nlines <- NROW(y)
  
    if (missing(ylab)) ylab <- "Absolute risk of event"
    if (missing(xlab)) xlab <- "Time"
    if (missing(ylim)) ylim <- c(0, 1)
    if (missing(lwd)) lwd <- rep(3,nlines)
    if (missing(col)) col <- rep(1,nlines)
    if (missing(density)){
        if (nlines>5){
            density <- pmax(33,100-nlines)
        }
        else density <- 100
    }
    ## print(density)
    if (density<100){
        col <- sapply(col,function(coli){
            ccrgb=as.list(col2rgb(coli,alpha=TRUE))
            names(ccrgb) <- c("red","green","blue","alpha")
            ccrgb$alpha=density
            cc=do.call("rgb",c(ccrgb,list(max=255)))
        })
    }
    if (missing(lty)) lty <- rep(1, nlines)
    if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
    if (length(lty) < nlines) lty <- rep(lty, nlines)
    if (length(col) < nlines) col <- rep(col, nlines)
  
  axis1.DefaultArgs <- list()
  
  axis2.DefaultArgs <- list(at=seq(0,1,.25),las=2,labels=paste(100*seq(0,1,.25),"%"))
         
  plot.DefaultArgs <- list(x=0,
                           y=0,
                           type = "n",
                           ylim = ylim,
                           xlim = xlim,
                           xlab = xlab,
                           ylab = ylab)
  
  
  
  legend.DefaultArgs <- list(legend=rownames(y),
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             x="topright")
  
  # }}}
  # {{{ smart controls
    
  if (match("legend.args",names(args),nomatch=FALSE)){
    legend.DefaultArgs <- c(args[[match("legend.args",names(args),nomatch=FALSE)]],legend.DefaultArgs)
    legend.DefaultArgs <- legend.DefaultArgs[!duplicated(names(legend.DefaultArgs))]
  }
  smartA <- prodlim::SmartControl(call=list(...),
                                   keys=c("plot","legend","axis1","axis2"),
                                   ignore=c("x", "newdata", "times", "xlim","ylim","xlab","ylab","col","lty","lwd","add","legend","percent","axes","legend.args"),
                                   defaults=list("plot"=plot.DefaultArgs,
                                     "legend"= legend.DefaultArgs,
                                     "axis1"=axis1.DefaultArgs,
                                     "axis2"=axis2.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),
                                     "axis1"=list(side=1),
                                     "axis2"=list(side=2)),
                                   verbose=TRUE)

  # }}} 
  # {{{ empty plot
  
  if (!add) {
    do.call("plot",smartA$plot)
    if (axes){
      do.call("axis",smartA$axis1)
      if (percent & is.null(smartA$axis1$labels))
        smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
      do.call("axis",smartA$axis2)
    }
  }

  # }}}
  # {{{ adding lines
  nix <- lapply(1:nlines, function(s) {
    lines(x = X, y = y[s,], type = "s", col = col[s], lty = lty[s], lwd = lwd[s])
  })

  # }}}
  # {{{ legend

  if(legend==TRUE && !add && !is.null(rownames(y))){
    save.xpd <- par()$xpd
    do.call("legend",smartA$legend)
    par(xpd=save.xpd)
  }
  
  # }}}
  invisible(x)
}

