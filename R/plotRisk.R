### plotRisk.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Mar 13 2017 (16:53) 
## Version: 
## Last-Updated: Jul  6 2022 (17:38) 
##           By: Thomas Alexander Gerds
##     Update #: 167
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' plot predicted risks 
##'
##' Two rival prediction models are applied to the same data.
##' @title plot predicted risks
##' @param x Object obtained with function \code{Score}
##' @param models Choice of two models to plot. The predicted risks of
##'     the first (second) are shown along the x-axis (y-axis).
##' @param times Time point specifying the prediction horizon.
##' @param xlim x-axis limits
##' @param ylim y-axis limits
##' @param xlab x-axis labels
##' @param ylab y-axis labels
##' @param col Colors used according to the outcome.
##' binary outcome (two colors: no event, event),
##' survival outcome (three colors: censored, event, no event)
##' competing risk outcome (4 or more colors: event, competing risk 1, ..., competing risk k, censored, no event)
##' @param pch Symbols used according to the outcome
##' binary outcome (two symbols: no event, event),
##' survival outcome (three symbols: censored, event, no event)
##' competing risk outcome (4 or more symbols: event, competing risk 1, ..., competing risk k, censored, no event)
##' @param cex point size
##' @param preclipse Value between 0 and 1 defining the preclipse area
##' @param preclipse.shade Logical. If \code{TRUE} shade the area of clinically meaningful change.
##' @param ... Used to control the subroutines: plot, axis, lines,
##'     barplot, legend. See \code{\link{SmartControl}}.
##' @return a nice graph
##' @examples
##' library(prodlim)
##' ## uncensored
##' set.seed(10)
##' learndat = sampleData(40,outcome="binary")
##' testdat = sampleData(40,outcome="binary")
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family="binomial")
##' lr2 = glm(Y~X3+X5+X6,data=learndat,family="binomial")
##' xb=Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,
##'          data=testdat,summary="risks",null.model=0L)
##' plotRisk(xb)
##' ## survival
##' library(survival)
##' set.seed(10)
##' learndat = sampleData(40,outcome="survival")
##' testdat = sampleData(40,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=learndat,x=TRUE)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=learndat,x=TRUE)
##' xs=Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),formula=Surv(time,event)~1,
##'          data=testdat,summary="risks",null.model=0L,times=c(3,5,6))
##' plotRisk(xs,times=5)
##' ## competing risk
##' \dontrun{
##' library(prodlim)
##' library(survival)
##' set.seed(8)
##' learndat = sampleData(80,outcome="competing.risk")
##' testdat = sampleData(140,outcome="competing.risk")
##' m1 = FGR(Hist(time,event)~X2+X7+X9,data=learndat,cause=1)
##' m2 = CSC(Hist(time,event)~X2+X7+X9,data=learndat,cause=1)
##' xcr=Score(list("FGR"=m1,"CSC"=m2),formula=Hist(time,event)~1,
##'          data=testdat,summary="risks",null.model=0L,times=c(3,5))
##' plotRisk(xcr,times=1)
##' }
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
plotRisk <- function(x,
                     models,
                     times,
                     xlim=c(0,1),
                     ylim=c(0,1),
                     xlab,
                     ylab,
                     col,
                     pch,
                     cex=1,
                     preclipse=0,
                     preclipse.shade=FALSE,
                     ...){
    model=ReSpOnSe=risk=status=event=NULL
    if (is.null(x$risks$score)) stop("No predicted risks in object. You should set summary='risks' when calling Score.")
    if (!is.null(x$null.model)){
        pframe <- x$risks$score[model!=x$null.model]
        pframe[,model:=factor(model)]
    } else{
        pframe <- x$risks$score
    }
    if (x$response.type!="binary"){
        if (missing(times)){
            tp <- max(pframe[["times"]])
            if (length(unique(pframe$times))>1)
                warning("Time point not specified, use max of the available times: ",tp)
        } else{ ## can only do one time point
            tp <- times[[1]]
            if (!(tp%in%unique(pframe$times)))
                stop(paste0("Requested time ",times[[1]]," is not in object"))
        }
        pframe <- pframe[times==tp]
    }else tp <- NULL
    if (missing(models)){
        models <- pframe[,levels(model)[1:2]]
    }else{
        if (is.na(models[2])) stop("Need two models for a scatterplot of predicted risks")
        fitted.models <- x$models
        if (is.numeric(models)) {
            if (any(notfitted <- !(models%in%fitted.models))){
                stop(paste("Cannot find all models. Fitted models are:\n",paste(paste0(x$models,": ",names(x$models)),collapse="\n")))
            } else{
                models <- names(fitted.models[models])
            }
        }else{
            if (any(notfitted <- !(models%in%names(fitted.models)))){
                stop(paste("Cannot find all models. Fitted models are:\n",paste(paste0(x$models,": ",names(x$models)),collapse="\n")))
            }
        }
    }
    pframe <- pframe[model%in%models]
    modelnames <- models
    pframe[,model:=factor(model,levels=models)]
    data.table::setkey(pframe,model)
    if (x$response.type=="binary"){
        R <- pframe[model==modelnames[1],ReSpOnSe]
    }
    else{
        if  (x$response.type=="survival"){
            R <- pframe[model==modelnames[1],
            {
                r <- status # 0,1
                r[time>times] <- 2
                r
            }]
        }else{ ## competing risks
            R <- pframe[model==modelnames[1],
            {
                # respect the event of interest
                r <- rep("competing.risk",.N)
                r[event == x$cause] <- "event"
                r[status == 0] <- "censored"
                r[time>times] <- "event-free"
                r <- as.numeric(as.character(factor(r,levels = c("event","competing.risk","event-free","censored"),labels = c(0,1,2,3))))
                r
            }]
        }
    }
    ## pframe[,model:=factor(model)]
    m1 <- levels(pframe$model)[[1]]
    m2 <- levels(pframe$model)[[2]]
    if (missing(xlab)) xlab <- paste0("Risk (%, ",modelnames[1],")")
    if (missing(ylab)) ylab <- paste0("Risk (%, ",modelnames[2],")")
    if (x$response.type=="binary"){
        ppframe <- data.table::dcast(pframe,ID~model,value.var="risk")
    }else{
        ppframe <- data.table::dcast(pframe,times+ID~model,value.var="risk")
    }
    pred.m1 <- ppframe[[m1]]
    pred.m2 <- ppframe[[m2]]
    if (preclipse>0){
        diff <- abs(pred.m1-pred.m2)
        which <- diff>quantile(diff,preclipse)
        message(paste(sprintf("\nShowing absolute differences > %1.2f",
                              100*quantile(diff,preclipse)),"%"))
        pred.m1 <- pred.m1[which]
        pred.m2 <- pred.m2[which]
        R <- R[which]
    }
    nR <- length(unique(R))
    if (x$response.type=="competing.risks") nCR <- nR-3
    Rfactor <- factor(R,
                      levels=sort(unique(R)),
                      labels=switch(x$response.type,
                                    "binary"={c("event-free","event")},
                                    "survival"={c("right-censored","event","event-free")},
                                    "competing.risks"={if (nR==4)
                                                           c("event","competing-risk","right-censored","event-free")
                                                       else
                                                           c("event",paste0("competing-risk",1:nCR),"right-censored","event-free")}))
    if (missing(col)){
        colcode <- switch(x$response.type,
                          "binary"={c("#52da58","red")},
                          "survival"={c("gray55","red","#52da58")},
                          "competing.risks"={c("red",rep("purple",nCR),"gray55","#52da58")})
    }else{
        if (nR!=length(col)){
            warning(paste0("Need exactly ",nR," colors (argument col) in the following order: ",paste0(levels(Rfactor),collapse=", ")))
        }
        colcode <- col
    }
    if (missing(pch)){
        pchcode <- switch(x$response.type,
                          "binary"={c(1,16)},
                          "survival"={c(8,16,1)},
                          "competing.risks"={c(16,rep(17,nCR),8,1)})
    } else{
        if (nR!=length(pch)){
            warning(paste0("Need exactly ",nR," symbols (argument pch) in the following order: ",paste0(labels(Rfactor),
                                                                                                        collapse=", ")))
        }
        pchcode <- pch
    }
    pch=as.numeric(as.character(factor(Rfactor,labels=pchcode)))
    col=as.character(factor(Rfactor,labels=colcode))
    # {{{ smart argument control
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,ylab=ylab,xlab=xlab)
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4))
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),mgp=c(4,1,0))
    this.legend <- switch(x$response.type,"binary"={paste0(c("No event","Event")," (n=",c(sum(R==0),sum(R==1)),")")},
                          "survival"={paste0(c("Censored","Event","No event")," (n=",c(sum(R==0),sum(R==1),sum(R==2)),")")},
                          "competing.risks"={
                              if (nCR==1){
                                  paste0(c("Event","Competing risk","Censored","No event")," (n=",c(sum(R==0),sum(R==1),sum(R==2),sum(R==3)),")")
                              }else{
                                  paste0(c("Event",paste0("Competing risk ",1:nCR),"Censored","No event")," (n=",sapply(1:nR,function(r){sum(R==r)}),")")
                              }
                          })
    legend.DefaultArgs <- list(legend=this.legend,
                               pch=pchcode,
                               col=colcode,
                               cex=cex,
                               bty="n",
                               y.intersp=1.3,
                               x="topleft")
    points.DefaultArgs <- list(x=pred.m1,y=pred.m2,pch=pch,cex=cex,col=col)
    abline.DefaultArgs <- list(a=0,b=1,lwd=1,col="gray66")
    control <- prodlim::SmartControl(call= list(...),
                                     keys=c("plot","points","legend","axis1","axis2","abline"),
                                     ignore=NULL,
                                     ignore.case=TRUE,
                                     defaults=list("plot"=plot.DefaultArgs,
                                                   "points"=points.DefaultArgs,
                                                   "abline"=abline.DefaultArgs,
                                                   "legend"=legend.DefaultArgs,
                                                   "axis1"=axis1.DefaultArgs,
                                                   "axis2"=axis2.DefaultArgs),
                                     forced=list("plot"=list(axes=FALSE),"legend"=list(col=colcode)),
                                     verbose=TRUE)
    # }}}
    do.call("plot",control$plot)
    if (preclipse.shade==TRUE){
        ## rect(xleft=0,xright=1,ybottom=0,ytop=1,col="white",border=0)
        plotrix::draw.ellipse(x=0.5,y=.5, a=.75, b=.1, border=0, angle=c(45), lty=3,col="gray88")
    }
    do.call("abline",control$abline)
    do.call("points",control$points)
    control$axis2$labels <- paste(100*control$axis2$at,"%")
    control$axis1$labels <- paste(100*control$axis1$at,"%")
    do.call("axis",control$axis1)
    do.call("axis",control$axis2)
    do.call("legend",control$legend)
}


######################################################################
### plotRisk.R ends here
