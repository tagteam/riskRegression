### plotCalibration.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb 23 2017 (11:15) 
## Version: 
## last-updated: Jan 29 2019 (11:13) 
##           By: Thomas Alexander Gerds
##     Update #: 277
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Plot Calibration curve 
##'
##' @title Plot Calibration curve
##' @export
##' @param x Object obtained with function \code{Score}
##' @param models Choice of models to plot
##' @param times Time point specifying the prediction horizon.
##' @param method The method for estimating the calibration curve(s):
##' 
##' \code{"nne"}: The expected event status is obtained in the nearest
##' neighborhood around the predicted event probabilities.
##' 
##' \code{"quantile"}: The expected event status is obtained in groups
##' defined by quantiles of the predicted event probabilities.
##' @param round If \code{TRUE} predicted probabilities are rounded to
##'     two digits before smoothing. This may have a considerable
##'     effect on computing efficiency in large data sets.
##' @param bandwidth The bandwidth for \code{method="nne"}
##' @param q The number of quantiles for \code{method="quantile"} and
##'     \code{bars=TRUE}.
##' @param bars If \code{TRUE}, use barplots to show calibration.
##' @param hanging Barplots only. If \code{TRUE}, hang bars
##'     corresponding to observed frequencies at the value of the
##'     corresponding prediction.
##' @param names Barplots only. Names argument passed to
##'     \code{names.arg} of \code{barplot}.
##' @param pseudo If \code{TRUE} show pseudo values (only for right
##'     censored data).
##' @param rug If \code{TRUE} show rug plot at the predictions
##' @param show.frequencies Barplots only. If \code{TRUE}, show
##'     frequencies above the bars.
##' @param plot If \code{FALSE}, do not plot the results, just return
##'     a plottable object.
##' @param add If \code{TRUE} the line(s) are added to an existing
##'     plot.
##' @param diag If \code{FALSE} no diagonal line is drawn.
##' @param legend Logical. If \code{TRUE} draw legend.
##' @param auc.in.legend Logical. If \code{TRUE} add AUC to legend.
##' @param brier.in.legend Logical. If \code{TRUE} add Brier score to
##'     legend.
##' @param axes If \code{FALSE} no axes are drawn.
##' @param xlim Limits of x-axis.
##' @param ylim Limits of y-axis.
##' @param xlab Label for y-axis.
##' @param ylab Label for x-axis.
##' @param col Vector with colors, one for each element of
##'     object. Passed to \code{\link{lines}}.
##' @param lwd Vector with line widths, one for each element of
##'     object. Passed to \code{\link{lines}}.
##' @param lty lwd Vector with line style, one for each element of
##'     object.  Passed to \code{\link{lines}}.
##' @param pch Passed to \code{\link{lines}}.
##' @param type Passed to \code{\link{lines}}.
##' @param cause For competing risks models, the cause of failure or
##'     event of interest
##' @param percent If TRUE axes labels are multiplied by 100 and thus
##'     interpretable on a percent scale.
##' @param na.action what to do with NA values. Passed to
##'     \code{\link{model.frame}}
##' @param cex Default cex used for legend and labels.
##' @param ... Used to control the subroutines: plot, axis, lines,
##'     barplot, legend, addtable2plot, points (pseudo values), rug. See
##'     \code{\link{SmartControl}}.
##' @examples
##' # binary 
##' db=sampleData(100,outcome="binary")
##' fb1=glm(Y~X1+X5+X7,data=db,family="binomial")
##' fb2=glm(Y~X1+X3+X6+X7,data=db,family="binomial")
##' xb=Score(list(model1=fb1,model2=fb2),Y~1,data=db,
##'           plots="cal")
##' plotCalibration(xb)
##' plotCalibration(xb,bars=TRUE,model="model1")
##' plotCalibration(xb,models=1,bars=TRUE,names.cex=1.3)
##' 
##' # survival
##' library(survival)
##' ds=sampleData(100,outcome="survival")
##' fs1=coxph(Surv(time,event)~X1+X5+X7,data=ds,x=1)
##' fs2=coxph(Surv(time,event)~X1+X3+X6+X7,data=ds,x=1)
##' xs=Score(list(Cox1=fs1,Cox2=fs2),Surv(time,event)~1,data=ds,
##'           plots="cal",metrics=NULL)
##' plotCalibration(xs)
##'
##'
##' # competing risks
##' \dontrun{
##' data(Melanoma)
##' f1 <- CSC(Hist(time,status)~age+sex+epicel+ulcer,data=Melanoma)
##' f2 <- CSC(Hist(time,status)~age+sex+logthick+epicel+ulcer,data=Melanoma)
##' x <- Score(list(model1=f1,model2=f2),Hist(time,status)~1,data=Melanoma,
##'            cause= 2,times=5*365.25,plots="cal")
##' plotCalibration(x)
##' }
##' 
plotCalibration <- function(x,
                            models,
                            times,
                            method="nne",
                            round=TRUE,
                            bandwidth=NULL,
                            q=10,
                            bars=FALSE,
                            hanging=FALSE,
                            names="quantiles",
                            pseudo,
                            rug,
                            show.frequencies=FALSE,
                            plot=TRUE,
                            add=FALSE,
                            diag=!add,
                            legend=!add,
                            auc.in.legend=TRUE,
                            brier.in.legend=TRUE,
                            axes=!add,
                            xlim=c(0,1),
                            ylim=c(0,1),
                            xlab=ifelse(bars,"Risk groups","Predicted risk"),
                            ylab="Observed frequency",
                            col,
                            lwd,
                            lty,
                            pch,
                            type,
                            cause=1,
                            percent=TRUE,
                            na.action=na.fail,
                            cex=1,
                            ...){
    # {{{ plot frame
    model=risk=NULL
    if (missing(pseudo) & missing(rug))
        if (x$cens.type=="rightCensored"){
            showPseudo <- TRUE
            showRug <- FALSE
        } else{
            showPseudo <- FALSE
            showRug <- TRUE
        }
    else{ 
            if (missing(pseudo) || pseudo[[1]]==FALSE) showPseudo <- FALSE else showPseudo <- TRUE
            if (missing(rug) || rug[[1]]==FALSE) showRug <- FALSE else showRug <- TRUE
        }
    pframe <- x$Calibration$plotframe
    if (is.null(pframe))
        stop("Object has no information for calibration plot.\nYou should call the function \"riskRegression::Score\" with plots=\"calibration\".")
    Rvar <- grep("^(ReSpOnSe|pseudovalue)$",names(pframe),value=TRUE)
    if (!missing(models)){
        fitted.models <- pframe[,unique(model)]
        if (all(models%in%fitted.models)){
            pframe <- pframe[model%in%models]
        } else{
            if (all(is.numeric(models)) && (max(models)<=length(fitted.models))){
                models <- fitted.models[models]
                pframe <- pframe[model%in%models]
            }else{
                stop(paste0("The requested models: ",
                            models,
                            "\ndo not match the fitted models: ",
                            paste0(fitted.models,collapse=", ")))
            }
        }
    }
    data.table::setkey(pframe,model)
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
    ## plot(0,0,type="n",ylim = 0:1,xlim = 0:1,axes=FALSE,xlab = xlab,ylab = ylab)
    NF <- pframe[,length(unique(model))]
    if (bars){
        method="quantile"
        if (!(NF==1)) stop(paste0("Barplots work only for one risk prediction model at a time. Provided are ",NF, "models."))
    }

    # }}}
    # {{{ lines 

    if (missing(lwd)) lwd <- rep(3,NF)
    if (missing(col)) {
        if (bars)
            col <- c("grey90","grey30")
        else
            col <- 1:NF
    }
    if (missing(type)){
        if (method=="quantile"){
            type <- rep("b",length.out=NF)
        } else{
            type <- rep("l",length.out=NF)
        }
    }
    if (missing(lty)) lty <- rep(1, length.out=NF)
    if (missing(pch)) pch <- rep(1, length.out=NF)
    if (length(lwd) < NF) lwd <- rep(lwd, length.out=NF)
    if (length(type) < NF) type <- rep(type,length.out=NF)
    if (length(lty) < NF) lty <- rep(lty,length.out=NF)
    if (length(col) < NF) col <- rep(col,length.out=NF)
    if (length(pch) < NF) pch <- rep(pch,length.out=NF)

    # }}}
    # {{{ SmartControl

    modelnames <- pframe[,unique(model)]
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(0,xlim[2],xlim[2]/4))
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(0,ylim[2],ylim[2]/4),mgp=c(4,1,0))
    if (is.character(legend[[1]])|| legend[[1]]==TRUE){
        legend.data <- getLegendData(object=x,
                                     models=models,
                                     times=tp,
                                     auc.in.legend=auc.in.legend,
                                     brier.in.legend=brier.in.legend,
                                     drop.null.model=TRUE)

        if (is.character(legend))
            legend.text <- legend
        else
            legend.text <- unlist(legend.data[,1])
        nrows.legend <- NROW(legend.data)
        if (nrows.legend==1){
            legend.lwd <- NA
        } else{ 
            legend.lwd <- lwd
        }
        legend.DefaultArgs <- list(legend=legend.text,lwd=legend.lwd,col=col,ncol=1,lty=lty,cex=cex,bty="n",y.intersp=1,x="topleft",title="")
        if (NCOL(legend.data)>1){
            addtable2plot.DefaultArgs <- list(yjust=1.18,cex=cex, table=legend.data[,-1,drop=FALSE])
        }else{
            addtable2plot.DefaultArgs <- NULL
        }
    }else{
        legend.DefaultArgs <- NULL
        addtable2plot.DefaultArgs <- NULL
    }
    if (bars){
        legend.DefaultArgs <- list(legend=modelnames,col=col,cex=cex,bty="n",x="topleft")
        names.DefaultArgs <- list(cex=.7*par()$cex,y=c(-abs(diff(ylim))/15,-abs(diff(ylim))/25))
        frequencies.DefaultArgs <- list(cex=.7*par()$cex,percent=FALSE,offset=0)
    } else{
        if (length(modelnames)<=1){
            legend=FALSE
        }
    }
    if(bars){
        legend.DefaultArgs$legend <- c("Predicted risks","Observed frequencies")
    }
    lines.DefaultArgs <- list(pch=pch,type=type,cex=cex,lwd=lwd,col=col,lty=lty)
    abline.DefaultArgs <- list(lwd=1,col="red")
    if (missing(ylim)){
        if (showPseudo[1] && !bars[1]){
            ylim <- range(pframe[[Rvar]])
        }
        else
            ylim <- c(0,1)
    }
    if (missing(xlim)){
        xlim <- c(0,1)
    }
    plot.DefaultArgs <- list(x=0,
                             y=0,
                             type = "n",
                             ylim = ylim,
                             xlim = xlim,
                             ylab=ylab,
                             xlab=xlab)
    rug.DefaultArgs <- list(col=col,lwd=lwd,side=1,ticksize=0.03)
    pseudo.DefaultArgs <- list(col=col,cex=cex,density=55)
    barplot.DefaultArgs <- list(ylim = ylim,
                                col=col,
                                axes=FALSE,
                                ylab=ylab,
                                xlab=xlab,
                                beside=TRUE,
                                legend.text=NULL,
                                cex.axis=cex,
                                cex.lab=par()$cex.lab,
                                cex.names=cex)
    if (bars){
        control <- prodlim::SmartControl(call= list(...),
                                         keys=c("barplot","legend","addtable2plot","axis2","abline","names","frequencies"),
                                         ignore=NULL,
                                         ignore.case=TRUE,
                                         defaults=list("barplot"=barplot.DefaultArgs,
                                                       "addtable2plot"=addtable2plot.DefaultArgs,
                                                       "abline"=abline.DefaultArgs,
                                                       "legend"=legend.DefaultArgs,
                                                       "names"=names.DefaultArgs,
                                                       "frequencies"=frequencies.DefaultArgs,
                                                       "axis2"=axis2.DefaultArgs),
                                         forced=list("abline"=list(h=0)),
                                         verbose=TRUE)
    }else{
        control <- prodlim::SmartControl(call= list(...),
                                         keys=c("plot","rug","pseudo","lines","legend","addtable2plot","axis1","axis2"),
                                         ignore=NULL,
                                         ignore.case=TRUE,
                                         defaults=list("plot"=plot.DefaultArgs,"pseudo"=pseudo.DefaultArgs,"rug"=rug.DefaultArgs,
                                                       "addtable2plot"=addtable2plot.DefaultArgs,
                                                       "lines"=lines.DefaultArgs,
                                                       "legend"=legend.DefaultArgs,
                                                       "axis1"=axis1.DefaultArgs,
                                                       "axis2"=axis2.DefaultArgs),
                                         forced=list("plot"=list(axes=FALSE),
                                                     "axis1"=list(side=1)),
                                         verbose=TRUE)
    }

    # }}}
    # {{{ smoothing
    method <- match.arg(method,c("quantile","nne"))
    getXY <- function(f){
        risk=NULL
        p <- pframe[model==f,risk]
        jackF <- pframe[model==f][[Rvar]]
        switch(method,
               "quantile"={
                   if (length(q)==1)
                       groups <- quantile(p,seq(0,1,1/q))
                   else{
                       groups <- q
                   }
                   xgroups <- (groups[-(length(groups))]+groups[-1])/2
                   pcut <- cut(p,groups,include.lowest=TRUE)
                   ## if (x$cens.type=="rightCensored"){
                   plotFrame=data.frame(Pred=tapply(p,pcut,mean),Obs=pmin(1,pmax(0,tapply(jackF,pcut,mean))))
                   attr(plotFrame,"quantiles") <- groups
                   plotFrame
               },
               "nne"={
                   ## Round probabilities to 2 digits
                   ## to avoid memory explosion ...
                   ## a difference in the 3 digit should
                   ## not play a role for the single subject.
                   if (round==TRUE){
                       if (!is.null(bandwidth) && bandwidth[1]>=1){
                           ## message("No need to round predicted probabilities to calculate calibration in the large")
                       } else{
                           p <- round(p,2)
                       }
                   }
                   p <- na.omit(p)
                   if (no <- length(attr(p,"na.action")))
                       warning("calPlot: removed ",no," missing values in risk prediction.",call.=FALSE,immediate.=TRUE)
                       if (is.null(bandwidth)){
                           ## if (length(p)>length(apppred[,f+1])){
                           ## bw <- prodlim::neighborhood(apppred[,f+1])$bandwidth
                           ## }else{
                           if (length(unique(p))==1)
                               bw <- 1
                           else 
                               bw <- prodlim::neighborhood(p)$bandwidth
                           ## }
                       } else{
                           bw <- bandwidth
                       }                       
                       if (bw>=1){
                           ## calibration in the large
                           plotFrame <- data.frame(Pred=mean(p),Obs=mean(jackF))
                       } else{
                           ## nn <- prodlim::neighborhood(x=p,bandwidth=bw)
                           nbh <- prodlim::meanNeighbors(x=p,y=jackF,bandwidth=bw)
                           plotFrame <- data.frame(Pred=nbh$uniqueX,Obs=nbh$averageY)
                       }
                       attr(plotFrame,"bandwidth") <- bw
                       plotFrame
               })
    }
    plotFrames <- lapply(modelnames,function(f){getXY(f)})
    names(plotFrames) <- modelnames
    # }}}
    # {{{ plot and/or invisibly output the results

    if (bars){
        if ((is.logical(names[[1]]) && names[[1]]==TRUE)|| names[[1]] %in% c("quantiles.labels","quantiles")){
            qq <- attr(plotFrames[[1]],"quantiles")
            if (names[1]=="quantiles.labels"){
                pp <- seq(0,1,1/q)
                names <- paste0("(",
                                sprintf("%1.0f",100*pp[-length(pp)]),",",
                                sprintf("%1.0f",100*pp[-1]),
                                ")\n",
                                sprintf("%1.1f",100*qq[-length(qq)])," - ",
                                sprintf("%1.1f",100*qq[-1]))
            }
            else 
                names <- paste0(sprintf("%1.1f",100*qq[-length(qq)])," - ",
                                sprintf("%1.1f",100*qq[-1]))
        }
    }
    out <- list(plotFrames=plotFrames,
                times=tp,
                cause=cause,
                control=control,
                legend=legend,
                bars=bars,
                diag=diag,
                add=add,
                names=names,
                method=method,
                axes=axes,
                percent=percent,
                hanging=hanging,
                show.frequencies=show.frequencies,
                col=col,
                ylim=ylim,
                xlim=xlim,
                ylab=ylab,
                xlab=xlab,
                lwd=lwd,
                lty=lty,
                pch=pch,
                lty=lty,
                type=type,
                NF=NF)
                
    if (method=="nne")
        out <- c(out,list(bandwidth=sapply(plotFrames,
                                           function(x)attr(x,"bandwidth"))))
    # }}}
    # {{{ do the actual plot

    if (plot){
        if (out$add[1]==FALSE && !out$bars[1]){
            do.call("plot",control$plot)
        }
        if (out$diag[1] && !out$bars[1]){
            segments(x0=0,y0=0,x1=1,y1=1,col="gray77",lwd=2,xpd=FALSE)
        }
        if (out$diag[1] && !out$bars[1]){
            segments(x0=0,y0=0,x1=1,y1=1,col="gray77",lwd=2,xpd=FALSE)
        }
        if (out$bars[1]) {
            stopifnot(NF==1)
            pf <- na.omit(plotFrames[[1]])
            Pred <- pf$Pred
            Obs <- pf$Obs
            if (is.character(legend[[1]]) || legend[[1]]==FALSE){
                control$barplot$legend.text <- NULL
            }else{
                if (is.null(control$barplot$legend.text)){
                    control$barplot$legend.text <- control$legend$legend
                }
                ## }else{
                control$barplot$args.legend <- control$legend
                ## }
            }
            if (is.null(control$barplot$space))
                control$barplot$space <- rep(c(1,0),length(Pred))
            PredObs <- c(rbind(Pred,Obs))
            control$barplot$height <- PredObs
            if (out$hanging){
                control$barplot$offset <- c(rbind(0,Pred-Obs))
                minval <- min(Pred-Obs)
                if (minval<0)
                    negY.offset <- 0.05+seq(0,1,0.05)[prodlim::sindex(jump.times=seq(0,1,0.05),eval.times=abs(minval))]
                else
                    negY.offset <- 0
                control$barplot$ylim[1] <- min(control$barplot$ylim[1],-negY.offset)
                control$names$y <- control$names$y-negY.offset
            }
            coord <- do.call("barplot",control$barplot)
            if (length(out$names)>0 && (out$names[[1]]!=FALSE) && is.character(out$names)){
                if (out$names[[1]]!=FALSE && length(out$names)==(length(coord)/2)){
                    mids <- rowMeans(matrix(coord,ncol=2,byrow=TRUE))
                    text(x=mids,
                         ## x=coord,
                         y=control$names$y,
                         ## c(rbind(out$names,rbind(rep("",length(coord)/2)))),
                         out$names,
                         xpd=NA,
                         cex=control$names$cex)
                }
            }
            ## if (out$legend) print(control$barplot$args.legend)n
            ## message(paste0("Bars are located at ",paste(coord,collapse=",")))
            if (out$hanging){
                do.call("abline",control$abline)
            }
            if (out$show.frequencies){
                if(out$hanging){
                    text(x=coord,
                         cex=control$frequencies$cex,
                         pos=3,
                         y=(as.vector(rbind(Pred,Pred)) +rep(control$frequencies$offset,times=length(as.vector(coord))/2)),
                         paste(round(100*c(rbind(Pred,Obs)),0),ifelse(control$frequencies$percent,"%",""),sep=""),xpd=NA)
                }else{
                    text(coord,
                         pos=3,
                         c(rbind(Pred,Obs))+control$frequencies$offset,
                         cex=control$frequencies$cex,
                         paste(round(100*c(rbind(Pred,Obs)),0),ifelse(control$frequencies$percent,"%",""),sep=""),xpd=NA)
                }
            }
            coords <- list(xcoord=coord[,1],ycoord=PredObs,offset=control$barplot$offset)
            out <- c(out,coords)
        }else{
            nix <- lapply(1:NF,function(f){
                pf <- out$plotFrames[[f]]
                pf <- na.omit(pf)
                ## rug plot predictions
                if (showRug){
                    do.call(graphics::rug,c(list(x=pf$Pred,col=control$rug$col[f],lwd=control$rug$lwd[f]*0.5)))
                }
                ## show pseudo values 
                if (showPseudo) {
                    if (!is.null(control$pseudo$density) & control$pseudo$density>0){
                        control$pseudo$col <- prodlim::dimColor(control$pseudo$col[f],
                                                                control$pseudo$density)
                    }
                    if ((gotcha <- match("density",names(control$pseudo),nomatch=0))>0){
                        control$pseudo <- control$pseudo[-gotcha]
                    }
                    do.call(points, c(list(x=pframe[model==modelnames[f]][,risk],y=pframe[model==modelnames[f]][[Rvar]]),control$pseudo))
                }
                ## add the lines
                lineF <- lapply(control$lines,function(x)x[[min(f,length(x))]])
                lineF$x <- pf$Pred
                lineF$y <- pf$Obs
                do.call("lines",lineF)
            })
            if (is.character(out$legend[[1]]) || out$legend[[1]]==TRUE){
                legend.coords <- do.call("legend",control$legend)
                if (!is.null(addtable2plot.DefaultArgs)){
                    if (is.null(control$addtable2plot[["x"]]))
                        control$addtable2plot[["x"]] <- legend.coords$rect$left+legend.coords$rect$w
                    ## strange error: $y does not work as yjust is matched, thus use [["y"]]
                    if (is.null(control$addtable2plot[["y"]]))
                        control$addtable2plot[["y"]] <- legend.coords$rect$top-legend.coords$rect$h
                    do.call(plotrix::addtable2plot,control$addtable2plot)
                }
            }
        }
        if (out$axes){
            if (out$percent){
                if (is.null(control$axis1$labels)){
                    control$axis1$labels <- paste(100*control$axis1$at,"%")
                }
                if (is.null(control$axis2$labels)){
                    control$axis2$labels <- paste(100*control$axis2$at,"%")
                }
            }
            if (!out$bars)
                do.call("axis",control$axis1)
            do.call("axis",control$axis2)
        }
    }

    # }}}
    class(out) <- "calibrationPlot"
    invisible(out)
}

#----------------------------------------------------------------------
### plotCalibration.R ends here
