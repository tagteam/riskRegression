#########################################
# Function 'predictEventProb.selectFGR' #
#########################################

#Authors:  Thomas A. Gerds & Rob C.M. van Kruijsdijk
#Date: 20-02-2013

##Arguments:
# FGR.obj :  A fitted FGR-model from which to extract predicted event probabilities, including all candidate covariates (cencode must be 0)
# event:     Integer variable that denotes type of failure for each person, also see crrstep
# data:      A data frame containing predictor variable combinations for which to compute predicted event probabilities.
# rule:      Rule to pass on to crrstep ("AIC", "BIC", "BICcr"), also see crrstep
# direction: Direction of model selection ("backward" or "forward")

selectFGR <- function(FGR.obj,event,data,rule="AIC",direction="backward"){
  
  require(crrstep)
  require(prodlim)
  crrstep.form <- reformulate(paste(FGR.obj$call[[2]][3]),response=all.vars(as.formula(FGR.obj$call[2]))[1])
  
  crrstep.fit <- do.call("crrstep",list(formula=crrstep.form,data=substitute(data),etype=substitute(event), 
                                       failcode=FGR.obj$cause, cencode=0, direction = direction, criterion = rule, 
                                       crr.object = FALSE, trace = FALSE))
  
  #The results of crrstep are printed after this, which might not be desirable in a CV-loop, but I do not know how to prevent this, without changing the crrstep-code..
  #..Setting crr.object to TRUE works, but then the length of $coef is 1 when all variables are dropped, which complicates the next step
  
  if (length(crrstep.fit$coefficients)==0){
    newform <- reformulate("1",FGR.obj$call[[2]][[2]])
    newfit <- prodlim(newform,data=data) #Is this correct for competing risks (does prodlim apply CIF or KM in this case)??
  }
  else{
    newform <- reformulate(paste(rownames(crrstep.fit$coefficients),collapse="+"),response=FGR.obj$call[[2]][[2]])
    newfit <- FGR(newform,data=data,cause=FGR.obj$cause)
    newfit$call$formula <- newform
  }
  out <- list(fit=newfit,In=rownames(crrstep.fit$coefficients))
  out$call <- match.call()
  class(out) <- "selectFGR"
  out
}

predictEventProb.selectFGR <- function(object,newdata,times,...){
  predictEventProb(object[[1]],newdata=newdata,times=times,...)
}

###########
# TESTING #
###########
## library(riskRegression)
## library(cmprsk)
#Simulation of data (adjusted from FGR documentation):
## d <- prodlim:::SimCompRisk(100)
## newvars <- as.data.frame(matrix(runif(8*100),nrow=100))
## colnames(newvars) <- c('X3','X4','X5','X6','X7','X8','X9','X10')
## d <- cbind(d,newvars)
## evaltime <- quantile(d$time)[3]
## mod.1<- FGR(formula = Hist(time, cause) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, cause=1, data = d)
## test.1 <- selectFGR(FGR.obj=mod.1,data=d,event=cause) #Seems to work fine..
## class(test.1[[1]]) #class "FGR"
## predictEventProb.selectFGR(test.1,newdata=d,times=evaltime) #I get: Error in ff[[2]] : subscript out of bounds
## cindex(list(test.1),formula=Hist(time,cause)~1,data=d)

