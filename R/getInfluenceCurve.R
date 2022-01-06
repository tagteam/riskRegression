## IMPORTANT : data have to be ordered by time (and for ties by reverse status,
## such that events (status=1,2,...) come before censored (status=0))
getInfluenceCurve.AUC.survival <- function(t,n,time,risk,Cases,Controls,ipcwControls,ipcwCases,MC){
    F01t <- sum(ipcwCases)
    St <- sum(ipcwControls)
    nbCases <- sum(Cases)
    nbControls <- sum(Controls)
    if (nbCases==0 || nbControls ==0 || length(unique(risk))==1) return(rep(0,n))
    # mcase <- matrix(risk[Cases],nrow=nbCases,ncol=nbControls)
    # mcontrol <- matrix(risk[Controls],nrow=nbCases,ncol=nbControls,byrow=TRUE)
    # wcase <- matrix(ipcwCases[Cases],nrow=nbCases,ncol=nbControls)
    # wcontrol <- matrix(ipcwControls[Controls],nrow=nbCases,ncol=nbControls,byrow=TRUE)
    #htij1 <- (1*(mcase>mcontrol)+.5*(mcase==mcontrol))*wcase*wcontrol*n*n
    htij1 <- htijCalculationHelper(risk[Cases],risk[Controls],ipcwCases[Cases],ipcwControls[Controls],n,nbCases,nbControls)
    ht <- (sum(htij1))/(n*n)
    fi1t <- Cases*ipcwCases*n
    colSumshtij1 <- rep(0,n) # initialise at 0
    colSumshtij1[Cases] <- rowSums(htij1) 
    rowSumshtij1 <- rep(0,n) # initialize at 0
    rowSumshtij1[Controls] <- colSums(htij1)
    hathtstar <- (sum(htij1))/(n*n)
    vectTisupt <- n*Controls/sum(Controls)
    # Integral_0^T_i dMC_k/S for i %in% Cases
    MC.Ti.cases <- MC[sindex(eval.times=time[Cases],jump.times=unique(time)),,drop=FALSE]
    T1 <- rowSumsCrossprod(htij1,1+MC.Ti.cases,0)/n
    ## print(system.time(T1 <- colSums(crossprod(htij1,1+MC.Ti.cases))/n))
    ## the matrix MC has as many rows as there are unique time points
    ## and as many colums as there are subjects
    ## in case of ties need to expand MC to match the dimension of fi1t
    #MC.all <- MC
    if (all(match(time,unique(time)) != 1:ncol(MC))) {
        MC <- MC[match(time,unique(time)),]
    }
    #T3 <- colSums(hathtstar*(vectTisupt + (fi1t*(1+MC)-F01t)/F01t))
    #T3 <- hathtstar*sum(vectTisupt) + hathtstar/F01t * (T3CalculationHelper(fi1t,MC)-nrow(MC)*F01t)
    T3 <- hathtstar*(sum(vectTisupt) + 1/F01t * T3CalculationHelper(fi1t,MC)-nrow(MC))
    
    Term.ijak <- (T1-T3)/(F01t*St)
    Term.ikaj <- (rowSumshtij1 - n*hathtstar)/(F01t*St)
    Term.jkai <-  (colSumshtij1 - n*hathtstar*(vectTisupt+(1/F01t)*(fi1t-F01t)))/(F01t*St)
    ## the influence function according to Blanche et al. 2013, DOI: 10.1002/sim.5958, Statistics in Medicine, Appendix A
    (Term.ijak + Term.ikaj + Term.jkai)/(n)
}


## NTC <- NCOL(MC.Ti.cases)
## T1 <- numeric(NTC)
## for (j in 1:NTC){
## T1[j] <- sum(cprod(htij1,(1+MC.Ti.cases[,j,drop=FALSE]),0,1,0)/n)
## }
## system.time(T1a <- sapply(1:NCOL(MC.Ti.cases),function(j){sum(cprod(htij1,(1+MC.Ti.cases[,j,drop=FALSE]),0,1,0)/n)}))
## system.time(T1 <- colSums(crossprod(htij1,1+MC.Ti.cases))/n)
## system.time(T1 <- (crossprod(crossprod(one,t(htij1)),1+MC.Ti.cases))/n)
## A <- big.matrix(NROW(htij1),NCOL(htij1))
## A[,] <- htij1
## B <- big.matrix(NROW(MC.Ti.cases),NCOL(MC.Ti.cases))
## B[,] <- 1+MC.Ti.cases
## T1 <- t(A)%*%B
## T1 <- cprod(htij1,1+MC.Ti.cases,1,1,0)/n
## system.time(T1 <- cprod(htij1,1+MC.Ti.cases,1,1,0)/n)
## system.time(T1b <- cprod(A,1+MC.Ti.cases,1,0,0)/n)


## IMPORTANT : data have to be ordered by time (and for ties by reverse status)
getInfluenceCurve.AUC.competing.risks <- function(t,n,time,risk,Cases,Controls1,Controls2,ipcwControls1,ipcwControls2,ipcwCases,MC){
    F01t <- sum(ipcwCases)
    ## 1-F01t <- sum(ipcwControls1+ipcwControls2)
    nbCases <- sum(Cases)
    nbControls1 <- sum(Controls1)
    nbControls2 <- sum(Controls2)
    if (nbCases==0 || (nbControls1 + nbControls2) ==0  || length(unique(risk))==1) return(rep(0,n))
    # #a(i,j)=a(i,j')
    # mcase1 <- matrix(risk[Cases],nrow=nbCases,ncol=nbControls1)
    # mcase2 <- matrix(risk[Cases],nrow=nbCases,ncol=nbControls2)
    # #a(i,j) = a(i',j)
    # mcontrol1 <- matrix(risk[Controls1],nrow=nbCases,ncol=nbControls1,byrow=TRUE)
    # mcontrol2 <- matrix(risk[Controls2],nrow=nbCases,ncol=nbControls2,byrow=TRUE)
    # wcase1 <- matrix(ipcwCases[Cases],nrow=nbCases,ncol=nbControls1)
    # wcase2 <- matrix(ipcwCases[Cases],nrow=nbCases,ncol=nbControls2)
    # wcontrol1 <- matrix(ipcwControls1[Controls1],nrow=nbCases,ncol=nbControls1,byrow=TRUE)
    # wcontrol2 <- matrix(ipcwControls2[Controls2],nrow=nbCases,ncol=nbControls2,byrow=TRUE)
    #htij1 <- (1*(mcase1>mcontrol1)+.5*(mcase1==mcontrol1))*wcase1*wcontrol1*n*n
    #htij1<-htijCalculationHelper(mcase1,mcontrol1,wcase1,wcontrol1,n)
    htij1<-htijCalculationHelper(risk[Cases],risk[Controls1],ipcwCases[Cases],ipcwControls1[Controls1],n,nbCases,nbControls1)
    #htij2 <- (1*(mcase2>mcontrol2) + .5*(mcase2==mcontrol2))*wcase2*wcontrol2*n*n
    htij2<-htijCalculationHelper(risk[Cases],risk[Controls2],ipcwCases[Cases],ipcwControls2[Controls2],n,nbCases,nbControls2)
    #htij2<-htijCalculationHelper(mcase2,mcontrol2,wcase2,wcontrol2,n)
    ht <- (sum(htij1)+sum(htij2))/(n*n)
    fi1t <- Cases*ipcwCases*n
    colSumshtij1 <- rep(0,n) # initialise at 0
    colSumshtij1[Cases] <- rowSums(htij1)
    colSumshtij2 <- rep(0,n) # initialise at 0
    colSumshtij2[Cases] <- rowSums(htij2) 
    rowSumshtij1 <- rep(0,n) # match(time,unique(time))initialize at 0
    rowSumshtij1[Controls1] <- colSums(htij1)
    rowSumshtij2 <- rep(0,n) # initialize at 0
    rowSumshtij2[Controls2] <- colSums(htij2)
    hathtstar <- (sum(htij1))/(n*n)  
    vectTisupt <- n*Controls1/nbControls1
    # Integral_0^T_i dMC_k/S for i %in% Cases
    MC.Ti.cases <- MC[sindex(eval.times=time[Cases],jump.times=unique(time),),,drop=FALSE]
    ## MC.Ti.cases <- rbind(0,MC)[1+prodlim::sindex(eval.times=time[Cases],jump.times=unique(time),strict=1),,drop=FALSE]
    # Integral_0^T_i dMC_k/S for i %in% Controls 2
    MC.Ti.controls2 <- MC[sindex(eval.times=time[Controls2],jump.times=unique(time)),,drop=FALSE]
    ## MC.Ti.controls2 <- rbind(0,MC)[1+prodlim::sindex(eval.times=time[Controls2],jump.times=unique(time),strict=1),,drop=FALSE]
    # Integral_0^t dMC_k/S for all i
    MC.t <- MC[prodlim::sindex(eval.times=t,jump.times=unique(time)),,drop=TRUE]
    ## MC.t <- rbind(0,MC)[1+prodlim::sindex(eval.times=t,jump.times=unique(time),strict=1),,drop=TRUE]
    # we compute \frac{1}{n}\sum{i=1}^n \sum{j=1}^n \sum{l=1}^n \Psi{ijkl}(t)
    T1 <- rowSumsCrossprod(htij1,1+MC.Ti.cases,0)
    ## T1 <- colSums(crossprod(htij1,1+MC.Ti.cases))
    T2 <- rowSumsCrossprod(htij2,1+MC.Ti.cases,0)
    ## T2 <- colSums(crossprod(htij2,1+MC.Ti.cases))
    ## in case of ties need to expand MC to match the dimension of fi1t
    # data does not have ties, this saves memory 
    if (all(match(time,unique(time)) != 1:ncol(MC))) {
        MC <- MC[match(time,unique(time)),]
    }
    #MC.all <- MC[match(time,unique(time)),]
    #should be quicker
    #T3 <- ht*(1-2*F01t)/(F01t*(1-F01t))*colSums(fi1t*(1+MC))-ht*(1-2*F01t)/(F01t*(1-F01t))*nrow(MC)*F01t
    T3 <- ht*(1-2*F01t)/(F01t*(1-F01t))*T3CalculationHelper(fi1t,MC)-ht*(1-2*F01t)/(F01t*(1-F01t))*nrow(MC)*F01t
    #T3 <- colSums((ht*(1-2*F01t)/(F01t*(1-F01t)))*(fi1t*(1+MC.all)-F01t))
    Term.ijlk <- ((T1 + T2) - n^2*ht - n*T3)/(F01t*(1-F01t))
    # we compute \frac{1}{n}\sum_{i=1}^n \sum_{j=1}^n \sum_{k=1}^n \Psi_{ijkl}(t)
    # Q1 <- sapply(1:n,function(i)sum(htij1*(1+MC.t[i])))
    Q1 <- sum(htij1)*(1+MC.t)
    ## Q2 <- colSums(crossprod(t(htij2),(1+MC.Ti.controls2)))
    Q2 <- colSumsCrossprod(htij2,(1+MC.Ti.controls2),0)
    Term.ijkl <- ((Q1 + Q2) - n^2*ht)/(F01t*(1-F01t))
    # we compute \frac{1}{n}\sum_{j=1}^n \sum_{k=1}^n \sum_{l=1}^n \Psi_{ijkl}(t)
    Term.jkli <-((colSumshtij1 + colSumshtij2)*n - n^2*ht - ( ht*n^2*(1-2*F01t) / (F01t*(1-F01t)) ) *(fi1t - F01t) )/(F01t*(1-F01t))
    # we compute \frac{1}{n}\sum{i=1}^n \sum{k=1}^n \sum{l=1}^n \Psi{ijkl}(t)
    Term.iklj<-((rowSumshtij1 + rowSumshtij2)*n - n^2*ht)/(F01t*(1-F01t))
    ## the influence function according to Blanche et al. 2013, DOI: 10.1002/sim.5958, Statistics in Medicine, Appendix A
    ## print(list(Term.jkli,Term.iklj,Term.ijkl,Term.ijlk))
    ## print(rbind(Term.jkli, Term.iklj, Term.ijkl, Term.ijlk))
    ## print((Term.jkli + Term.iklj + Term.ijkl + Term.ijlk)/(n*n))
    (Term.jkli + Term.iklj + Term.ijkl + Term.ijlk)/(n*n)
}

getInfluenceCurve.Brier1 <- function(t,time,Yt,residuals,MC){
    hit1=(Yt==0)*residuals
    hit2=(Yt==1)*residuals
    Brier <- mean(residuals)
    ## FIXME: make sure that sindex cannot be 0
    ## browser(skipCalls=1)
    ## if (length(dim(MC))==3) browser()
    Int0tdMCsurEffARisk <- MC[prodlim::sindex(jump.times=unique(time),eval.times=t),,drop=FALSE]
    IF.Brier=hit1+hit2-Brier + mean(hit1)*Int0tdMCsurEffARisk + colMeans(MC*hit2)
}

getInfluenceCurve.Brier <- function(t,
                                    time,
                                    IC0,
                                    residuals,
                                    WTi,
                                    Wt,
                                    IC.G,
                                    cens.model,
                                    nth.times=NULL){
    ##
    ## Compute influence function of Brier score estimator using weights of the reverse Cox model 
    ## This function evaluates the part of influence function which is related to the IPCW weights
    ## The other part is IC0.
    ## 
    ## \frac{1}{n}\sum_{i=1}^n 
    ## m_{t,n}^{(1)}(X_i) 
    ## [\frac{I_{T_i\leq t}\Delta_i}{G^2(T_i\vert Z_i)}IF_G(T_i,X_k; X_i)+\frac{I_{T_i>t}}{G^2(t|Z_i)}IF_G(t,X_k; X_i)] 
    ## with 
    ## IF_G(t,X_k; X_i)=-\exp(-\Lambda(t\vert Z_i))IF_{\Lambda}(t,X_k; X_i) 
    ##
    ## IC_G(t,z;x_k) is an array with dimension (nlearn=N, gtimes, newdata)
    ## where gtimes = subject.times (Weights$IC$IC.subject) or times (Weights$IC$IC.times)
    ## and subject.times=Y[(((Y<=max(times))*status)==1)]
    ##
    ## don't square the weights because they will be multiplied with the
    ## residuals that are already weighted 
    ## 
    N <- length(residuals)
    if (cens.model=="cox") {## Cox regression
        ic.weights <- matrix(0,N,N)
        k=0 ## counts subject-times with event before t
        for (i in 1:N){
            if (residuals[i]>0){
                if (time[i]<=t){ ## min(T,C)<=t, note that (residuals==0) => (status==0)  
                    k=k+1
                    ic.weights[i,] <- IC.G$IC.subject[i,k,]/(WTi[i])
                }else{## min(T,C)>t
                    ic.weights[i,] <- IC.G$IC.times[i,nth.times,]/(Wt[i])
                }
            }
        }
        IF.Brier <- ic.weights*residuals
        IF.Brier <- IC0-colMeans(IF.Brier)
        IF.Brier
    }else{
        ## Blanche et al. 2015 (joint models) web appendix equation (14)
        Yt <- 1*(time<=t)
        hit1=(time>t)*residuals ## equation (7)
        hit2=(time<=t)*residuals ## equation (8) 
        Brier <- mean(residuals)
        if (!is.null(IC.G)){
            if (prodlim::sindex(jump.times=unique(time),eval.times=t) > 1){
                #Int0tdMCsurEffARisk <- IC.G[prodlim::sindex(jump.times=unique(time),eval.times=t),,drop=FALSE]
                #IF.Brier <- hit1+hit2-Brier + mean(hit1)*Int0tdMCsurEffARisk+ colMeans(IC.G*hit2)
                IF.Brier <- hit1+hit2-Brier + mean(hit1)*IC.G[prodlim::sindex(jump.times=unique(time),eval.times=t),,drop=FALSE]+ columnMeanWeight(IC.G,hit2)
            }
            else {
                #Int0tdMCsurEffARisk <- rep(0,ncol(IC.G))
                IF.Brier <- hit1+hit2-Brier + colMeans(IC.G*hit2)
            }
            #Int0tdMCsurEffARisk <- rbind(0,IC.G)[1+prodlim::sindex(jump.times=unique(time),eval.times=t),,drop=FALSE]
            
        }else{# uncensored
            IF.Brier <- hit1+hit2-Brier
        }
        IF.Brier
    }
}


## now using cpp function ../src/IC_Nelson_Aalen_cens_time.cpp
getInfluenceCurve.NelsonAalen <- function(time,status){
    ##
    ## compute influence function for reverse Nelson-Aalen estimator
    ## to deal with ties we sort such that events come before censored
    ## but we do not collapse times at ties so that the at-risk set is
    ## slowly decreased "during" a tie and we sort such that events come
    ## before censored such that, "during" a tie, the censoring hazard first
    ## changes at the row (time in rows, subjects in columns) of the first
    ## censored subject. Also, the indicator min(T_i,C_i)<=t is set to zero
    ## if min(T_i,C_i)==t but t<i "during" a tie.
    ##
    ## i = 1,..., n are columns
    ## s = 1,..., s_tmax are rows
    ## ----- no need to order when called from Score
    ## neworder <- order(time,-status)
    ## time <- time[neworder]
    ## status <- status[neworder]
    ## ----- no need to order when called from Score
    n <- length(time)
    # compute hazard function of the censoring
    hazardC <- (status==0)/(n:1)
    # probability to be at risk
    atrisk <- (n:1) # atrisk/n = prob(min(T,C)>=t)
    # matrix with one column for each subject and one row for each time point
    hatMC <- do.call("cbind",lapply(1:n,function(i){
        (status[i]==0)*rep(c(0,1),c(i-1,n-i+1))*n/atrisk[i]-cumsum(c(hazardC[1:i], rep(0,(n-i)))*n/atrisk)
    }))
    hatMC
}

getInfluenceCurve.NelsonAalen.slow <- function(time,status){    
    time <- time[order(time)]
    status <- status[order(time)] 
    n <- length(time)
    mat.data<-cbind(time,as.numeric(status==0))
    colnames(mat.data)<-c("T","indic.Cens")
    # compute the empirical survival function corresponding to the counting process 1(\tilde{eta}=0, \tilde{T}<=t)
    hatSdeltaCensTc<-1-cumsum(mat.data[,c("indic.Cens")])/n  
    # Build the matrix required for computing  dM_C(u) for all time u (all observed times \tilde{T}_i)
    temp1 <- cbind(mat.data[,c("T","indic.Cens")],1-(1:n)/n,hatSdeltaCensTc)
    temp1 <- rbind(c(0,0,1,1),temp1) # Add the first row corresponding to time t=0
    colnames(temp1)<-c("T","indic.Cens","hatSTc","hatSdeltaCensTc")
    # compute hazard function of the censoring
    lambdaC<-(temp1[-1,"indic.Cens"])/(n:1)
    # Add the column of the hazard function of the censoring (equal to 0 at time t=0)
    temp1<-cbind(temp1,c(0,lambdaC))
    colnames(temp1)[ncol(temp1)]<-"lambdaC"
    # Cumulative hazard of censoring
    LambdaC<-cumsum(lambdaC)         
    # Add the column of the cumulative hazard function of the censoring (equal to 0 at time t=0)
    temp1 <- cbind(temp1,c(0,LambdaC))
    colnames(temp1)[ncol(temp1)]<-"LambdaC"
    temp2<-temp1[-1,]
    # compute  martingale of censoring \hat{M}_{C_i}(u) for all time u (all observed times \tilde{T}_i) using previous matrix
    # We obtain a matrix. Each column contains the vector of M_{C_i}(\tilde{T}_j) for  all j.
    hatMC<-matrix(NA,n,n)
    for (i in 1:n){
        hatMC[,i] <-temp2[i,2]*as.numeric(temp2[i,1]<=temp2[,"T"])- c(temp2[0:i,"LambdaC"], rep(temp2[i,6],(n-i)))
    }  
    # In order to draw martingale paths
    #matplot(mat.data[,"T"],hatMC,type="l")
    #lines(mat.data[,"T"],rowMeans(hatMC),lwd=5)  
    # Compute d \hat{M}_{C_i} (u) for all time u (all observed times \tilde{T}_i)
    dhatMC<-rbind(hatMC[1,],hatMC[-1,]-hatMC[-nrow(hatMC),])
    # Compute d \hat{M}_{C_i} (u)/(S_{\tilde{T}}(u)) for all time u (all observed times \tilde{T}_i)
    # We need this for integrals in the martingale representation of the Kaplan-Meier estimator of the censoring survival function
    # function to divide d \hat{M}_{C_i} (u) by (S_{\tilde{T}}(u))
    MulhatSTc<-function(v){
        n <- length(v)
        v/c(1,1-(1:(n-1))/n)      # c(1,1-(1:(n-1))/n) is the at risk probability (S_{\tilde{T}}(u))
    }
    # apply the function for each column (corresponding to the
    # vector M_{C_i}(u)  for all time u (all observed times \tilde{T}_i), 
    # time \tilde{T}_i corresponds to the i-th row of the matrix)
    dhatMCdivST<-apply(dhatMC,2,MulhatSTc)
    # Compute \int_0^{\tilde{T}_j} d{ \hat{M}_{C_l} (u) } / (S_{\tilde{T}}(u)) for each subject l, we compute for all time \tilde{T}_j.
    # l=column, j=row
    MatInt0TcidhatMCksurEff<-apply(dhatMCdivST,2,cumsum)  # (Remark : on of the row corresponds to the previous step...) 
    ## colnames(MatInt0TcidhatMCksurEff)<-paste("M_{C_",1:length(time),"}",sep="")
    ## rownames(MatInt0TcidhatMCksurEff)<-time  
    return(MatInt0TcidhatMCksurEff)  
}


getInfluenceCurve.KM <- function(time,status){
    ## compute influence function for reverse Nelson-Aalen
    ## i = 1,..., n are columns
    ## s = 1,..., s_tmax are rows
    N <- length(time)
    times <- unique(time)
    NU <- length(times)
    lagtime <- c(0,times[-NU])
    dd <- data.frame(time=time,status=status)
    F <- prodlim::prodlim(Hist(time,status)~1,data=dd,reverse=FALSE)
    G <- prodlim::prodlim(Hist(time,status)~1,data=dd,reverse=TRUE)
    Stilde.T <- prodlim::predictSurvIndividual(F,lag=1)*prodlim::predictSurvIndividual(G,lag=1)
    Stilde.s <- predict(F,times=lagtime)*predict(G,times=lagtime)
    out <- lapply(1:N,function(i){
        ((1-status[i])*(time[i]<=times))/Stilde.T[i] - cumsum((time[i]>=times)*(G$hazard)/Stilde.s)
    })
    do.call("cbind",out)
}


