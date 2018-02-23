### getVcov.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Feb 23 2018 (14:03) 
## Version: 
## Last-Updated: Feb 23 2018 (17:08) 
##           By: Thomas Alexander Gerds
##     Update #: 19
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getVcov <- function(data,IF.name,times=NULL){
    model=models=n=times=NULL
    models <- data[,unique(model)]
    if (!is.null(times)){
        times <- data[,unique(times)]
        N <- data[model==model[[1]] & times==times[[1]],.N]
        AllComb <- expand.grid(model=models,times=times)
        allnames <- apply(AllComb,1,function(x){paste0("model=",x[1],"times=",x[2],sep="")})
        matVoCov <- matrix(0,length(allnames),length(allnames))
        rownames(matVoCov) <- allnames
        colnames(matVoCov) <- allnames    
        for(i in 1:nrow(AllComb)){    # First pair of model and time
            themodel1 <- AllComb[i,"model"]
            thetimes1 <- AllComb[i,"times"]
            ## print(paste("Compute for model 1 =",themodel1,"and times 1=",thetimes1))
            for(j in i:nrow(AllComb)){ # Second pair of model and time
                themodel2 <- AllComb[j,"model"]
                thetimes2 <- AllComb[j,"times"]
                ## print(paste("Compute for model 2 =",themodel2,"and times 2=",thetimes2))
                # extract the two iid decompositions
                IF1 <- data[model==themodel1 & times==thetimes1][[IF.name]]
                IF2 <- data[model==themodel2 & times==thetimes2][[IF.name]]
                # compute and save covariance
                matVoCov[i,j] <- cov(IF1,IF2)/N       
            }
        }
    } else{
        N <- data[model==model[[1]],.N]
        AllComb <- expand.grid(model=models)
        allnames <- unique(AllComb[["model"]])
        matVoCov <- matrix(0,length(allnames),length(allnames))
        rownames(matVoCov) <- allnames
        colnames(matVoCov) <- allnames    
        for(i in 1:nrow(AllComb)){    # First pair of model and time
            themodel1 <- AllComb[i,"model"]
            for(j in i:nrow(AllComb)){ # Second pair of model and time
                themodel2 <- AllComb[j,"model"]
                # extract the two iid decompositions
                IF1 <- data[model==themodel1][[IF.name]]
                IF2 <- data[model==themodel2][[IF.name]]
                # compute and save covariance
                matVoCov[i,j] <- cov(IF1,IF2)/N       
            }
        }
    }
    matVoCov[lower.tri(matVoCov)] <- t(matVoCov)[lower.tri(matVoCov)]
    matVoCov
}

######################################################################
### getVcov.R ends here
