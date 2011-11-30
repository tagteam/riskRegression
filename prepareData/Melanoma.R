Melanoma <- read.csv2("~/research/SoftWare/timereg/pkg/riskRegression/prepareData/Melanoma.csv")
names(Melanoma) <- c("event","time","invasion","ici","epicel","ulcer","thick","sex","age")
Melanoma$event=factor(Melanoma$event,levels=c(2,1,3),labels=c("censored","death.malignant.melanoma","death.other.causes"))
Melanoma$status=as.numeric(Melanoma$event)-1
Melanoma=Melanoma[,c("time","status","event","invasion","ici","epicel","ulcer","thick","sex","age")]
Melanoma$thick=as.numeric(gsub(",",".",Melanoma$thick))
Melanoma$epicel=factor(Melanoma$epicel,levels=c(1,2),labels=c("not present","present"))
Melanoma$invasion=factor(Melanoma$invasion,levels=c(0,1,2),labels=c("level.0","level.1","level.2"))
Melanoma$ulcer=factor(Melanoma$ulcer,levels=c(0,1),labels=c("not present","present"))
Melanoma$sex=factor(Melanoma$sex,levels=c(0,1),labels=c("Female","Male"))
save(Melanoma,file="~/research/SoftWare/timereg/pkg/riskRegression/data/Melanoma.rda")
