library(stockassessment)
#fit<-fitfromweb("NEA_sei_21_v5Reca")
load("fit.RData")

dat<-list()
dat$obs <- exp(fit$data$logobs)
dat$aux <- fit$data$aux
dat$idx1 <- fit$data$idx1
dat$idx2 <- fit$data$idx2
dat$minYear <- min(fit$data$years)
dat$minAge <- min(fit$data$minAgePerFleet)
dat$fleetTypes <- fit$data$fleetTypes
dat$sampleTimes <- fit$data$sampleTimes
dat$year <- fit$data$years
dat$age <-  min(fit$data$minAgePerFleet):max(fit$data$maxAgePerFleet)   
dat$M <- fit$data$natMor
dat$SW <- fit$data$stockMeanWeight
dat$MO <- fit$data$propMat
dat$PF <- fit$data$propF
dat$PM <- fit$data$propM

dat$srmode <- 2
dat$fcormode <- 2
dat$keyF <- fit$conf$keyLogFsta[1,]
dat$keyQ <- fit$conf$keyLogFpar
dat$keySd <- fit$conf$keyVarObs
dat$keySd[dat$keySd<(-.1)] <- NA
dat$covType <- c(0,1,2) #as.integer(fit$conf$obsCorStruct)-1
dat$keyIGAR <- fit$conf$keyCorObs
dat$keyIGAR[fit$conf$keyCorObs==-1] <- NA
dat$keyIGAR[is.na(fit$conf$keyCorObs)] <- -1
dat$keyIGAR[2, 1:4]<-0
dat$noParUS <- sapply(1:length(dat$fleetTypes),
                      function(f){
                        A <- sum(!is.na(dat$keySd[f,]))
                        ifelse(dat$covType[f]==2, (A*A-A)/2, 0)
                      })
dat$mode<-2
#dat$obs[c(7,9,13)+100]<-NA

library(TMB)
compile("babycam.cpp", framework="TMBad")#, flags="-g -O0") 
dyn.load(dynlib("babycam"))

par <- list()
par$logsdR <-0
par$logsdS <- 0
par$logsdF <- numeric(max(dat$keyF)+1)
par$rickerpar <- if(dat$srmode==1){c(1,1)}else{numeric(0)}
par$transRhoF <- if(dat$fcormode==0){numeric(0)}else{0.1}
par$bhpar <- if(dat$srmode==2){c(1,1)}else{numeric(0)}
par$logQ <- numeric(max(dat$keyQ, na.rm=TRUE)+1)
par$logsd <- numeric(max(dat$keySd, na.rm=TRUE)+1)
par$logIGARdist <- numeric(max(dat$keyIGAR, na.rm=TRUE)+1)
par$parUS <- numeric(sum(dat$noParUS))
par$logN <- matrix(0, nrow=length(dat$year), ncol=length(dat$age))
par$logF <- matrix(0, nrow=length(dat$year), ncol=max(dat$keyF)+1)
par$missing <- numeric(sum(is.na(dat$obs)))

#
dat$mode<-0 # dense
obj <- MakeADFun(dat, par, DLL="babycam", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
obj$fn()
#g1<-obj$gr()

dat$mode<-1 # generator
obj <- MakeADFun(dat, par, DLL="babycam", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
obj$fn()
#g2<-obj$gr()

dat$mode<-2 # hand written
obj <- MakeADFun(dat, par, DLL="babycam", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
obj$fn()
#g3<-obj$gr()

#dat$mode<-1 # generator
#obj <- MakeADFun(dat, par, random=c("logN", "logF", "missing"), DLL="babycam", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
#TMB:::checkSparseHessian(obj)


                                        #stop("done")
dat$mode<-1 # Kasper
obj <- MakeADFun(dat, par, random=c("logN", "logF", "missing"), DLL="babycam", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000))

sdr<-sdreport(obj)
pl<-as.list(sdr, "Est")

#arg<-list(data=dat, parameters=par, random=c("logN", "logF", "missing"), DLL="babysam", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=TRUE)

