library(TMB)
compile("kkr.cpp", framework="TMBad")#, flags="-g -O0") 
dyn.load(dynlib("kkr"))

load("kkr.RData")

#
dat$mode<-0 # dense
obj <- MakeADFun(dat, par, DLL="kkr", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
print(paste("Dense jnll:", obj$fn()))

dat$mode<-1 # generator
obj <- MakeADFun(dat, par, DLL="kkr", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
print(paste("Generator jnll:", obj$fn()))

dat$mode<-2 # hand written
obj <- MakeADFun(dat, par, DLL="kkr", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
print(paste("Hand written jnll:", obj$fn()))


dat$mode<-2 # hand written
obj <- MakeADFun(dat, par, random=c("logN", "logF", "missing"), DLL="kkr", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
print(paste("Hand written nll:", obj$fn()))

#dat$mode<-1 # generator
#obj <- MakeADFun(dat, par, random=c("logN", "logF", "missing"), DLL="kkr", map=list(logsdF=as.factor(rep(0,length(par$logsdF)))), silent=FALSE)
#print(paste("Generator nll:", obj$fn()))



