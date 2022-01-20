library(bench)

bm1<-mark(source("babysam.R"))
print(bm1)


bm2<-TMB::benchmark(obj)
print(bm2)

## and also for each line 
bmX<-workout_expressions(as.list(parse("./babysam.R")))
bm3 <- round(as.numeric(bmX$real), 5)
names(bm3)<-bmX$exprs
print(as.matrix(bm3))

# in addition functions to get every part by itself and a "press" function to evaluate
# some expressions for a range of input parameters (e.g. data sizes)
