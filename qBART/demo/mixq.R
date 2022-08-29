library(qBART)
data(mixcure)
x.train <- mixcure[,1:10]
time <- mixcure$obstime
event <- mixcure$event
timegrid <- c(10,20)
fit_train <- mc.qbart(x.train1=x.train, x.train2=x.train, times=time, delta=event, grid=timegrid, sparse=TRUE, mc.core=8)
fit_all <- mc.qbart(x.train1=x.train, x.train2=x.train, times=time, delta=event, sparse=TRUE, mc.core=8)
fit_test <- mc.qbart(x.train1=x.train, x.train2=x.train, times=time, delta=event, grid=timegrid, x.test1=x.train, x.test2=x.train, sparse=TRUE, mc.core=8)
