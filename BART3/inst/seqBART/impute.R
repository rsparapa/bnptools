
system('make clean; make pbartall')

library(BART3)
library(MASS)
library(Rlab)
library(msm)

(seqBART.dir <- paste0(system.file(package = 'BART3'), '/seqBART/'))
(data.dir <- 'DIRECTORY/')

source(paste0(seqBART.dir, 'load-pbart.R'))

x.train <- readRDS(paste0(data.dir, 'ex10-x.train.rds'))
y.train <- readRDS(paste0(data.dir, 'ex10-y.train.rds'))

set.seed(21)
P <- ncol(x.train)

imputed<-BartMI(bartmpicstem, xx=x.train, yy=y.train, datatype=rep(0, P),
               type=1, numskip=199, burn=100)

saveRDS(imputed, paste0(data.dir, 'ex10-imputed.rds'))
