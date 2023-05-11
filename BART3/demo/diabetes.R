
options(mc.cores=8)
library(BART3)

## full data set
data(xdm.train)
data(ydm.train)

check = tsvs(xdm.train, ydm.train, type='pbart', keepevery=1, ndpost=2500,
             rds.file='diabetes.rds', pdf.file='diabetes.pdf')


