
library(BART3)

data(NCItopo)

##NCItopo <- read.csv(file = 'NCI-topo.csv')
##save(NCItopo, file = '~/git/bnptools/BART3/data/NCItopo.rda')

table(NCItopo$Potency) 

NCItopo$Potency[NCItopo$Potency == 2] <- 1 ## combine 1:2

table(NCItopo$Potency) 

(P <- ncol(NCItopo))

H = 10                          ## number of trees
T = 20                          ## number of steps
x.train <- NCItopo[ , -c(1, P)] ## the first column is an id
                                ## and the last col is y

sparse <- tsvs(x.train, NCItopo$Potency, keepevery = 1, ntree = H, T = T,
               type = 'pbart',
               rds.file = 'nci-topo.rds', pdf.file = 'nci-topo.pdf')

