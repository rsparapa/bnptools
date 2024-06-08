
library(BART3)
data(NCItopo)

set.seed(10)
NCItopo$Potency[NCItopo$Potency == 2] <- 1 ## combine 1:2
table(NCItopo$Potency) 
(P <- ncol(NCItopo))
active <- which(NCItopo$Potency == 1)
inactive <- which(NCItopo$Potency == 0)
(M <- length(active))
pick <- sample(inactive, M) 
train <- c(active, pick)
table(NCItopo$Potency[train])
x.train <- NCItopo[train, -c(1, P)] ## the first column is an id
                                    ## and the last col is y 

post <-  list()
H <- c(5, 10)
H. <- length(H)
##H <- c(5, 10, 20)
for(l in 1:2) 
    for(i in 1:H.) {
        m <- H.*(l-1)+i
        post[[m]] <- gbart(x.train, NCItopo$Potency[train], 
                       ntree = H[i], sparse = (l == 1), type = 'pbart')
        if(m == 1) {
            x.train <- post[[1]]$x.train
            (P <- ncol(x.train)) 
            ## due to the relatively small sample, some TII are constant
            ## by default, these are dropped since no splits are observed
            ## see the rm.const argument
        }
    }

sub = c('DART', 'BART')
ylim <- c(0.275, 0.10)
top <- 5
par(mfrow = c(2, 2))
for(l in 1:2) {
    for(i in 1:H.) {
        m <- H.*(l-1)+i
        if(l == 2) 
post[[m]]$varprob.mean <- post[[m]]$varcount.mean/sum(post[[m]]$varcount.mean)
## for BART, varprob is simply the constant 1/P
## so we will use varcount normalized between 0 and 1 
        if(i == 1) {
            plot(post[[m]]$varprob.mean, type = 'p', 
                 xlim = c(-10, 280), ylim = c(0, ylim[l]),
                 xlab = 'variable index', sub = sub[l],
                 ylab = 'selection probability')
            abline(h = (0:2)/P, col = 8)
        }
        print(summary(post[[m]]$varprob.mean))
        points(post[[m]]$varprob.mean, col = 2^(i-1))
        x <- which(post[[m]]$varprob.mean>(2/P))
        y <- post[[m]]$varprob.mean[x]
        ##for(j in 1:length(x)) lines(c(x[j], x[j]), c(1/P, y[j]), col = 2^(i-1))
        text(x, y, dimnames(x.train)[[2]][x], col = 2^(i-1), pos = i+1)
    }
    for(i in 1:H.) {
        m <- H.*(l-1)+i
        h <- order(post[[m]]$varprob.mean, decreasing = TRUE)[1:top]
        if(i == 1) { 
            plot(post[[m]]$varprob.mean[h], type = 'p', 
                 xlim = c(-top/5, 1.2*top), ylim = c(0, ylim[l]),
                 xlab = 'ordered variable index', sub = sub[l],
                 ylab = 'selection probability')
            abline(h = (0:2)/P, col = 8)
            if(m == 1)
                legend('topright', legend = H, col = 2^(0:2), pch = 1)
        }
        points(post[[m]]$varprob.mean[h], col = 2^(i-1))
        x <- which(post[[m]]$varprob.mean[h]>(2/P))
        y <- post[[m]]$varprob.mean[h][x]
        for(j in 1:length(x)) {
            ##lines(c(x[j], x[j]), c(1/P, y[j]), col = 2^(i-1))
            text(x[j], y[j], dimnames(x.train)[[2]][h][x[j]], 
                 col = 2^(i-1), pos = (j%%3)+1)
        }
    }
}
par(mfrow = c(1, 1))
dev.copy2pdf(file = 'nci-topo.pdf')

