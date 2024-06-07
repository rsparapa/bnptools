
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
(P <- P-2)

post <-  list()
H <- c(5, 10, 20)
for(i in 1:3)
    post[[i]] <- gbart(x.train, NCItopo$Potency[train], 
                       ntree = H[i], sparse = TRUE, type = 'pbart')

for(i in 1:3) {
    if(i == 1) 
        plot(post[[i]]$varprob.mean, type = 'p', 
             xlim = c(-10, 280), ylim = c(0, 0.275), ylab = 'probability')
    print(summary(post[[i]]$varprob.mean))
    points(post[[i]]$varprob.mean, col = 2^(i-1))
    x <- which(post[[i]]$varprob.mean>(1/P))
    y <- post[[i]]$varprob.mean[x]
    for(j in 1:length(x)) lines(c(x[j], x[j]), c(1/P, y[j]), col = 2^(i-1))
    text(x, y, dimnames(x.train)[[2]][x], col = 2^(i-1), pos = i+1)
}
abline(h = c(0, 1/P), col = 8)
legend('topleft', legend = H, col = 2^(0:2), pch = 1)
##dev.copy2pdf(file = 'nci-topo.pdf')

