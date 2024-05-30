
library(ggplot2)

options(mc.cores = 8)
library(BART3)
source(system.file('HNSCC/recode.R', package = 'BART3'))
##source("recode.R")

# fit bart

times <- hnscc$DFS_months
delta <- hnscc$DFS

x.train <-  hnscc[ , -(1:2)]

## chains run in parallel
## post = mc.surv.bart(x.train, times = times, delta = delta, 
## x.test = x.train, ndpost=1000, nskip=5000, keepevery=10, seed = 2025)

## a single serial chain
set.seed(2025)
post = surv.bart(x.train, times = times, delta = delta, x.test = x.train, 
                 ndpost=1000, nskip=5000, keepevery=10)

pre <- surv.pre.bart(times = times, delta = delta, x.train = x.train, 
                     x.test = x.train)

# predict from bart
 
K <- pre$K
M <- post$ndpost
N <- length(hnscc$DFS)

pre$tx.test <- rbind(pre$tx.test, pre$tx.test)
pre$tx.test[ , 2] <- c(rep(1, N*K), rep(0, N*K))

pred <- predict(post, newdata=pre$tx.test)

pd <- matrix(nrow=M, ncol=2*K)

for(j in 1:K) {
    h <- seq(j, N*K, by=K)
    pd[ , j] <- apply(pred$surv.test[ , h], 1, mean)
    pd[ , j+K] <- apply(pred$surv.test[ , h+N*K], 1, mean)
}

pd.mu  <- apply(pd, 2, mean)
pd.025 <- apply(pd, 2, quantile, probs=0.025)
pd.975 <- apply(pd, 2, quantile, probs=0.975)

trial <- 1:K
hist <- trial+K

data.trial = data.frame(cbind(c(0,pred$times), c(1,pd.mu[trial]),
c(1,pd.025[trial]), c(1,pd.975[trial])))
names(data.trial) = c("time","surv","lower","upper")

data.hist = data.frame(cbind(c(0,pred$times), c(1,pd.mu[hist]),
c(1,pd.025[hist]), c(1,pd.975[hist])))
names(data.hist) = c("time","surv","lower","upper")

# plot bart predicted survival 

jpeg("bart.jpg", width = 7, height = 7, units = 'in', res = 300)
ggplot() + 
  geom_step(data=data.trial,mapping=aes(x=time,y=surv),color='blue',lwd=1) + 
  geom_step(data=data.hist,mapping=aes(x=time,y=surv),color='red',lwd=1) + 
  geom_step(data=data.trial,mapping=aes(x=time,y=lower),color='blue', lty=2, lwd=.5) + 
  geom_step(data=data.hist,mapping=aes(x=time,y=lower),color='red', lty=2, lwd=.5) + 
  geom_step(data=data.trial,mapping=aes(x=time,y=upper),color='blue', lty=2, lwd=.5) + 
  geom_step(data=data.hist,mapping=aes(x=time,y=upper),color='red', lty=2, lwd=.5) + 
  theme_bw() +
     theme(
          plot.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.position="top",
          axis.line = element_line(color = 'black'),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=16, hjust=.5, vjust=.25),
          axis.title.y = element_text(size=16, hjust=.5, vjust=1.5),
          plot.title = element_text(size=16, hjust=-.1, vjust=1)) +
     scale_y_continuous("Disease-free Survival Probability", limits=c(0, 1)) + 
     scale_x_continuous("Months from Surgery", limits=c(0, 33), breaks=seq(0, 33, 3)) +
     geom_point(aes(x=10, y=1), shape=95, col= "blue", size=6) +
     annotate("text", x = 15, y = 1, label = "Trial Patients", size=5) +
     geom_point(aes(x=21, y=1), shape=95, col= "red", size=6) +
     annotate("text", x = 27, y = 1, label = "Historical Controls", size=5)
dev.off()

#2-year DFS rate (95%CI) for trial patients
trial = tail(subset(data.trial, time<=24),n=1)
c(trial$surv,trial$lower,trial$upper)

#2-year DFS rate (95%CI) for historical controls
hist = tail(subset(data.hist, time<=24),n=1)
c(hist$surv,hist$lower,hist$upper)

# difference in 2-year DFS rates

ind = nrow(subset(data.trial, time<=24)) - 1

diff24.mean = mean(pd[,ind] - pd[,ind+K])
diff24.bci = quantile(pd[,ind] - pd[,ind+K], prob=c(0.025,0.975))
c(diff24.mean,diff24.bci)
 
# waterfall plot for difference and plot superiority

pd24 <- matrix(nrow=M, ncol=2*K)

for(j in 1:K) {
    h <- seq(j*67, N*K, by=K)
    pd24[ , j] <- apply(pred$surv.test[ , h], 1, mean)
    pd24[ , j+K] <- apply(pred$surv.test[ , h+N*K], 1, mean)

}


trial <- 1:K
hist <- trial+K

pd24.trial = pd24[,trial]
pd24.hist = pd24[,hist]

diff24 = pd24.trial - pd24.hist
sup24 = I(pd24.trial - pd24.hist > 0.20)

diff24.mu  <- apply(diff24, 2, mean)
diff24.025 <- apply(diff24, 2, quantile, probs=0.025)
diff24.975 <- apply(diff24, 2, quantile, probs=0.975)
diff24.25 <- apply(diff24, 2, quantile, probs=0.25)
diff24.75 <- apply(diff24, 2, quantile, probs=0.75)


data.diff24 = cbind(diff24.mu,diff24.025,diff24.975,diff24.25,diff24.75)
data.diff24 = data.diff24[order(data.diff24[,1]), ]

pt = seq(1,108)

sup24.mu  <- sort(apply(sup24, 2, mean))


jpeg("wf.jpg", width = 7, height = 7, units = 'in', res = 300)

plot(pt,data.diff24[,1], type="l", ylim = c(0,.5), xlab = "Patient Number (Ordered)", 
ylab = "Difference in 2-year Disease-free Survival Probabilities")
polygon(c(pt,rev(pt)),c(data.diff24[,2],rev(data.diff24[,3])),col = "grey80", border = FALSE)
polygon(c(pt,rev(pt)),c(data.diff24[,4],rev(data.diff24[,5])),col = "grey70", border = FALSE)
lines(pt,data.diff24[,1], lwd=2)
legend("topleft", legend=c("95% Credible Interval", "Inter-quartile Range"), fill=c("gray80", "gray70"), bty="n") 
legend("topright", legend="Posterior Mean", lty = 1, lwd = 2, bty="n") 
dev.off()


jpeg("sup.jpg",
width = 7, height = 7, units = 'in', res = 300)

plot(pt,sup24.mu, type="l", ylim = c(0,1), lwd = 2,
xlab = "Patient Number (Ordered)", 
ylab = "Superiority in 2-year Disease-free Survival Probabilities")
dev.off()




