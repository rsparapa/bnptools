library(BART3)
library(ggplot2)
data(hnscc)

# fit bart

times <- hnscc$DFS_Time
delta <- hnscc$DFS

x.train <- data.frame(hnscc[, -(1:2)])

## Windows does NOT support parallel processing by forking
if(.Platform$OS.type == 'windows') {
    set.seed(2025)
    post = surv.bart(x.train, times = times, delta = delta, x.test = x.train, 
	ndpost=1000, nskip=5000, keepevery=10)
} else {
    options(mc.cores = 8)
    post = mc.surv.bart(x.train, times = times, delta = delta, x.test = x.train, 
	ndpost=1000, nskip=5000, keepevery=10, seed = 2025)
}

# predict from bart
 
K <- post$K
M <- post$ndpost
N <- length(hnscc$DFS)

tx.test <- rbind(post$tx.test, post$tx.test)
tx.test[ , 'Study'] <- c(rep(1, N*K), rep(0, N*K))

pred <- predict(post, newdata=tx.test)

saveRDS(pred, "bart.rds")

pd <- matrix(nrow=M, ncol=2*K)

for(j in 1:K) {
    h <- seq(j, N*K, by=K)
    pd[ , j] <- apply(pred$surv.test[ , h], 1, mean) ## treated
    pd[ , j+K] <- apply(pred$surv.test[ , h+N*K], 1, mean) ## control
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

#2-year DFS rate (95%CI) for trial patients
trial = tail(subset(data.trial, time<=24),n=1)
c(trial$surv,trial$lower,trial$upper)

#2-year DFS rate (95%CI) for historical controls
hist = tail(subset(data.hist, time<=24),n=1)
c(hist$surv,hist$lower,hist$upper)

# difference in 2-year DFS rates

ind = nrow(subset(data.trial, time<=25)) - 1

diff24.mean = mean(pd[,ind] - pd[,ind+K])
diff24.bci = quantile(pd[,ind] - pd[,ind+K], prob=c(0.025,0.975))
c(diff24.mean,diff24.bci)


# plot bart predicted survival 

pdf("bart.pdf", width = 7, height = 7)
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
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          axis.title.x = element_text(size=18, hjust=.5, vjust=.25),
          axis.title.y = element_text(size=18, hjust=.5, vjust=1.5),
          plot.title = element_text(size=16, hjust=-.1, vjust=1)) +
     scale_y_continuous("Disease-free Survival Probability", limits=c(0, 1)) + 
     scale_x_continuous("Months from Surgery", limits=c(0, 36), breaks=seq(0, 36, 3)) +
  geom_point(aes(x=5, y=1), shape=95, col= "red", size=6) +
  annotate("text", x = 13, y = 1, label = "Historical Controls", size=7) +
  geom_point(aes(x=21, y=1), shape=95, col= "blue", size=6) +
  annotate("text", x = 27, y = 1, label = "Trial Patients", size=7) +
  geom_segment(aes(x = 24, y = 0, xend = 24, yend = trial$surv)) + 
  geom_segment(aes(x = 0, y = trial$surv, xend = 24, yend = trial$surv)) +
  geom_segment(aes(x = 0, y = hist$surv, xend = 24, yend = hist$surv))
dev.off()

