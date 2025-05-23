
library(randomForestSRC)
library(ggplot2)

data(hnscc, package = 'BART3')

# fit rsf
rf <- rfsrc(Surv(DFS_Time, DFS) ~ Age + Gender + Race + Smoking + Alcohol + Prior_Chemo +         
              Primary_Disease_Site + ECS +  Positve_M + PNI + LVI + 
              Risk_Stratification + Lymph_Nodes_Involved_2 + 
              diag_year + Study, data = hnscc, seed = 2025)

times <- rf$time.interest

# create survival function for each unique combination of covariates
tdata1 = hnscc
tdata1$Study = 1

tdata0 = hnscc
tdata0$Study = 0
             
N <- length(hnscc$DFS)

newdata = rbind(tdata1,tdata0)
pred = predict(rf,newdata)

pd = pred$survival
times = pred$time.interest

trial <- 1:N
hist <- trial+N

pd.mu.trial  <- apply(pd[trial,], 2, mean)
pd.025.trial <- apply(pd[trial,], 2, quantile, probs=0.025)
pd.975.trial <- apply(pd[trial,], 2, quantile, probs=0.975)

pd.mu.hist  <- apply(pd[hist,], 2, mean)
pd.025.hist <- apply(pd[hist,], 2, quantile, probs=0.025)
pd.975.hist <- apply(pd[hist,], 2, quantile, probs=0.975)

# difference in 2-year DFS rates

## the 24 month time point is exceeded in interval 47 (24.671485)

#predicted counterfactuals at 24 months (ITE)
Y = c(pd[trial,46], pd[hist,46])
#treatment indicator
Z = c(rep(1,length(pd[trial,46])), rep(0,length(pd[hist,46])))

#marginal structural model for ATE
msm = lm(Y ~ Z)
round(c(summary(msm)$coef[2,1], confint(msm)[2,]),3)


data.trial = data.frame(cbind(c(0,times), c(1,pd.mu.trial),
c(1,pd.025.trial), c(1,pd.975.trial)))
names(data.trial) = c("time","surv","lower","upper")

data.hist = data.frame(cbind(c(0,times), c(1,pd.mu.hist),
c(1,pd.025.hist), c(1,pd.975.hist)))
names(data.hist) = c("time","surv","lower","upper")

#2-year DFS rate (95%CI) for trial patients
trial = tail(subset(data.trial, time<=24),n=1)
c(trial$surv,trial$lower,trial$upper)

#2-year DFS rate (95%CI) for historical controls
hist = tail(subset(data.hist, time<=24),n=1)
c(hist$surv,hist$lower,hist$upper)

# difference
##c(pd.diff.mu.trial, pd.diff.025.trial, pd.diff.975.trial)


pdf("rsf.pdf", width = 7, height = 7)
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


## calculate ATE at 24 months ##

## the 24 month time point is exceeded in interval 47 (24.671485)
surv.est <- pd[,46]

## start with survival
ITE1 <- surv.est[1:N]
ITE0 <- surv.est[-(1:N)]

observed <- cbind(hnscc$Study == 1, hnscc$Study == 0)

## account for time=24 of observed potential outcomes
failure24 <- (hnscc$DFS_Time <= 24 & hnscc$DFS == 1)
table(failure24, useNA = 'ifany')

success24 <- (hnscc$DFS_Time >= 24)
table(success24, useNA = 'ifany')

table(failure24 & success24, useNA = 'ifany') ## you cannot be both

ITE1[observed[ , 1] & failure24] <- 0
ITE1[observed[ , 1] & success24] <- 1

ITE0[observed[ , 2] & failure24] <- 0
ITE0[observed[ , 2] & success24] <- 1

#predicted counterfactuals at 24 months (ITE)
Y = c(ITE1, ITE0)
#treatment indicator
Z = c(rep(1,length(ITE1)), rep(0,length(ITE0)))

#marginal structural model for ATE
msm = lm(Y ~ Z)
round(c(summary(msm)$coef[2,1], confint(msm)[2,]),3)







