library(survival)
library(ggplot2)

data(hnscc, package = 'BART3')

cfit <- coxph(Surv(DFS_Time, DFS) ~ 
                Age + Gender + Race + Smoking + Alcohol + Prior_Chemo +         
                Primary_Disease_Site + ECS +  Positve_M + PNI + LVI + 
                Risk_Stratification + Lymph_Nodes_Involved_2 + 
                diag_year + Study, data=hnscc)

#summary(cfit)   
cox.zph(cfit) 


# create survival function for each unique combination of covariates

tdata1 = hnscc
tdata1$Study = 1

tdata0 = hnscc
tdata0$Study = 0

est1 = survfit(cfit, tdata1)
est1.time = est1$time
est1.surv = rowMeans(matrix(est1$surv,nrow=length(est1$time)))


est0 = survfit(cfit, tdata0)
est0.time = est0$time
est0.surv = rowMeans(matrix(est0$surv,nrow=length(est0$time)))


trial.se = rep(0,90)
trial.lower = rep(0,90)
trial.upper = rep(0,90)
hist.se = rep(0,90)
hist.lower = rep(0,90)
hist.upper = rep(0,90)

for(i in 1:90) {
  trial.se[i] <- sd(est1$surv[i,])/sqrt(length(est1$surv[i,]))
  trial.lower[i] <- rowMeans(est1$surv)[i] - 1.96*trial.se[i]
  trial.upper[i] <- rowMeans(est1$surv)[i] + 1.96*trial.se[i]
  
  hist.se[i] <- sd(est0$surv[i,])/sqrt(length(est0$surv[i,]))
  hist.lower[i] <- rowMeans(est0$surv)[i] - 1.96*hist.se[i]
  hist.upper[i] <- rowMeans(est0$surv)[i] + 1.96*hist.se[i]
  
}

data.trial = data.frame(cbind(c(0,est1.time), c(1,est1.surv), c(1,trial.lower), c(1,trial.upper), c(0,trial.se)))
names(data.trial) = c("time","surv","lower","upper","se")

data.hist = data.frame(cbind(c(0,est0.time), c(1,est0.surv), c(1,hist.lower), c(1,hist.upper), c(0,hist.se)))
names(data.hist) = c("time","surv","lower","upper","se")

#2-year DFS rate (95%CI) for trial patients
trial = tail(subset(data.trial, time<=24),n=1)
c(trial$surv,trial$lower,trial$upper)

#2-year DFS rate (95%CI) for historical controls
hist = tail(subset(data.hist, time<=24),n=1)
c(hist$surv,hist$lower,hist$upper)

# difference in 2-year DFS rates

## the 24 month time point is exceeded in interval 66 (24.0473062)

#predicted counterfactuals at 24 months (ITE)
Y = c(est1$surv[65,], est0$surv[65,])
#treatment indicator
Z = c(rep(1,length(est1$surv[65,])), rep(0,length(est1$surv[65,])))

#marginal structural model for ATE
msm = lm(Y ~ Z)
round(c(summary(msm)$coef[2,1], confint(msm)[2,]),3)


pdf("coxph.pdf", width = 7, height = 7)
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

## the 24 month time point is exceeded in interval 66 (24.0473062)

surv.est1 <- est1$surv[65,]
surv.est0 <- est1$surv[65,]

## start with survival
ITE1 <- surv.est1
ITE0 <- surv.est0

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


