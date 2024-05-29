
library(survival)
library(ggplot2)

library(BART3)
source(system.file('HNSCC/recode.R', package = 'BART3'))
##source("recode.R")

cfit <- coxph(Surv(DFS_months, DFS) ~ 
Age + White + Black + Asian + Male +
	Smoker + Alcohol + PriorChemo + Larynx + OralCavity + Oropharynx +
	P16P + P16N + P16M + ECSY + ECSN + ECSM + Positive_MY +
	Positive_MN + Positive_MM + Close_MY + Close_MN +
	Close_MM + PNIY + PNIN + PNIM + LVIY + LVIN + LVIM +LN2Y + LN2N +
	LN2M + RSY + RSN + RSM + Study + diag_year, data=hnscc)

summary(cfit)   
cox.zph(cfit) 
          
# create survival function for each unique combination of covariates
tdata <- unique(hnscc[c("Age", "White", "Black", "Asian", "Male", 
	"Smoker", "Alcohol", "PriorChemo", "Larynx", "OralCavity",
	"Oropharynx", "P16P", "P16N", "P16M", "ECSY", "ECSN", "ECSM",
	"Positive_MY", "Positive_MN", "Positive_MM", "Close_MY", "Close_MN", 
	"Close_MM", "PNIY", "PNIN", "PNIM", "LVIY", "LVIN", "LVIM",
	"LN2Y", "LN2N", "LN2M", "RSY", "RSN", "RSM", "diag_year")])
tdata$count = 1

tdata1 = tdata
tdata1$Study = 1

tdata0 = tdata
tdata0$Study = 0

tdata2 = rbind(tdata1,tdata0)

est1 = survfit(cfit, tdata1)
est1.time = rowMeans(matrix(est1$time,nrow=length(est1$time)))
est1.surv = rowMeans(matrix(est1$surv,nrow=length(est1$time)))

est0 = survfit(cfit, tdata0)
est0.time = rowMeans(matrix(est0$time,nrow=length(est0$time)))
est0.surv = rowMeans(matrix(est0$surv,nrow=length(est0$time)))

trial.se = rep(0,108)
trial.lower = rep(0,108)
trial.upper = rep(0,108)
hist.se = rep(0,108)
hist.lower = rep(0,108)
hist.upper = rep(0,108)

for(i in 1:108) {
trial.se[i] <- sd(est1$surv[i,])/sqrt(length(est1$surv[i,]))
trial.lower[i] <- rowMeans(est1$surv)[i] - 1.96*trial.se[i]
trial.upper[i] <- rowMeans(est1$surv)[i] + 1.96*trial.se[i]

hist.se[i] <- sd(est0$surv[i,])/sqrt(length(est0$surv[i,]))
hist.lower[i] <- rowMeans(est0$surv)[i] - 1.96*hist.se[i]
hist.upper[i] <- rowMeans(est0$surv)[i] + 1.96*hist.se[i]

}

data.trial = data.frame(cbind(c(0,est1.time), c(1,est1.surv), c(0,trial.se), c(1,trial.lower), c(1,trial.upper)))
names(data.trial) = c("time","surv","se","lower","upper")

data.hist = data.frame(cbind(c(0,est0.time), c(1,est0.surv), c(0,hist.se), c(1,hist.lower), c(1,hist.upper)))
names(data.hist) = c("time","surv","se","lower","upper")


jpeg("coxph.jpg", width = 7, height = 7, units = 'in', res = 300)
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
diff = trial$surv - hist$surv
diffSE <- sqrt(trial$se^2 + hist$se^2)
diffL = diff - 1.96 *diffSE
diffU = diff + 1.96 *diffSE
c(diff, diffL, diffU)

  


