
library(randomForestSRC)
library(ggplot2)

library(BART3)
source(system.file('HNSCC/recode.R', package = 'BART3'))
##source("recode.R")

# fit rsf

rfdata = hnscc[,c("DFS_months", "DFS", "Age", "White", "Black", "Asian", "Male", 
"Smoker", "Alcohol", "PriorChemo", "Larynx", "OralCavity", "Oropharynx", "P16P", "P16N", "P16M", 
"ECSY", "ECSN", "ECSM", "Positive_MY", "Positive_MN", "Positive_MM", "Close_MY", "Close_MN", "Close_MM", 
"PNIY", "PNIN", "PNIM", "LVIY", "LVIN", "LVIM", "LN2Y", "LN2N", "LN2M", 
"RSY", "RSN", "RSM", "diag_year","Study")]


rf <- rfsrc(Surv(DFS_months, DFS) ~ Age + White + Black + Asian + Male +
	Smoker + Alcohol + PriorChemo + Larynx + OralCavity + Oropharynx +
	P16P + P16N + P16M + ECSY + ECSN + ECSM + Positive_MY + Positive_MN + Positive_MM + 
	Close_MY + Close_MN + Close_MM + PNIY + PNIN + PNIM + LVIY + LVIN + LVIM +LN2Y + LN2N +
	LN2M + RSY + RSN + RSM + diag_year + Study, data = rfdata, 
	nsplit = 20, ntree = 200, seed = 2025)

times <- rf$time.interest

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
             
N <- length(hnscc$DFS)

newdata = rbind(tdata1,tdata0)
pd = predict(rf,newdata)$survival

trial <- 1:N
hist <- trial+N

pd.mu.trial  <- apply(pd[trial,], 2, mean)
pd.025.trial <- apply(pd[trial,], 2, quantile, probs=0.025)
pd.975.trial <- apply(pd[trial,], 2, quantile, probs=0.975)

pd.mu.hist  <- apply(pd[hist,], 2, mean)
pd.025.hist <- apply(pd[hist,], 2, quantile, probs=0.025)
pd.975.hist <- apply(pd[hist,], 2, quantile, probs=0.975)


data.trial = data.frame(cbind(c(0,times), c(1,pd.mu.trial),
c(1,pd.025.trial), c(1,pd.975.trial)))
names(data.trial) = c("time","surv","lower","upper")

data.hist = data.frame(cbind(c(0,times), c(1,pd.mu.hist),
c(1,pd.025.hist), c(1,pd.975.hist)))
names(data.hist) = c("time","surv","lower","upper")

jpeg("rsf.jpg", width = 7, height = 7, units = 'in', res = 300)
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
trial.24 = tail(subset(data.trial, time<=24),n=1)
c(trial.24$surv,trial.24$lower,trial.24$upper)

#2-year DFS rate (95%CI) for historical controls
hist.24 = tail(subset(data.hist, time<=24),n=1)
c(hist.24$surv,hist.24$lower,hist.24$upper)

# difference
pd.diff = pd[trial,] - pd[hist,]
pd.diff.mu.trial  <- mean(pd.diff)
pd.diff.025.trial <- quantile(pd.diff, probs=0.025)
pd.diff.975.trial <- quantile(pd.diff, probs=0.975)
c(pd.diff.mu.trial, pd.diff.025.trial, pd.diff.975.trial)






