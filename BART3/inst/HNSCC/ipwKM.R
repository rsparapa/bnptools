
library(MatchIt)
library(dbarts)
library(ggplot2)

library(BART3)
source(system.file('HNSCC/recode.R', package = 'BART3'))
##source("recode.R")

match = matchit(Study ~ Age + White + Black + Asian + Male +
	Smoker + Alcohol + PriorChemo + Larynx + OralCavity + Oropharynx +
	P16P + P16N + P16M + ECSY + ECSN + ECSM + Positive_MY +
	Positive_MN + Positive_MM + Close_MY + Close_MN +
	Close_MM + PNIY + PNIN + PNIM + LVIY + LVIN + LVIM +
        LN2Y + LN2N + LN2M + RSY + RSN + RSM + diag_year, 
        data = hnscc, method = "nearest", 
	distance = "bart", ratio = 2, caliper = 0.25, 
	distance.options = list(n.samples = 10000, n.burn = 5000, 
	n.thin = 10, seed = 2025))

mdata = match.data(match, drop.unmatched=F)

# modified adjusted.KM function from RISCA package to include se
source(system.file('HNSCC/adjKM.R', package = 'BART3'))

# fit adjusted KM estimator by IPTW
ps <- mdata$distance
w <- (mdata$Study==1)*(1/ps) + (mdata$Study==0)*(1)/(1-ps)

res.akm <- adjusted.KM.se(times=mdata$DFS_months, failures=mdata$DFS,
 variable=mdata$Study, weights=w)


# plot adjusted KM curve

data.trial = data.frame(res.akm$times[res.akm$variable==1], res.akm$survival[res.akm$variable==1],
res.akm$survival[res.akm$variable==1]+1.96*res.akm$std.err[res.akm$variable==1],
res.akm$survival[res.akm$variable==1]-1.96*res.akm$std.err[res.akm$variable==1])
names(data.trial) = c("time","surv","upper","lower")

data.hist = data.frame(res.akm$times[res.akm$variable==0], res.akm$survival[res.akm$variable==0],
res.akm$survival[res.akm$variable==0]+1.96*res.akm$std.err[res.akm$variable==0],
res.akm$survival[res.akm$variable==0]-1.96*res.akm$std.err[res.akm$variable==0])
names(data.hist) = c("time","surv","upper","lower")

jpeg("iptw.jpg", width = 7, height = 7, units = 'in', res = 300)
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
trial = tail(subset(res.akm, variable==1 & times<=24),n=1)
trial.surv = trial$survival
trial.lower = trial$survival - 1.96*trial$std.err
trial.upper = trial$survival + 1.96*trial$std.err
c(trial.surv,trial.lower,trial.upper)

#2-year DFS rate (95%CI) for historical controls
hist = tail(subset(res.akm, variable==0 & times<=24),n=1)
hist.surv = hist$survival
hist.lower = hist$survival - 1.96*hist$std.err
hist.upper = hist$survival + 1.96*hist$std.err
c(hist.surv,hist.lower,hist.upper)

# difference in 2-year DFS rates 

diff = trial.surv - hist.surv
diffSE <- sqrt(trial$std.err^2 + hist$std.err^2)
diffL = diff - 1.96 *diffSE
diffU = diff + 1.96 *diffSE
c(diff, diffL, diffU)


  
