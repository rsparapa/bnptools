library(survival)
library(ggplot2)
library(BART3)
source(system.file('HNSCC/recode.R', package = 'BART3'))
write.csv(hnscc, 'hnscc.csv')
## run SAS program diradj.sas before continuing

diradj = read.csv("diradj.csv")

diradj$time = diradj$DFS_Time

diradj0000 = subset(diradj, P16N == 0 & P16M == 0 & PNIY == 0 & Study == 0)
diradj1000 = subset(diradj, P16N == 1 & P16M == 0 & PNIY == 0 & Study == 0)
diradj0100 = subset(diradj, P16N == 0 & P16M == 1 & PNIY == 0 & Study == 0)
diradj1010 = subset(diradj, P16N == 1 & P16M == 0 & PNIY == 1 & Study == 0)
diradj0110 = subset(diradj, P16N == 0 & P16M == 1 & PNIY == 1 & Study == 0)

diradj0001 = subset(diradj, P16N == 0 & P16M == 0 & PNIY == 0 & Study == 1)
diradj1001 = subset(diradj, P16N == 1 & P16M == 0 & PNIY == 0 & Study == 1)
diradj0101 = subset(diradj, P16N == 0 & P16M == 1 & PNIY == 0 & Study == 1)
diradj1011 = subset(diradj, P16N == 1 & P16M == 0 & PNIY == 1 & Study == 1)
diradj0111 = subset(diradj, P16N == 0 & P16M == 1 & PNIY == 1 & Study == 1)

jpeg("diradj000.jpg", width = 7, height = 7, units = 'in', res = 300)
ggplot() + 
  geom_step(data=diradj0001,mapping=aes(x=time,y=surv),color='blue',lwd=1) + 
  geom_step(data=diradj0000,mapping=aes(x=time,y=surv),color='red',lwd=1) + 
  geom_step(data=diradj0001,mapping=aes(x=time,y=lower),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj0000,mapping=aes(x=time,y=lower),color='red', lty=2, lwd=.5) + 
  geom_step(data=diradj0001,mapping=aes(x=time,y=upper),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj0000,mapping=aes(x=time,y=upper),color='red', lty=2, lwd=.5) + 
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
     geom_point(aes(x=25, y=0.4), shape=95, col= "blue", size=6) +
     annotate("text", x = 27, y = 0.35, label = "Trial Patients", size=4) +
     geom_point(aes(x=25, y=0.3), shape=95, col= "red", size=6) +
     annotate("text", x = 28, y = 0.25, label = "Historical Controls", size=4)
dev.off()

jpeg("diradj100.jpg", width = 7, height = 7, units = 'in', res = 300)
ggplot() + 
  geom_step(data=diradj1001,mapping=aes(x=time,y=surv),color='blue',lwd=1) + 
  geom_step(data=diradj1000,mapping=aes(x=time,y=surv),color='red',lwd=1) + 
  geom_step(data=diradj1001,mapping=aes(x=time,y=lower),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj1000,mapping=aes(x=time,y=lower),color='red', lty=2, lwd=.5) + 
  geom_step(data=diradj1001,mapping=aes(x=time,y=upper),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj1000,mapping=aes(x=time,y=upper),color='red', lty=2, lwd=.5) + 
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
     geom_point(aes(x=25, y=0.4), shape=95, col= "blue", size=6) +
     annotate("text", x = 27, y = 0.35, label = "Trial Patients", size=4) +
     geom_point(aes(x=25, y=0.3), shape=95, col= "red", size=6) +
     annotate("text", x = 28, y = 0.25, label = "Historical Controls", size=4)
dev.off()

jpeg("diradj010.jpg", width = 7, height = 7, units = 'in', res = 300)
ggplot() + 
  geom_step(data=diradj0101,mapping=aes(x=time,y=surv),color='blue',lwd=1) + 
  geom_step(data=diradj0100,mapping=aes(x=time,y=surv),color='red',lwd=1) + 
  geom_step(data=diradj0101,mapping=aes(x=time,y=lower),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj0100,mapping=aes(x=time,y=lower),color='red', lty=2, lwd=.5) + 
  geom_step(data=diradj0101,mapping=aes(x=time,y=upper),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj0100,mapping=aes(x=time,y=upper),color='red', lty=2, lwd=.5) + 
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
     geom_point(aes(x=25, y=0.4), shape=95, col= "blue", size=6) +
     annotate("text", x = 27, y = 0.35, label = "Trial Patients", size=4) +
     geom_point(aes(x=25, y=0.3), shape=95, col= "red", size=6) +
     annotate("text", x = 28, y = 0.25, label = "Historical Controls", size=4)
dev.off()

jpeg("diradj101.jpg", width = 7, height = 7, units = 'in', res = 300)
ggplot() + 
  geom_step(data=diradj1011,mapping=aes(x=time,y=surv),color='blue',lwd=1) + 
  geom_step(data=diradj1010,mapping=aes(x=time,y=surv),color='red',lwd=1) + 
  geom_step(data=diradj1011,mapping=aes(x=time,y=lower),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj1010,mapping=aes(x=time,y=lower),color='red', lty=2, lwd=.5) + 
  geom_step(data=diradj1011,mapping=aes(x=time,y=upper),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj1010,mapping=aes(x=time,y=upper),color='red', lty=2, lwd=.5) + 
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
     geom_point(aes(x=25, y=0.4), shape=95, col= "blue", size=6) +
     annotate("text", x = 27, y = 0.35, label = "Trial Patients", size=4) +
     geom_point(aes(x=25, y=0.3), shape=95, col= "red", size=6) +
     annotate("text", x = 28, y = 0.25, label = "Historical Controls", size=4)
dev.off()

jpeg("diradj011.jpg", width = 7, height = 7, units = 'in', res = 300)
ggplot() + 
  geom_step(data=diradj0111,mapping=aes(x=time,y=surv),color='blue',lwd=1) + 
  geom_step(data=diradj0110,mapping=aes(x=time,y=surv),color='red',lwd=1) + 
  geom_step(data=diradj0111,mapping=aes(x=time,y=lower),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj0110,mapping=aes(x=time,y=lower),color='red', lty=2, lwd=.5) + 
  geom_step(data=diradj0111,mapping=aes(x=time,y=upper),color='blue', lty=2, lwd=.5) + 
  geom_step(data=diradj0110,mapping=aes(x=time,y=upper),color='red', lty=2, lwd=.5) + 
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
     geom_point(aes(x=25, y=0.4), shape=95, col= "blue", size=6) +
     annotate("text", x = 27, y = 0.35, label = "Trial Patients", size=4) +
     geom_point(aes(x=25, y=0.3), shape=95, col= "red", size=6) +
     annotate("text", x = 28, y = 0.25, label = "Historical Controls", size=4)
dev.off()

## when P16N=P16M=0 there is only one patient with PNIY=1
## therefore, move that patient into the same stratum as PNIY=0
hnscc$PNIY[hnscc$P16N == 0 & hnscc$P16M == 0 & hnscc$PNIY == 1] <- 0

cfit <- coxph(Surv(DFS_Time, DFS) ~ 
	Age + White + Black + Asian + Male +
	Smoker + Alcoholic + PriorChemo + Larynx + OralCavity + Oropharynx +
	P16P + ECSY + ECSN + ECSM + Positve_MY +
	Positve_MN + Positve_MM + Close_MY + Close_MN +
	Close_MM + PNIN + PNIM + LVIY + LVIN + LVIM +LN2Y + LN2N +
	LN2M + RSY + RSN + RSM + Study + diag_year +
	strata(P16N) + strata(P16M) + strata(PNIY), data=hnscc)

#summary(cfit)   
cox.zph(cfit) 
  
          
# create survival function for each unique combination of covariates
tdata <- unique(hnscc[c("Age", "White", "Black", "Asian", "Male", 
	"Smoker", "Alcoholic", "PriorChemo", "Larynx", "OralCavity",
	"Oropharynx", "P16P", "P16N", "P16M", "ECSY", "ECSN", "ECSM",
	"Positve_MY", "Positve_MN", "Positve_MM", "Close_MY", "Close_MN", 
	"Close_MM", "PNIY", "PNIN", "PNIM", "LVIY", "LVIN", "LVIM",
	"LN2Y", "LN2N", "LN2M", "RSY", "RSN", "RSM", "diag_year")])
tdata$count = 1

tdata1 = tdata
tdata1$Study = 1

tdata0 = tdata
tdata0$Study = 0

est1 = survfit(cfit, tdata1)
pred1.24 = summary(est1, times = 24)
est1.surv.24 = mean(pred1.24$surv)
est1.lower.24 = mean(pred1.24$lower)
est1.upper.24 = mean(pred1.24$upper)
est1.se.24 = mean(pred1.24$std.err)

est0 = survfit(cfit, tdata0)
pred0.24 = summary(est0, times = 24)
est0.surv.24 = mean(pred0.24$surv)
est0.lower.24 = mean(pred0.24$lower)
est0.upper.24 = mean(pred0.24$upper)
est0.se.24 = mean(pred0.24$std.err)


#2-year DFS rate (95%CI) for trial patients
c(est1.surv.24,est1.lower.24,est1.upper.24)

#2-year DFS rate (95%CI) for historical controls
c(est0.surv.24,est0.lower.24,est0.upper.24)

# difference in 2-year DFS rates 
diff = est1.surv.24 - est0.surv.24
diffSE <- sqrt(est1.se.24^2 + est0.se.24^2)
diffL = diff - 1.96 *diffSE
diffU = diff + 1.96 *diffSE
c(diff, diffL, diffU)

###########################################

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

  


