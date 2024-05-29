
library(survival)
library(survminer)

library(BART3)
source(system.file('HNSCC/recode.R', package = 'BART3'))
##source("recode.R")

# fit and plot KM curve

fitDFS <- survfit(Surv(DFS_months, DFS) ~ Study, data = hnscc)

jpeg("km.jpg", width = 7, height = 7, units = 'in', res = 300)
ggsurvplot(fitDFS, hnscc, 
         xlab = "Months from Surgery", cex.xlab = .4, legend.title="",
	legend.labs=c("Trial Patients", "Historical Controls"), xlim=c(0,33),
     	ylab = "Disease-free Survival Probability", cex.ylab = 1,
	font.x =  16, font.y = 16, font.tickslab = 14, font.legend = 14,
	risk.table = F, pval = FALSE, 
	    conf.int = TRUE, conf.int.style = "step", palette = c("blue","red"), break.time.by = 3)
dev.off()

# obtain 2-year DFS rates and difference

twoyr = summary(fitDFS,time=24)
twoyr

diff = twoyr$surv[2] - twoyr$surv[1]
diffSE <- sqrt(twoyr$std.err[2]^2 + twoyr$std.err[1]^2)
diffL = diff - 1.96 *diffSE
diffU = diff + 1.96 *diffSE
c(diff, diffL, diffU)


