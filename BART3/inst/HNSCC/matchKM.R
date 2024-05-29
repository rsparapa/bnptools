
library(MatchIt)
library(dbarts)
library(survival)
library(survminer)

library(BART3)
source(system.file('HNSCC/recode.R', package = 'BART3'))
##source("recode.R")

match = matchit(Study ~ Age + White + Black + Asian + Male +
	Smoker + Alcohol + PriorChemo + Larynx + OralCavity + Oropharynx +
	P16P + P16N + P16M + ECSY + ECSN + ECSM + hnscc$Positive_MY +
	hnscc$Positive_MN + hnscc$Positive_MM + hnscc$Close_MY + hnscc$Close_MN +
	hnscc$Close_MM + PNIY + PNIN + PNIM + LVIY + LVIN + LVIM +LN2Y + LN2N +
	LN2M + RSY + RSN + RSM + diag_year, data = hnscc, method = "nearest", 	
	distance = "bart", ratio = 2, caliper = 0.25, 
	distance.options = list(n.samples = 10000, n.burn = 5000, 
	n.thin = 10, seed = 2025))

mdata2 =  match.data(match)

# fit and plot KM curve

fitDFS <- survfit(Surv(DFS_months, DFS) ~ Study, data = mdata2)

jpeg("psmatch.jpg", width = 7, height = 7, units = 'in', res = 300)
ggsurvplot(fitDFS, mdata, 
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


