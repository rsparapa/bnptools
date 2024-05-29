
library(MatchIt)
library(dbarts)
library(randomForest)

library(BART3)
source(system.file('HNSCC/recode.R', package = 'BART3'))
##source("recode.R")

#use BART for calculating propensity score

match = matchit(Study ~ Age + White + Black + Asian + Male +
	Smoker + Alcohol + PriorChemo + Larynx + OralCavity + Oropharynx +
	P16P + P16N + P16M + ECSY + ECSN + ECSM + Positive_MY +
	Positive_MN + Positive_MM + Close_MY + Close_MN +
	Close_MM + PNIY + PNIN + PNIM + LVIY + LVIN + LVIM +LN2Y + LN2N +
	LN2M + RSY + RSN + RSM + diag_year, data = hnscc, method = "nearest", 	
        distance = "bart", 	ratio = 2, caliper = 0.25,
	distance.options = list(n.samples = 10000, n.burn = 5000, 
	n.thin = 10, seed = 2025))

mdata =  match.data(match, drop.unmatched = FALSE)

wbart0 = 1/mdata$distance[mdata$Study==0]
wbart1 = 1/mdata$distance[mdata$Study==1]

jpeg("WHC.BART.jpeg", width = 7, height = 7, units = 'in', res = 300)
hist(wbart0, main = "Historical Controls", xlab = "Weights")
dev.off()

jpeg("WTP.BART.jpeg", width = 7, height = 7, units = 'in', res = 300)
hist(wbart1, main = "Trial Patients", xlab = "Weights")
dev.off()

#use logistic regression for calculating propensity score

match = matchit(Study ~ Age + White + Black + Asian + Male +
	Smoker + Alcohol + PriorChemo + Larynx + OralCavity + Oropharynx +
	P16P + P16N + P16M + ECSY + ECSN + ECSM + Positive_MY +
	Positive_MN + Positive_MM + Close_MY + Close_MN +
	Close_MM + PNIY + PNIN + PNIM + LVIY + LVIN + LVIM +LN2Y + LN2N +
	LN2M + RSY + RSN + RSM + diag_year, data = hnscc, method = "nearest", 	distance = "glm", ratio = 2, caliper = 0.25)

mdata =  match.data(match, drop.unmatched = FALSE)

wlogit0 = 1/mdata$distance[mdata$Study==0]
wlogit1 = 1/mdata$distance[mdata$Study==1]

jpeg("WHC.LOGIT.jpeg", width = 7, height = 7, units = 'in', res = 300)
hist(wlogit0, main = "Historical Controls", xlab = "Weights")
dev.off()

jpeg("WTP.LOGIT.jpeg", width = 7, height = 7, units = 'in', res = 300)
hist(wlogit1, main = "Trial Patients", xlab = "Weights")
dev.off()

#use random forests for calculating propensity score

match = matchit(Study ~ Age + White + Black + Asian + Male +
	Smoker + Alcohol + PriorChemo + Larynx + OralCavity + Oropharynx +
	P16P + P16N + P16M + ECSY + ECSN + ECSM + Positive_MY +
	Positive_MN + Positive_MM + Close_MY + Close_MN +
	Close_MM + PNIY + PNIN + PNIM + LVIY + LVIN + LVIM +LN2Y + LN2N +
	LN2M + RSY + RSN + RSM + diag_year, data = hnscc, method = "nearest", 	
        distance = "randomforest", ratio = 2, caliper = 0.25)

mdata =  match.data(match, drop.unmatched = FALSE)

wrf0 = 1/mdata$distance[mdata$Study==0]
wrf1 = 1/mdata$distance[mdata$Study==1]

jpeg("WHC.RF.jpeg", width = 7, height = 7, units = 'in', res = 300)
hist(wrf0, main = "Historical Controls", xlab = "Weights")
dev.off()

jpeg("WTP.RF.jpeg", width = 7, height = 7, units = 'in', res = 300)
hist(wrf1, main = "Trial Patients", xlab = "Weights")
dev.off()

