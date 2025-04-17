library(survival)
library(survminer)

data(hnscc, package = 'BART3')

# fit and plot KM curve

fitDFS <- survfit(Surv(DFS_Time, DFS) ~ Study, data = hnscc)

# obtain 2-year DFS rates and difference

twoyr = summary(fitDFS,time=24)
twoyr

diff = twoyr$surv[2] - twoyr$surv[1]
diffSE <- sqrt(twoyr$std.err[2]^2 + twoyr$std.err[1]^2)
diffL = diff - 1.96 *diffSE
diffU = diff + 1.96 *diffSE
c(diff, diffL, diffU)


pdf("km.pdf", width = 7, height = 7)
myplot <- ggsurvplot(fitDFS, hnscc, censor = F, 
                     xlab = "Months from Surgery", cex.xlab = .4, legend.title="",
                     legend.labs=c("Historical Controls", "Trial Patients"), xlim=c(0,36),
                     ylab = "Disease-free Survival Probability", cex.ylab = 1,
                     font.x =  18, font.y = 18, font.tickslab = 16, font.legend = 20,
                     risk.table = F, pval = FALSE, 
                     conf.int = TRUE, conf.int.style = "step", palette = c("red","blue"), 
                     break.time.by = 3)

myplot$plot +  
  geom_segment(aes(x = 24, y = 0, xend = 24, yend = twoyr$surv[2])) + 
  geom_segment(aes(x = 0, y = twoyr$surv[2], xend = 24, yend = twoyr$surv[2])) +
  geom_segment(aes(x = 0, y = twoyr$surv[1], xend = 24, yend = twoyr$surv[1]))
dev.off()





