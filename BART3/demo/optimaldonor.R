## Example code to demonstrate how to use the BART object to
## produce individual predictions for different donor ages

library(BART3)
stopifnot(require(ggplot2))

##options(mc.cores=8)
(B <- getOption('mc.cores', 1))

# Read in X.test object consisting of 100 pseudo patients
# with characteristics drawn randomly and independently
# from the validation cohort covariate distributions
# This original X.test object includes randomly sampled actual
# donor ages (in decades)
data(optimaldonor)
X.test <- optimaldonor

# Read in BART posterior samples
post <- readRDS(file=system.file("optimaldonor/post.rds", package="BART3"))
# Produce predictions of event probabilities
# Recall higher probabilities indicates worse outcome
# Produce predictions on X.test 
pred = predict(post, X.test, mc.cores=B)  

# Consider a new X.test object with donor age set to 18 (or 1.8 decades)
X.test18 <- X.test
X.test18[ , 33] <- rep(1.8, 100)
# Produce predictions on new X.test corresponding to use of an 18 year
# old donor instead of the original actual donor
# Compute posterior mean prediction of event probability along with 95
# percent intervals
pred18 = predict(post, X.test18, mc.cores=B)  

# Produce differences in predictions between age 18 donor and actual donor
trtdiff <- pred18$prob.test-pred$prob.test
# Compute posterior mean probabilities of better outcome with age 18 donor
trtp <- apply(1*(trtdiff<0), 2, mean)
# Compute posterior mean difference in event probabilities between age 18
# and actual donor, along with 95 percent intervals
trtdiff.mean <- apply(trtdiff, 2, mean)
trtdiff.lower <- apply(trtdiff, 2, quantile, probs=0.025)
trtdiff.upper <- apply(trtdiff, 2, quantile, probs=0.975)

# Plot posterior mean difference between age 18 donor and actual donor
# along with 95 percent posterior intervals
plotdata <- data.frame(pat.index=c(1:100), trtdiff.mean,
                       trtdiff.lower, trtdiff.upper, trtp)

p1 <- ggplot(plotdata, aes(x=pat.index, y=trtdiff.mean)) +
    geom_pointrange(data=plotdata, aes(ymax=plotdata$trtdiff.lower,
                                       ymin=plotdata$trtdiff.upper),
                    colour="black") +
    theme(axis.text=element_text(size=16, face="bold"),
          axis.title=element_text(size=18, face="bold"))
p1

# Plot probability of better outcome with age 18 donor vs. actual donor
p2 <- ggplot(plotdata, aes(x=pat.index, y=trtp)) +
    theme(axis.text=element_text(size=16, face="bold"),
          axis.title=element_text(size=18,face="bold")) +
    geom_point(color="black", size=2) + ylim(0, 1) 
p2

# Obtain population level predictions for age 18 donor outcome and for
# difference between age 18 donor and actual donor at each MCMC sample
pop18 <- apply(pred18$prob.test, 1, mean)
popdiff <- apply(trtdiff, 1, mean)
popdf <- data.frame(popdiff, pop18)

# Histogram of population outcomes if all patients
# had a matched donor of age 18
p3 <- ggplot() +
    theme(axis.text=element_text(size=16, face="bold"),
          axis.title=element_text(size=18, face="bold")) +
    labs(y="Probability Density", x="Risk probability") +
    scale_colour_manual(values=c("black", "black")) +
    geom_density(aes(x=pop18), colour="black", fill="black", size=2,
                 data=popdf, alpha=0.4)
p3

# Histogram of population level differences in outcomes if all patients
# utilized a matched donor of age 18 vs. the actual donor in the dataset
p4 <-  ggplot() + theme(axis.text=element_text(size=16, face="bold"),
                        axis.title=element_text(size=18, face="bold")) +
    labs(y="Probability Density", x="Risk Difference") +
    scale_colour_manual(values=c("black", "black")) +
    geom_density(aes(x=popdiff), colour="black", fill="black", size=2,
                 data=popdf, alpha=0.4)
p4
