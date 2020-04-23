
library(BART3)

##options(mc.cores=8)

B <- getOption('mc.cores', 1)
figures = getOption('figures', default='NONE')

data(transplant)

pfit <- survfit(Surv(futime, event) ~ abo, transplant)

# competing risks for type O
plot(pfit[4,], xscale=7, xmax=735, col=1:3, lwd=2, ylim=c(0, 1),
       xlab='t (weeks)', ylab='Aalen-Johansen (AJ) CI(t)')
    legend(450, .4, c("Death", "Transplant", "Withdrawal"), col=1:3, lwd=2)
## plot(pfit[4,], xscale=30.5, xmax=735, col=1:3, lwd=2, ylim=c(0, 1),
##        xlab='t (months)', ylab='Aalen-Johansen (AJ) CI(t)')
##     legend(450, .4, c("Death", "Transplant", "Withdrawal"), col=1:3, lwd=2)

table(transplant$event)
delta <- (as.numeric(transplant$event)-1)
## recode so that delta=1 is cause of interest; delta=2 otherwise
delta[delta==1] <- 4
delta[delta==2] <- 1
delta[delta==4] <- 2
table(delta, transplant$event)

times <- pmax(1, ceiling(transplant$futime/7)) ## weeks
##times <- pmax(1, ceiling(transplant$futime/30.5)) ## months
table(times)

typeO <- 1*(transplant$abo=='O')
typeA <- 1*(transplant$abo=='A')
typeB <- 1*(transplant$abo=='B')
typeAB <- 1*(transplant$abo=='AB')
table(typeA, typeO)

x.train <- cbind(typeO, typeA, typeB, typeAB)

x.test <- cbind(1, 0, 0, 0)
dimnames(x.test)[[2]] <- dimnames(x.train)[[2]]

## run one long MCMC chain in one process
## set.seed(99)
## post <- crisk3.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)

## in the interest of time, consider speeding it up by parallel processing
## run "mc.cores" number of shorter MCMC chains in parallel processes
post <- mc.crisk3.bart(x.train=x.train, times=times, delta=delta,
                       x.test=x.test, seed=99, mc.cores=B)

K <- post$K

typeO.ltx.mean <- apply(post$cif.test, 2, mean)
typeO.ltx.025 <- apply(post$cif.test, 2, quantile, probs=0.025)
typeO.ltx.975 <- apply(post$cif.test, 2, quantile, probs=0.975)

typeO.dth.mean <- apply(post$cif.test2, 2, mean)
typeO.dth.025 <- apply(post$cif.test2, 2, quantile, probs=0.025)
typeO.dth.975 <- apply(post$cif.test2, 2, quantile, probs=0.975)

typeO.wth.mean <- apply(post$cif.test3, 2, mean)
typeO.wth.025 <- apply(post$cif.test3, 2, quantile, probs=0.025)
typeO.wth.975 <- apply(post$cif.test3, 2, quantile, probs=0.975)

plot(pfit[4,], xscale=7, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
       xlab='t (weeks)', ylab='CI(t)', lty=3)
points(c(0, post$times)*7, c(0, typeO.ltx.mean), col=2, type='s', lwd=2)
points(c(0, post$times)*7, c(0, typeO.ltx.025), col=2, type='s', lwd=2, lty=2)
points(c(0, post$times)*7, c(0, typeO.ltx.975), col=2, type='s', lwd=2, lty=2)
points(c(0, post$times)*7, c(0, typeO.dth.mean), col=1, type='s', lwd=2)
points(c(0, post$times)*7, c(0, typeO.dth.025), col=1, type='s', lwd=2, lty=2)
points(c(0, post$times)*7, c(0, typeO.dth.975), col=1, type='s', lwd=2, lty=2)
points(c(0, post$times)*7, c(0, typeO.wth.mean), col=3, type='s', lwd=2)
points(c(0, post$times)*7, c(0, typeO.wth.025), col=3, type='s', lwd=2, lty=2)
points(c(0, post$times)*7, c(0, typeO.wth.975), col=3, type='s', lwd=2, lty=2)
     legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
                       "Death(AJ)", "Withdrawal(AJ)"),
            col=c(2, 2, 1, 3), lwd=2, lty=c(1, 3, 3, 3))
if(figures!='NONE')
    dev.copy2pdf(file=paste(figures, 'liver3-BART.pdf', sep='/'))

##    dev.copy2pdf(file='liver3-BART.pdf')

## checking predict function, but no missing data here
pre <- crisk3.pre.bart(x.train=x.train, times=times, delta=delta,
                       x.test=x.test)

pred <- predict(post, pre$tx.test, mc.cores=B)
K <- pred$K

max(post$cif.test-pred$cif.test)
max(post$cif.test2-pred$cif.test2)
max(post$cif.test3-pred$cif.test3)

typeO.ltx.mean <- apply(pred$cif.test, 2, mean)
typeO.ltx.025 <- apply(pred$cif.test, 2, quantile, probs=0.025)
typeO.ltx.975 <- apply(pred$cif.test, 2, quantile, probs=0.975)

typeO.dth.mean <- apply(pred$cif.test2, 2, mean)
typeO.dth.025 <- apply(pred$cif.test2, 2, quantile, probs=0.025)
typeO.dth.975 <- apply(pred$cif.test2, 2, quantile, probs=0.975)

typeO.wth.mean <- apply(pred$cif.test3, 2, mean)
typeO.wth.025 <- apply(pred$cif.test3, 2, quantile, probs=0.025)
typeO.wth.975 <- apply(pred$cif.test3, 2, quantile, probs=0.975)

plot(pfit[4,], xscale=7, xmax=735, col=1:3, lwd=2, ylim=c(0, 0.8),
       xlab='t (weeks)', ylab='CI(t)', lty=3)
points(c(0, pred$times)*7, c(0, typeO.ltx.mean), col=2, type='s', lwd=2)
points(c(0, pred$times)*7, c(0, typeO.ltx.025), col=2, type='s', lwd=2, lty=2)
points(c(0, pred$times)*7, c(0, typeO.ltx.975), col=2, type='s', lwd=2, lty=2)
points(c(0, pred$times)*7, c(0, typeO.dth.mean), col=1, type='s', lwd=2)
points(c(0, pred$times)*7, c(0, typeO.dth.025), col=1, type='s', lwd=2, lty=2)
points(c(0, pred$times)*7, c(0, typeO.dth.975), col=1, type='s', lwd=2, lty=2)
points(c(0, pred$times)*7, c(0, typeO.wth.mean), col=3, type='s', lwd=2)
points(c(0, pred$times)*7, c(0, typeO.wth.025), col=3, type='s', lwd=2, lty=2)
points(c(0, pred$times)*7, c(0, typeO.wth.975), col=3, type='s', lwd=2, lty=2)
     legend(450, .4, c("Transplant(BART)", "Transplant(AJ)",
                       "Death(AJ)", "Withdrawal(AJ)"),
            col=c(2, 2, 1, 3), lwd=2, lty=c(1, 3, 3, 3))
