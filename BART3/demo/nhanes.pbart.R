
options(mc.cores=8)
library(BART3)

B <- getOption('mc.cores', 1)
## figures = getOption('figures', default='NONE')

data(arq)
str(arq)

(N <- length(arq$riagendr))
table(arq$riagendr, useNA='ifany')
quantile(arq$bmxwaist, (1:19)/20)
pick=c(5:6, 10)
x.train=arq[ , pick]
print(cor(x.train))

post1 <- mc.gbart(x.train=x.train, y.train=arq$arq010de, ## lower
                  seed=21, type='pbart')

post2 <- mc.gbart(x.train=x.train, y.train=arq$arq010a, ## neck
                  seed=11, type='pbart')

(bmxwaist <- seq(80, 120, by=5))
H <- length(bmxwaist)
x.test = x.train[1:H, ]
x.test[ , 'bmxwaist']=bmxwaist
x.test=rbind(x.test, x.test)
x.test[ , 'riagendr']=rep(1:2, each=H)
print(x.test)

file.=c('nhanes-pbart-FPD.pdf', 'nhanes-pbart-FPDK.pdf') 
for(i in 1:2) {
    if(i==1) {
        pred1=FPD(post1, x.test, S=c(1, 3))
        pred2=FPD(post2, x.test, S=c(1, 3))
    } else {
        pred1=FPDK(post1, x.test, S=c(1, 3), mult.impute=20)
        pred2=FPDK(post2, x.test, S=c(1, 3), mult.impute=20)
    }
    pdf(file=file.[i])
    par(mfrow=c(1, 2))
    plot(bmxwaist, pred1$prob.test.mean[1:H], type='l', col='blue',
         ylim=c(0, 0.35), xlab='Waist circumference (cm)', ylab=expression(p(x)),
         sub='Low-back pain: M(blue) vs. F(red)')
    lines(bmxwaist, pred1$prob.test.lower[1:H], type='l', col='blue', lty=2)
    lines(bmxwaist, pred1$prob.test.upper[1:H], type='l', col='blue', lty=2)
    lines(bmxwaist, pred1$prob.test.mean[H+1:H], type='l', col='red')
    lines(bmxwaist, pred1$prob.test.lower[H+1:H], type='l', col='red', lty=2)
    lines(bmxwaist, pred1$prob.test.upper[H+1:H], type='l', col='red', lty=2)
    abline(h=0:1)
    lines(bmxwaist, rep(mean(pred1$prob.test.mean[1:H]), H), type='l', lty=3, col='blue')
    lines(bmxwaist, rep(mean(pred1$prob.test.mean[H+1:H]), H), type='l', lty=3, col='red')
    plot(bmxwaist, pred2$prob.test.mean[1:H], type='l', col='blue',
         ylim=c(0, 0.35), xlab='Waist circumference (cm)', ylab='', ##expression(p(x)),
         sub='Neck pain: M(blue) vs. F(red)')
    lines(bmxwaist, pred2$prob.test.lower[1:H], type='l', col='blue', lty=2)
    lines(bmxwaist, pred2$prob.test.upper[1:H], type='l', col='blue', lty=2)
    lines(bmxwaist, pred2$prob.test.mean[H+1:H], type='l', col='red')
    lines(bmxwaist, pred2$prob.test.lower[H+1:H], type='l', col='red', lty=2)
    lines(bmxwaist, pred2$prob.test.upper[H+1:H], type='l', col='red', lty=2)
    abline(h=0:1)
    lines(bmxwaist, rep(mean(pred2$prob.test.mean[1:H]), H), type='l', lty=3, col='blue')
    lines(bmxwaist, rep(mean(pred2$prob.test.mean[H+1:H]), H), type='l', lty=3, col='red')
    par(mfrow=c(1, 1))
    dev.off()
}
