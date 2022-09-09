
options(mc.cores=8)
library(nftbart)
library(lattice)

data(lung)
N=length(lung$status)

##lung$status: 1=censored, 2=dead
##delta: 0=censored, 1=dead
delta=lung$status-1

## this study reports time in days rather than weeks or months
times=lung$time
times=times/7  ## weeks

## matrix of covariates
x.train=cbind(lung[ , -c(1:3)])

joint=tsvs(x.train=x.train, times=times, delta=delta)
dual =tsvs2(xftrain=x.train, xstrain=x.train, times=times, delta=delta)

(names.=dimnames(x.train)[[2]])
df=data.frame(steps=1:20, TSVS=rep(c('both', 'f', 'sd'), each=20))
df=cbind(df, rbind(joint$prob, dual$probf, dual$probs))
str(df)

df.=data.frame(steps=df$steps, TSVS=df$TSVS, prob=df$age, x='age')
for(var in names.[-1])
    eval(parse(
        text=paste0("df.=rbind(df.,",
"data.frame(steps=df$steps, TSVS=df$TSVS, prob=df$", var, ", x='", var, "'))")))
df.$TSVS=factor(df.$TSVS)
df.$x=factor(df.$x)
str(df.)

xyplot(prob~steps|x, df., groups=TSVS, type='l', as.table=TRUE,
       col=c(1, 2, 4), lwd=2,
       panel=function(...){
           if(panel.number()==5) {
               ltext(19, df.$prob[df.$steps==20 & df.$TSVS=='both' &
                                  df.$x=='ph.karno'], 'both', pos=1)
               ltext(19, df.$prob[df.$steps==20 & df.$TSVS=='f' &
                                  df.$x=='ph.karno'], 'f', col=2, pos=3)
               ltext(19, df.$prob[df.$steps==20 & df.$TSVS=='sd' &
                                  df.$x=='ph.karno'], 'sd', col=4, pos=3)
           }
           panel.xyplot(...)
       })
