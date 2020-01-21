### AIDS EXAMPLE
library(mxBART)
library(lattice)
library(rpart)
library(rpart.plot)
library(magrittr)
library(lme4)
library(grid)
library(gridExtra)
data(ACTG175)

nd <- 2000 ## number of posterior draws
nb <- 1000 ## number of burn-in draws
ke <- 3    ## keep every third draw
nt <- 100  ## 100 tree
## exclude those who do not have CD4 count at 96 weeks
## or their baseline CD4 is ineligible
## (inclusion criteria: CD4 counts between 200 and 500)
ex <- is.na(ACTG175$cd496) | ACTG175$cd40<200 | ACTG175$cd40>500
train <- as.data.frame(ACTG175)[!ex, -c(10, 14:15, 17, 18, 22:26)]
train <- cbind(1*(ACTG175$strat[!ex]==1), 1*(ACTG175$strat[!ex]==2),
               1*(ACTG175$strat[!ex]==3), train)
dimnames(train)[[2]][1:3] <- paste0('strat', 1:3)
train <- cbind(1*(ACTG175$arms[!ex]==0), 1*(ACTG175$arms[!ex]==1),
               1*(ACTG175$arms[!ex]==2), 1*(ACTG175$arms[!ex]==3), train)
dimnames(train)[[2]][1:4] <- paste0('arm', 0:3)
N <- nrow(train)
train <- train[rep(1:N,each=2),]
cd4 <- (train$cd420-train$cd40)/train$cd40
cd4[seq(2,nrow(train),2)] <- (train$cd496[seq(2,nrow(train),2)]-train$cd40[seq(2,nrow(train),2)])/train$cd40[seq(2,nrow(train),2)]
train$time <- rep(c(20,96),nrow(train)/2)
train$arms <- factor(train$arms)

## Create test data
test <- train[seq(2,nrow(train),2),]
test <- test[rep(1:nrow(test),each=4),]
test$arms <- rep(0:3,nrow(test)/4)
test$arms <- factor(test$arms)
set.seed(21)

frm <- as.formula(paste0('cd4~',paste0('train$',names(train)[-c(1:8,22:23)],collapse='+'),'+train$arms*train$time+(1|train$pidnum)'))
re.sd <- summary(lmer(frm))$sigma
setts <- list(list(prior=1,df=1,scale=re.sd/qt(.95,df=1)))

## Run mixedBART
post <- gbmm(y.train=cd4,
             x.train=train[,-c(1:8,22:23)],
             x.test=test[,-c(1:8,22:23)],
             id.train=list(train$pidnum),
             mxps=setts,
             sparse=FALSE,ntree=nt,
             ndpost=nd,nskip=nb,keepevery=ke)

trt0 <- seq(1,ncol(post$yhat.test),4)
trt1 <- trt0+1
trt2 <- trt0+2
trt3 <- trt0+3

all.diff1 <- cbind(colMeans(post$yhat.test[,trt1]-post$yhat.test[,trt0]),
                  apply(post$yhat.test[,trt1]-post$yhat.test[,trt0],2,quantile,.25),
                  apply(post$yhat.test[,trt1]-post$yhat.test[,trt0],2,quantile,.75),
                  1:1104
                  )
all.diff2 <- cbind(colMeans(post$yhat.test[,trt3]-post$yhat.test[,trt0]),
                   apply(post$yhat.test[,trt3]-post$yhat.test[,trt0],2,quantile,.25),
                   apply(post$yhat.test[,trt3]-post$yhat.test[,trt0],2,quantile,.75),
                   1:1104
                   )
all.diff1 <- all.diff1[order(all.diff1[,1]),]
all.diff2 <- all.diff2[order(all.diff2[,1]),]
scl <- list(tck=c(1,0),
            y=list(at=c(0,.05,.10,.15,.20,.25),
                   labels=c('0.00','0.05','0.10','0.15','0.20','0.25'),
                   cex=.5),
            x=list(at=seq(0,1,.2),
                   labels=paste0(seq(0,100,20),'%'),
                   cex=.5))
## Create waterfall plots for HIV Example, Figure 5.
plt1 <- xyplot(all.diff1[,1]~seq(0,1,length.out=1104),type='l',col='black',lwd=1.5,
               xlim=c(-.01,1.01),
               ylim=c(-.06,.275),xlab=list('Ordered Patients\n(Out of 1,104)',cex=.5),
               scales=scl,
               ylab=list('Predicted Relative \nCD4 Trt. Difference',cex=.5),
               #main=list('(a) Didanosine Monotherapy',cex=.4),
               panel=function(x,y,...){
                   panel.rect(xleft=-1000,xright=2000,ybottom=-50,ytop=30,col='grey95')
                   panel.abline(v=c(.2,.4,.6,.8),col='grey55',lty=3,lwd=.7)
                   panel.abline(h=c(-.05,.05,.1,.15,.2,.25),col='grey55',lty=3,lwd=.7)
                   panel.rect(xleft=x-.001,xright=x+.001,ybottom=all.diff1[,2],ytop=all.diff1[,3],border='grey',col='grey')
                   panel.xyplot(x,y,...)
                   panel.abline(h=0,lty=2,col='black',lwd=1)
               })
plt2 <- xyplot(all.diff2[,1]~seq(0,1,length.out=1104),type='l',col='black',lwd=1.5,
               xlim=c(-.01,1.01),
               ylim=c(-.06,.275),xlab=list('Ordered Patients\n(Out of 1,104)',cex=.5),
               scales=scl,
               ylab=list('Predicted Relative \nCD4 Trt. Difference',cex=.5),
               #main=list('(b) Zudovudine Monotherapy',cex=.4),
               panel=function(x,y,...){
                   panel.rect(xleft=-1000,xright=2000,ybottom=-50,ytop=30,col='grey95')
                   panel.abline(v=c(.2,.4,.6,.8),col='grey55',lty=3,lwd=.7)
                   panel.abline(h=c(-.05,.05,.1,.15,.2,.25),col='grey55',lty=3,lwd=.7)
                   panel.rect(xleft=x-.001,xright=x+.001,ybottom=all.diff2[,2],ytop=all.diff2[,3],border='grey',col='grey')
                   panel.xyplot(x,y,...)
                   panel.abline(h=0,lty=2,col='black',lwd=1)
               })

## TO output or view plot (Figure 5), uncomment this
#png('Fig5.png',height=4.5,width=3,units='in',res=800)
#grid.arrange(grobs=list(plt1,plt2),ncol=1)
#grid.text(c('(a) ddl + AZT Combined vs. AZT Monotherapy','(b) ddl Monotherapy vs. AZT Monotherapy'),x=.9,y=c(.95,.45),just='right',
#          gp=gpar(fontsize=6,fontfamily='Times'))
#dev.off()

all.diff1 <- all.diff1[order(all.diff1[,4]),]
train2 <- train[seq(1,nrow(train),2),]
train2$drugs <- as.factor(train2$drugs)
levels(train2$drugs) <- c('No','Yes')
rrdf <- cbind(out=all.diff1[,1],train2)
colnames(rrdf)[10] <- 'Age'
colnames(rrdf)[14] <- 'IV_User'
colnames(rrdf)[22] <- 'Baseline_CD4'

all.diff2 <- all.diff2[order(all.diff2[,4]),]
rrdf2 <- cbind(out=all.diff2[,1],train2)
colnames(rrdf2)[10] <- 'Age'
colnames(rrdf2)[14] <- 'IV_User'
colnames(rrdf2)[22] <- 'Baseline_CD4'

frm1 <- as.formula(paste0('out~',paste0('',names(rrdf)[-c(1:9,23:26)],collapse='+')))
rr1 <- rpart(frm1,method='anova',control=rpart.control(cp=.1),data=rrdf)
rr2 <- rpart(frm1,method='anova',control=rpart.control(cp=.1),data=rrdf2)
print(rr1)
print(rr2)

## TO output or view plot (Figure 6), uncomment this
## png('Fig6.png',height=4.5,width=3,units='in',res=800)
## layout(matrix(1:2,ncol=1),heights=c(1,1))
## par(family='',font=1,cex=.8)
## prp(rr1,,cex=.8,
##     fallen.leaves = FALSE, type=4, extra=1, varlen=-14, faclen=0, yesno.yshift=-1)
## par(family='Times',font=1,cex=4)
## mtext(side=1,line=4,text='(a) ddl + AZT Combined vs. AZT Monotherapy',cex=.5)
## par(family='',font=1,cex=.8)
## prp(rr2,,cex=.8,
##     fallen.leaves = FALSE, type=4, extra=1, varlen=-14, faclen=0, yesno.yshift=-1)
## par(family='Times',font=1,cex=.4)
## mtext(side=1,line=4,text='(b) ddl Monotherapy vs. AZT Monotherapy',cex=.5)
## dev.off()
