
library(DPM)

N=20000

## check rcat
p=1:4
rcat(p*0)
rcat(-p)
set.seed(12)
y=integer(N)
for(i in 1:N) y[i]=rcat(p)
table(y)/N
p

## check rtnorm
T=5
set.seed(12)
y=rtnorm(N, T)
x=seq(5, 6, length.out=1001)
plot(x, dnorm(x)/pnorm(5, lower.tail=FALSE), col='blue', type='l')
lines(density(y, from=5))
abline(v=5, h=0)

## check lgamma
set.seed(12)
y=rlgamma(N, 0.01)
summary(y)
Y=exp(y)
table(Y==0)/N
summary(Y)
## underflow at exp(-746)
y=rgamma(N, 0.01)
table(y==0)/N

## check gamma
shape=0.01
set.seed(12)
y=rgamma(N, shape, 1)
table(y==0)/N
set.seed(12)
Y=.rgamma(N, shape, 1)
table(Y==0)/N
quantile(Y)
R=range(Y)
R=c(0, 0.05)
x=seq(R[1], R[2], length.out=1001)

plot(x, dgamma(x, shape, 1), col='blue', type='l')
lines(density(Y, from=0, to=0.05))
lines(density(y, from=0, to=0.05), col='red')
a=density(Y, from=0)

