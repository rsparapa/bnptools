
library(BART3)

f = function(x)
    -15+20*x[ , 1]+10*x[ , 2]+10*sin(pi*x[ , 4]*x[ , 5])+20*(x[ , 6]-0.5)^3

P = 50
N = 500
set.seed(12)
X = matrix(runif(P*N), nrow=N, ncol=P)
x = t(rmultinom(N, 1, c(0.17, 0.49, 0.34)))
X[ , 1:3] = x
y = rnorm(N, f(X))

summary(y[x[ , 1]==1])
summary(y[x[ , 2]==1])
summary(y[x[ , 3]==1])

miss = rbinom(N, 1, 0.4)
x. = x
x.[which(miss==1), 1:3]=NA
X. = X
X.[ , 1:3] = x.

impute.prob=matrix(0.05, nrow=N, ncol=3)
for(i in 1:N) 
    impute.prob[i, which(x[i, ]==1)]=0.9

post = mc.gbart(X., 1*(y>0), type='pbart', impute.mult=1:3,
             impute.prob=impute.prob, seed=34,
             sparse=TRUE)

post2 = mc.gbart(X., 1*(y>0), type='pbart', sparse=TRUE, seed=22)

post$varprob.mean[1:6]  ## missing imputation
post2$varprob.mean[1:6] ## vs. hot-decking

##cor(post$yhat.train.mean, y)^2
##cor(post2$yhat.train.mean, y)^2

A = 1*(y>0)==0
D = 1*(y>0)==1
sum(outer(post$prob.train.mean[A], post$prob.train.mean[D], "<"))/(sum(A)*sum(D))
sum(outer(post2$prob.train.mean[A], post2$prob.train.mean[D], "<"))/(sum(A)*sum(D))

