
library(BART3)

f = function(x)
    5+10*sin(pi*x[ , 1]*x[ , 2]) +  3*x[ , 3]^3

N = 2000
sigma = 1.0 ##y = f(x) + sigma*z where z~N(0, 1)
P = 4       ##number of covariates

set.seed(12)
x.train=matrix(runif(N*P, -2, 2), N, P)
dimnames(x.train)[[2]] <- paste0('x', 1:P)
x.train[ , 4]= 1*(x.train[ , 4]>0)
y.train=(f(x.train)+sigma*rnorm(N))

H=100
x=seq(-2, 2, length.out=H+1)[-(H+1)]
x.test=matrix(0, nrow=H, ncol=P)
x.test[ , 3]=x

B=8
post1=mc.gbart(x.train, y.train, mc.cores=B, seed=12)

post2=mc.gbart(x.train, y.train, mc.cores=B, seed=21, rfinit=TRUE)
write(post2$trees, file='trees.txt')

##pdf('read-forest.pdf')
plot(post2$sigma[ , 1], type='n', ylim=c(0, max(c(post2$sigma))))
for(i in 1:B) {
    lines(post2$sigma[ , i], lty=1)
    lines(post1$sigma[ , i], lty=2)
}
abline(h=sigma, lty=2)
dev.copy2pdf(file='read-forest.pdf')
##dev.off()

