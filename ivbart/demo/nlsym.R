
library(ivbart)
data(nlsym)

Y=nlsym$lwage76
subset=!is.na(Y)
table(subset)
Y=Y[subset]
T=nlsym$ed76[subset]
Z=cbind(nlsym$nearc2, nlsym$nearc4)[subset, ]
dimnames(Z)[[2]]=c('nearc2', 'nearc4')
X=cbind(nlsym$exp76, (nlsym$exp76-mean(nlsym$exp76))^2,
        nlsym$black, nlsym$smsa76r, nlsym$reg76r)[subset, ]
dimnames(X)[[2]]=c('exp76', 'exp762', 'black',
                   'smsa76r', 'reg76r')

sd.Y = sd(Y)
sd.T = sd(T)
T. = (T-mean(T))/sd.T ## ending in dot for rescaled
Y. = (Y-mean(Y))/sd.Y
T = T - mean(T)
Y = Y - mean(Y)
ratio = sd.Y/sd.T

ols = lm(Y~T+., data.frame(Y=Y, T=T, X)) ## original scale
print(summary(ols))
beta.ols = ols$coef[2]

ols. = lm(Y.~T.+., data.frame(Y.=Y., T.=T., X)) ## rescaled
print(summary(ols.))

cat('check:', beta.ols, ols.$coef[2]*ratio,
    'these should be the same\n')

## TSLS
## first stage manually
tsls1 = lm(T.~., data.frame(T.=T., Z, X))
## second stage manually
tsls2 = lm(Y.~T.+.,
           data.frame(Y.=Y., T. = tsls1$fitted.values, X))
print(summary(tsls2))

## and if sem installed
if(require(sem)) {
    tsls. = tsls(Y.~T.+exp76+exp762+black+smsa76r+reg76r,
                   ~nearc2+nearc4+exp76+exp762+black+smsa76r+reg76r,
                 data.frame(Y.=Y., T.=T., Z, X))
    print(summary(tsls.))
    cat("check:", tsls2$coef[2], tsls.$coeff[2],
        'these should be the same\n')
} else { tsls. = tsls2 }

(beta.tsls = tsls.$coeff[2])   ## rescaled
(beta.tsls = beta.tsls * ratio)## original scale

cat("beta on original scale: OLS=", beta.ols,
    "& TSLS=", beta.tsls, "\n")

## IVBART
L = t(chol(var(cbind(tsls1$resid, tsls.$resid))))
X. = X[ , -2] ## without the square BART will not need
set.seed(99)
post = ivbart(Z, X., T., Y.,
              nd=1000, sigmaf=1, sigmah=1,
              betas=tsls.$coeff[2], ## beta for rescaled data
              sTs=L[1, 1], gammas=L[2, 1], sYs=L[2, 2],
              Imin=2, Imax=floor(0.5*length(Y))+1, gs=500)

beta.ivbart = post$dbeta * ratio

max.=1.2*exp(max(beta.ols, beta.tsls, beta.ivbart))
plot(density(exp(beta.ivbart)),
     xlim=c(1/max., max.), ylim=c(0, 50),
     log='x', lwd=2, main=expression(beta),
     xlab='Additional year of schooling: multiple of wages')
abline(h=0)
abline(v=1, lty=2)
abline(v=exp(beta.ols), lwd=2, col="blue")
abline(v=exp(beta.tsls), lwd=2, col="red")
legend('topright', lwd=2, legend=c('IVBART', 'OLS', 'TSLS'),
       lty=1, col=c('black', 'blue', 'red'))

