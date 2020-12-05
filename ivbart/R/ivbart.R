
## ivbart: Instrumental Variable BART with Dirichlet Process Mixtures
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

ivbart = function(Z, X, T, Y,
   burn=1000, nd=2000, burnf=1000, burnh=1000,
   keepevery=10,
   betabar=0, Abeta=0.01,
   m=200, nc=100,
   power=2, base=0.95,
   k=2, sigmaf=NA, sigmah=NA,
   v=0.17, nu=2.004, a=0.016,
   n = length(Y),
   Imin=1, Imax=floor(0.1*n)+1, psi=0.5, gs=100,
   centermeans=TRUE,
   fs = rep(0, n), hs = rep(0, n),
   mTs=0, mYs=0,
   betas=NA, sTs=NA, sYs=NA, gammas=NA,
## betas=0, mTs=0, mYs=0, sTs=1, sYs=1, gammas=0,
   printevery=100
)
{
### data
if(length(T)!=n) stop('The length of T must equal the length of Y')
if(length(X)==1) X = rep(X, n)
X = cbind(X)
Z = cbind(Z)
ZX = cbind(Z, X)
if(nrow(Z)!=n) stop('The rows of Z must equal the length of Y')
if(nrow(X)!=n) stop('Either X is a constant (no measured confounders)\n or the rows of X must equal the length of Y')

## tau
if(is.na(sigmaf)) {
   tauf = (max(T)-min(T))/(2*k*sqrt(m));
} else {
   tauf = sigmaf/sqrt(m)
}
if(is.na(sigmah)) {
   tauh = (max(Y)-min(Y))/(2*k*sqrt(m));
} else {
   tauh = sigmah/sqrt(m)
}

### prior for alpha
priold = getaprior(length(Y),Imin,Imax,psi,gs)
priag = priold$p
ag = priold$a

### unless provided starting values from TSLS
## first stage manually
tsls1 = lm(T~., data.frame(T = T, Z, X))
## second stage manually
tsls2 = lm(Y~T+., data.frame(Y = Y, T = tsls1$fitted.values, X))
L = t(chol(var(cbind(tsls1$resid, tsls2$resid))))
if(is.na(betas)) betas=tsls2$coeff[2]
if(is.na(sTs)) sTs=L[1, 1]
if(is.na(sYs)) sYs=L[2, 2]
if(is.na(gammas)) gammas=L[2, 1]

res = .Call("cbiv",
   t(ZX),
   t(X),
   T,
   Y,
   burn, nd*keepevery,
   burnf,
   burnh,
   m, nc,
   power, base,
   tauf,
   tauh,
   betabar, Abeta,
   v, nu, a, #base prior
   ag, priag, #alpha prior
   centermeans,
   fs, hs,
   betas, mTs, mYs, sTs, gammas, sYs,
   printevery
)
res$check = NULL
res$ag = ag
res$priag = priag
thin = seq(keepevery, nd, keepevery)
res$dnpart = res$dnpart[thin]
res$dalpha = res$dalpha[thin]
res$dmu1 = res$dmu1[thin, ]
res$dsigma1 = res$dsigma1[thin, ]
res$dmu2 = res$dmu2[thin, ]
res$dsigma2 = res$dsigma2[thin, ]
res$dbeta = res$dbeta[thin]
res$dh = res$dh[thin, ]
res$df = res$df[thin, ]
res$betas = betas
res$sTs = sTs
res$sYs = sYs
res$gammas = gammas
res$tsls1 = tsls1
res$tsls2 = tsls2
return(res)
}
