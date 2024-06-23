
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2024 Robert McCulloch and Rodney Sparapani
## Rhat.R

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

oldRhat = function(y, C) {
    m = C
    n = length(y)/m
    x = rep(1:m, each=n)
    fit = lm(y~factor(x))
    a = anova(fit)
    B = a$'Mean Sq'[1]
    W = a$'Mean Sq'[2]
    V = ((n-1)*W+B)/n
    Rhat = sqrt(V/W)
    return(list(oldRhat=Rhat))
}

splitRhat = function(y, C) {
    m = C*2
    l <- length(y)
    n = l/m
    x = rep(1:m, each=n)
    fit = lm(y~factor(x))
    a = anova(fit)
    B = a$'Mean Sq'[1]
    W = a$'Mean Sq'[2]
    V = ((n-1)*W+B)/n
    Rhat = sqrt(V/W)

    L <- l-1
    if((L%%2) == 0) L <- L-1
    L <- min(c(99, L))
    for(j in 1:m) {
        i <- which(x == j)
        rho. <- var(y[i])*c(acf(y[i], lag.max=L, plot=FALSE)$acf)[-1]
        if(j == 1) {
            L <- length(rho.)
            if((L%%2) == 0) L <- L-1
            rho <- rho.[1:L]
        } else { rho <- cbind(rho, rho.) }
    }
    rho <- 1-(W-apply(rho, 1, mean))/V
    k <- L%/%2
    P <- 0
    for(i in 1:k) {
        j <- 2*i-1
        P[i] <- rho[j]+rho[j+1]
    }
    if(min(P)<0) k <- max(1, min(which(P<0)))
    P <- P[1:k]
    tau <- -1+2*sum(P)
    if(tau<0) tau <- 1
    ## t = 0
    ## rho = 1
    ## while(t<2 || rho[t]>0 || rho[t-1]>0) {
    ##     t = t+1
    ##     rho[t] = 0
    ##     for(j in 1:m) {
    ##         a = acf(y[x==j], lag.max=t, plot=FALSE)
    ##         rho[t] = rho[t]+a$acf[t+1, 1, 1]/m
    ##     }
    ##     if(t>1 && rho[t]<0 && rho[t-1]<0)
    ##         rho[(t-1):t] = c(0, 0)
    ## }

    ## rho.t = 1 - (W-rho)/V
    ## tau = 1 + 2*sum(rho.t)

    Seff = n*m/tau
    return(list(splitRhat=Rhat, splitSeff=Seff, splitSpct=1/tau,
                B=B, W=W, V=V, rho=rho, P.t = P))
}

maxRhat = function(y, C) {
    l = length(y)
    n = l/C
    if((n %% 2)==1) y=y[-seq(n, l, n)] ## odd length of chain
    m = C*2
    l = length(y)
    n = l/m
    y. = qnorm((rank(y)-0.5)/l)
    a = splitRhat(y., C)
    y. = abs(y-median(y))
    y. = qnorm((rank(y.)-0.5)/l)
    b = splitRhat(y., C)
    maxRhat = max(c(a$splitRhat, b$splitRhat))
    minSeff = min(c(a$splitSeff, b$splitSeff))
    split.rho=a$rho
    folded.rho=b$rho
    k=min(length(split.rho), length(folded.rho))
    ##rho=pmax(split.rho[1:k], folded.rho[1:k])
    ## if(acfPlot) {
    ##     k=length(split.rho)
    ##     plot(0:(k-1), split.rho, ylim=c(-1, 1))
    ##     k=length(folded.rho)
    ##     plot(0:(k-1), folded.rho, ylim=c(-1, 1))
    ## }
    return(list(maxRhat=maxRhat, minSeff=minSeff, minSpct=minSeff/l,
                splitSeff = a$splitSeff, foldedSeff = b$splitSeff,
                splitRhat=a$splitRhat, foldedRhat=b$splitRhat,
                splitrho=split.rho, foldedrho=folded.rho,
                splitP.t = a$P, foldedP.t = b$P)) ## rho=rho)) 
}
