
 ## DPM: Dirichlet Process Mixtures With Low Information Omnibus Priors
 ## Copyright (C) 2019 Prakash Laud and Rodney Sparapani
 ##
 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 2 of the License, or
 ## (at your option) any later version.
 ##
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 ##
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, a copy is available at
 ## https://www.R-project.org/Licenses/GPL-2

mutau <-
    function(y, delta=rep(1, N), states=c(N, rep(0, N-1)), C=rep(1, N),
             phi=matrix(c(mean(y), 1/var(y)),
                        nrow=N+neal.m, ncol=2, byrow=TRUE),
             neal.m=2, m0=0, k0=0.2, k0.a=1.5, k0.b=7.5, k0.draw=1,
             a0=1.5, b0=0.5, b0.a=0.5, b0.b=1, b0.draw=1,
             alpha=1, alpha.a=0.1, alpha.b=0.1, alpha.draw=1,
             burn=500, keep=2000, thin=10, standard=FALSE,
             N=length(y) ## for convenience, not an argument
             )
{

    stopifnot(neal.m>=0 & neal.m==as.integer(neal.m))
    
    if(sum(delta)==N) {
        Q=quantile(y, probs=(1:3)/4)
        ##Q=quantile(y, probs=c((1:3)/4, 0.95))
        attr(Q, 'names')=NULL
        scale=Q[3]-Q[1]
        ##scale=(Q[4]-Q[2])/2
        offset=Q[2]/scale
        Y=cbind(y)
    }
    else {
        fit <- survreg(Surv(time=exp(y), event=delta)~1)
        scale=fit$scale
        attr(fit$coefficients, 'names')=NULL
        offset=(fit$coefficients/scale)-0.5772156649015328655494
        ## Extreme Value dist correction with the Euler-Mascheroni constant
        Y=cbind(y, y, delta)
    }

    if(!standard) {
        scale2=scale^2
        m0=m0+offset
        a0=a0+0.5
        k0.b=k0.b/scale2
        b0.b=b0.b/scale2
    }

    dimnames(phi)[[2]]=c('mu', 'tau')
    prior=list(m=as.integer(neal.m),
               m0=m0, k0.a=k0.a, k0.b=k0.b,
               a0=a0, b0.a=b0.a, b0.b=b0.b,
               alpha.a=alpha.a, alpha.b=alpha.b)
    hyper=list(alpha=alpha, alpha.draw=alpha.draw,
               k0=k0, b0=b0, k0.draw=k0.draw, b0.draw=b0.draw)
    C=as.integer(C-1) ## convert from R to C/C++ indexing
    states=as.integer(states)
    mcmc=list(burn=as.integer(burn), keep=as.integer(keep),
              thin=as.integer(thin))
    a=proc.time()

    ## return(list(Y=Y, phi=phi, C=C, states=states, prior=prior,
    ##             hyper=hyper, mcmc=mcmc))
    
    fit=.Call("call_mutau", Y, phi, C, states, prior, hyper, mcmc,
              PACKAGE="DPM")

    a=proc.time()[c(1, 3)]-a[c(1, 3)]
    print(c(seconds=a, minutes=a/60))

    ## fit[[1]]$offset=offset
    ## fit[[1]]$scale=scale
    ## fit[[1]]$prior=prior
    ## fit[[1]]$mcmc=mcmc
    ## return(fit)
    
    ret = list()

    ret$C = t(rbind(sapply(fit, function(a) a$C)))
    ret$Cmax = apply(ret$C, 1, max)
    Cmaxmax = max(c(ret$Cmax))
    ret$states = t(rbind(sapply(fit, function(a) a$states)))
    check = 1:Cmaxmax
    ##check = which(apply(ret$states, 2, sum)>0)
    ret$states = ret$states[ , check]

    ret$phi = array(dim=c(keep, Cmaxmax, 2),
                    dimnames=list(NULL, NULL, dimnames(fit[[1]]$phi)[[2]]))
    for(i in 1:keep) {
        ret$phi[i, , 1]=fit[[i]]$phi[check, 1]
        ret$phi[i, , 2]=fit[[i]]$phi[check, 2]
    }
    
    if(alpha.draw==1)
        ret$alpha=sapply(fit, function(a) a$hyper$alpha)
    else ret$alpha = alpha

    if(k0.draw==1)
        ret$k0=sapply(fit, function(a) a$hyper$k0)
    else ret$k0 = k0

    if(b0.draw==1)
        ret$b0=sapply(fit, function(a) a$hyper$b0)
    else ret$b0 = b0

    ret$offset=offset
    ret$scale=scale
    ret$prior=prior
    ret$mcmc=mcmc                     

    return(ret)
}
