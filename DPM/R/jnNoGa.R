
 ## DPM: Dirichlet Process Mixtures With Low Information Omnibus Priors
 ## Copyright (C) 2019-2026 Prakash Laud and Rodney Sparapani
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

jnNoGa <-
    function(y, states=c(N, rep(0, N-1)), C=rep(1, N),
             phi=NULL,
             neal.m=2, m0=0, k0=0.2, k0.a=1.5, k0.b=7.5, k0.draw=0,
             a0=1.5, b0=0.5, b0.a=0.5, b0.b=1, b0.draw=0,
             alpha=1, alpha.a=0.1, alpha.b=0.1, alpha.draw=0,
             burn=500, keep=2000, thin=10, standard=FALSE,
             N=nrow(y) ## for convenience, not an argument
)
{
    stopifnot(neal.m>=0 & neal.m==as.integer(neal.m))
    
    p = ncol(y) 
    if(length(phi) == 0) { 
        phi <- matrix(c(apply(y, 2, mean), 1/apply(y, 2, var)), 
                      nrow=N+neal.m, ncol=2*p, byrow = TRUE)
        dimnames(phi)[[2]]=c(paste0('mu', 1:p), paste0('tau', 1:p))
    }
    if(length(m0) == 1) m0 <- rep(m0, p) 
    if(length(k0) == 1) k0 <- rep(k0, p)

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
    
    fit=.Call("call_jnNoGa", y, phi, C, states, prior, hyper, mcmc,
              PACKAGE="DPM")

    a=proc.time()[c(1, 3)]-a[c(1, 3)]
    print(c(seconds=a, minutes=a/60))

    ## fit[[1]]$offset=offset
    ## fit[[1]]$scale=scale
    ## fit[[1]]$prior=prior
    ## fit[[1]]$mcmc=mcmc
    ##return(fit)
    
    ret = list()

    ret$C = t(rbind(sapply(fit, function(a) a$C)))
    ret$Cmax = apply(ret$C, 1, max)
    Cmaxmax = max(c(ret$Cmax))
    h = 1:Cmaxmax

    if(Cmaxmax==1) ret$states = matrix(N, ncol=1, nrow=keep)
    else {
        ret$states = t(rbind(sapply(fit,
                                    function(a)
                                        c(a$states,
                                          rep(0,
                                              Cmaxmax-length(a$states))))))
        ret$states = ret$states[ , h]
    }

    ret$phi = array(dim=c(keep, Cmaxmax, 2*p),
                    dimnames=list(NULL, NULL, dimnames(fit[[1]]$phi)[[2]]))
    for(i in 1:keep) {
        k = nrow(fit[[i]]$phi)
        for(j in 1:p) {
        ret$phi[i, , j]  =c(fit[[i]]$phi[ , j],   rep(NA, Cmaxmax-k))[h]
        ret$phi[i, , p+j]=c(fit[[i]]$phi[ , p+j], rep(NA, Cmaxmax-k))[h]
        }
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

    ret$prior=prior
    ret$hyper=hyper
    ret$mcmc=mcmc                     

    return(ret)
}
