## Copyright (C) 2021-2022 Rodney A. Sparapani

## This file is part of nftbart.
## nft.R

## nftbart is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.

## nftbart is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Author contact information
## Rodney A. Sparapani: rsparapa@mcw.edu

nft = function(## data
               x.train, times, delta=NULL, 
               x.test=matrix(nrow=0, ncol=0),
               impute.bin=NULL, impute.prob=NULL,
               ## multi-threading
               tc=1, ##OpenMP thread count
               ##MCMC
               nskip=1000, ndpost=2000, 
               nadapt=1000, adaptevery=100, 
               chv = cor(x.train, method="spearman"),
               pbd=c(0.7, 0.7), pb=c(0.5, 0.5),
               stepwpert=c(0.1, 0.1), probchv=c(0.1, 0.1),
               minnumbot=c(5, 5),
               ## BART and HBART prior parameters
               ntree=c(50, 10), numcut=100, xicuts=NULL,
               power=c(2, 2), base=c(0.95, 0.95),
               ## f function
               k=5, sigmaf=NA,
               dist='weibull', 
               ## s function
               sigmav=NULL, total.lambda=NA, total.nu=10, mask=NULL,
               ## survival analysis 
               K=100, events=NULL, 
               ##impute.mult=NULL, impute.prob=NULL, impute.miss=NULL,
               ## DPM LIO
               drawDPM=1L, 
               alpha=1, alpha.a=1, alpha.b=0.1, alpha.draw=1,
               neal.m=2, constrain=1, 
               m0=0, k0.a=1.5, k0.b=7.5, k0=1, k0.draw=1,
               ##a0=1.5, b0.a=0.5,
               a0=3, b0.a=2, b0.b=1, b0=1, b0.draw=1,
               ## misc
               na.rm=FALSE, probs=c(0.025, 0.975), printevery=100
               )
{
    n=length(times)
    if(length(delta)==0) delta=rep(1, n)
    if(length(sigmav)==0) sigmav=rep(1, n)
    x = x.train
    xp = x.test
    ##if(!transposed) {
        impute=CDimpute(x.train=cbind(x),
                        x.test=cbind(xp),
                        impute.bin=impute.bin)
        x=t(impute$x.train)
        xp=t(impute$x.test)
        transposed=TRUE
    ##}

    take.logs = (dist!='extreme')
    if(length(events)==0 && K>0) {
        events <- unique(quantile(times, probs=(0:(K-1))/(K-1)))
        attr(events, 'names') <- NULL
    }
    K <- length(events)

    burn=nskip
    nd=ndpost
    nc=numcut
    
    p=nrow(x)
    np=ncol(xp)

    if(length(xicuts)==0) {
        xicuts=xicuts(x, transposed=transposed, numcut=numcut)
        ##return(xicuts)
    }

    df = data.frame(times, delta, t(x))
    ##type='interval' not supported by survreg
    aft = survreg(Surv(time=times, event=delta)~1,
                  data=subset(df, delta<2), dist=dist)

    fmu=coef(aft)[1]
    attr(fmu, 'names')=NULL
    if(take.logs) {
        y=log(times)-fmu
        if(K>0) events=log(events)
    } else {
        ##  fmu=mean(times)
        y=times-fmu
    }

    if(is.na(total.lambda)) total.lambda = (aft$scale^2)
        ##total.lambda = (aft$scale^2)*qchisq(1-total.q, total.nu)/total.nu

    ## tau
    if(is.na(sigmaf)) {
        tau=(max(y)-min(y))/(2*k*sqrt(ntree[1]))
    } else {
        tau = sigmaf/sqrt(ntree[1])
    }

    if(drawDPM==1) {
        C=rep(1, n)
        C=as.integer(C-1) ## C/C++ indexing
        states=as.integer(c(n, rep(0, n-1)))
            prior=list(m=as.integer(neal.m), constrain=as.integer(constrain),
                       m0=m0, k0.a=k0.a, k0.b=k0.b,
                       a0=a0, b0.a=b0.a, b0.b=b0.b,
                       ##a0=3, b0.a=2, b0.b=1,
                       ##a0=1.5, b0.a=0.5, b0.b=1,
                       alpha.a=alpha.a, alpha.b=alpha.b)
        phi=matrix(c(0, prior$b0.b), nrow=n, ncol=2, byrow=TRUE)
        hyper=list(alpha=alpha, alpha.draw=alpha.draw,
               k0=k0, b0=b0, k0.draw=k0.draw, b0.draw=b0.draw)
    }

    check = length(impute.bin)
    if(!(check %in% 0:1))
        stop("The number of imputed binomial columns must be 0 or 1")
    if(check==1) {
        impute.flag=TRUE
        if(length(impute.prob)==0) {
            impute.prob=x[impute.bin, ] ## transposed
            impute.prob[is.na(impute.prob)]=
                mean(x[impute.bin, ]==1, na.rm = TRUE) ## transposed
        }
        impute.bin=as.integer(impute.bin-1) ## convert from R index to C/C++
    } else {
        impute.flag=FALSE
        impute.bin=-1
        impute.prob=0
    }

    ## check = length(impute.mult)
    ## if(check==1)
    ##     stop("The number of multinomial columns must be greater than 1\nConvert a binary into two columns")
    ## if(check>1) {
    ##     impute.flag=TRUE
    ##     if(length(impute.miss)==0) {
    ##         impute.miss=integer(n)
    ##         for(j in 1:check) {
    ##             i=impute.mult[j]
    ##             impute.miss = pmax(impute.miss, is.na(x[i, ])) ## transposed
    ##         }
    ##     }
    ##     if(length(impute.prob)==0) {
    ##         impute.prob=double(check)
    ##         for(j in 1:check) {
    ##             i=impute.mult[j]
    ##             impute.prob[j]=sum(x[i, ]==1, na.rm = TRUE) ## transposed
    ##         }
    ##         impute.prob=impute.prob/sum(impute.prob)
    ##         impute.prob=matrix(impute.prob,
    ##                            nrow=n, ncol=check, byrow=TRUE)
    ##     }
    ##     impute.mult=as.integer(impute.mult-1) ## convert from R index to C/C++
    ## } else {
    ##     impute.flag=FALSE
    ##     impute.mult=integer(0) ## integer vector of column indicators for missing covariates
    ##     impute.miss=integer(0) ## integer vector of row indicators for missing values
    ##     impute.prob=matrix(nrow=0, ncol=0)
    ## }

    ptm <- proc.time()

    res=.Call("cnft",
              x,
              y,
              as.integer(delta),
              xp,
              c(ntree[1], ntree[2]),
              nd,
              burn,
              nadapt,
              adaptevery,
              tau,
              total.lambda,
              total.nu,
              base,
              power,
              tc,
              sigmav,
              chv,
              pbd,
              pb,
              stepwpert,
              probchv,
              minnumbot,
              printevery,
              xicuts,
              ##summarystats,
              ##alphao, betao,
              ##mstart, sstart,
              hyper, C, states, phi, prior,
              ##draws,
              drawDPM,
              impute.bin,
              impute.prob,
              ## impute.mult,
              ## impute.miss,
              ## impute.prob,
              PACKAGE="nftbart")

    ## res$elapsed <- (proc.time()-ptm)['elapsed']
    ## attr(res$elapsed, 'names')=NULL
    res$K=K

if(K>0) {
    if(take.logs) res$events=exp(events)
    else res$events=events
    ## if(take.logs) res$times=exp(events)
    ## else res$times=events
    res$take.logs=take.logs
}
    
    s.train.max.=NULL
    mask.=length(mask)
    if(mask.==0) s.train.mask=1:nd
    else if(mask.==1) {
        s.train.max =apply(res$s.train, 1, quantile, probs=1)
        s.train.max.=quantile(s.train.max, probs=mask)
        s.train.mask=(s.train.max<=s.train.max.)
        nd=sum(s.train.mask)
        s.train.mask=which(s.train.mask)
        res$f.train=res$f.train[s.train.mask, ]
        res$s.train=res$s.train[s.train.mask, ]
        res$z.train=res$z.train[s.train.mask, ]
    } else stop('mask should be of length 0 or 1')
    s.train.burn=c(1:burn, s.train.mask)
   
    if(drawDPM>0) {
        res$dpalpha=res$dpalpha[s.train.mask]
        res$dpn=res$dpn[s.train.burn]
        res$dpmu=res$dpmu[s.train.mask, ]
        res$dpsd=res$dpsd[s.train.mask, ]
        res$dpC=res$dpC[s.train.mask, ]
        if(burn>0) res$dpn.=res$dpn[-(1:burn)]
        else res$dpn.=res$dpn
    }

    ##if(take.logs) res$z.train=exp(res$z.train)
    res$f.train=res$f.train+fmu
    res$f.train.mean=apply(res$f.train, 2, mean)
    res$z.train=res$z.train+fmu
    res$z.train.mean=apply(res$z.train, 2, mean)
    res$e.train.mean=res$z.train.mean-res$f.train.mean
    res$MSE=mean(res$e.train.mean^2)
        res$s.train.mean=apply(res$s.train, 2, mean)
        res$ssd=sqrt(mean(res$s.train.mean^2))
    res$mask=mask
    res$s.train.mask=s.train.mask
    res$s.train.max=s.train.max.
        
    if(np>0) {
        res$f.test=res$f.test+fmu
        res$f.test.mean=apply(res$f.test, 2, mean)
        res$s.test.mean=apply(res$s.test, 2, mean)
        res$x.test=t(xp)
    }

    res$times=times
    res$delta=delta
    res$x.train=t(x)
    res$ntree=ntree[1]
    res$ntree[2]=ntree[2]
    res$ndpost=ndpost
    res$xicuts=xicuts
    res$fmu = fmu

    summarystats=(nd>=2)
    if(summarystats) {
        dimnames(res$f.varcount)[[2]]=names(xicuts)
        res$f.varcount[2:ndpost, ]=res$f.varcount[2:ndpost, ]-res$f.varcount[1:(ndpost-1), ]
        res$f.varcount=res$f.varcount[s.train.mask, ]
        res$f.varcount.mean=apply(res$f.varcount, 2, mean)
        res$f.varprob=res$f.varcount.mean/sum(res$f.varcount.mean)
            dimnames(res$s.varcount)[[2]]=names(xicuts)
            res$s.varcount[2:ndpost, ]=res$s.varcount[2:ndpost, ]-res$s.varcount[1:(ndpost-1), ]
        res$s.varcount=res$s.varcount[s.train.mask, ]
            res$s.varcount.mean=apply(res$s.varcount, 2, mean)
            res$s.varprob=res$s.varcount.mean/sum(res$s.varcount.mean)
    }
    
    if(drawDPM>0) {
        res$prior = prior
        res$hyper = hyper
        res$C = as.integer(C+1)
        res$states = states
        res$phi = phi

        H=max(res$dpn.)
        dpwt. = matrix(0, nrow=nd, ncol=H)
        for(h in 1:H) dpwt.[ , h]=apply(res$dpC==h, 1, sum)
        res$dpwt.=dpwt./n

        if(H<n) { ## it has to be, but just in case 
            res$dpmu.=res$dpmu.[ , -((H+1):n)]
            res$dpsd.=res$dpsd.[ , -((H+1):n)]
            for(i in 1:nd) {
                h=res$dpn.[i]
                if(h<H) for(j in (h+1):H) {
                            res$dpmu.[i, j]=0
                            res$dpsd.[i, j]=1
                        }
            }
        }
        
    }
    attr(res, 'class') <- 'nft'
    
    mu. = res$dpmu*res$s.train+res$f.train
    sd. = res$dpsd*res$s.train
    ## cpo = 1/apply(1/dnorm(res$z.train, mu., sd.), 2, mean)
    ## res$lpml = sum(log(cpo))
    
    z=matrix(log(times), nrow=nd, ncol=n, byrow=TRUE)
    delta=matrix(delta, nrow=nd, ncol=n, byrow=TRUE)
    cpo=(delta==2)*pnorm(z, mu., sd.)+(delta==0)*pnorm(z, mu., sd., FALSE)+
        (delta==1)*dnorm(z, mu., sd.)
    cpo = 1/apply(1/cpo, 2, mean)
    res$LPML = sum(log(cpo))

    res$pred=predict(res, res$x.train, tc=tc, XPtr=FALSE,
                     soffset=0, probs=probs, na.rm=na.rm) 
    res$soffset=0.5*log(mean(res$pred$s.test.mean^2)/
                        mean(res$s.train.mean^2)) ## for stability
    
    res$drawDPM=drawDPM
    res$aft=aft
    ##res$total.lambda=total.lambda
    res$elapsed <- (proc.time()-ptm)['elapsed']
    attr(res$elapsed, 'names')=NULL
    return(res)
}
