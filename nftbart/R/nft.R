## Copyright (C) 2021 Rodney A. Sparapani

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
               sigmav=NULL, overalllambda=NA, overallnu=10,
               ## survival analysis 
               K=100, events=NULL, 
               ##impute.mult=NULL, impute.prob=NULL, impute.miss=NULL,
               ## DPM LIO
               drawMuTau=1L, 
               alpha=1, alpha.a=1, alpha.b=0.1, alpha.draw=1,
               neal.m=2, constrain=1, 
               ##m0=0, k0=0, k0.a=1.5, k0.b=7.5,
               ##a0=1.5, b0=0, b0.a=0.5, b0.b=1,
               k0.draw=1, b0.draw=1,
               ## misc
               printevery=100
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
    ## if(length(xicuts)==0) {
    ##     xicuts=list()
    ##     minx=apply(x, 1, min, na.rm=TRUE)
    ##     maxx=apply(x, 1, max, na.rm=TRUE)
    ##     for(i in 1:p) {
    ##         xinc=(maxx[i]-minx[i])/(numcut+1)
    ##         xicuts[[i]]=(1:numcut)*xinc+minx[i]
    ##     }
    ##     class(xicuts)="BARTcutinfo"
    ## }

    df = data.frame(times, delta, t(x))
    ##type='interval' not supported by survreg
    lm1 = survreg(Surv(time=times, event=delta)~1,
                  data=subset(df, delta<2), dist=dist)

    ##nuo = alphao*2

    ##if(draws) {
        ##betao=alphao
    ## } else {
    ##     if(p > n) {
    ##         ## we use stepwise selection to reduce p<n
    ##         step=srstepwise(t(x)[delta<2, ], times[delta<2], delta[delta<2], dist=dist)
    ##         if(length(step)==0) lm2 = lm1
    ##         else lm2 = survreg(Surv(time=times[delta<2], event=delta[delta<2])~
    ##                                I(t(x)[delta<2, step]), dist=dist)
    ##     }
    ##     else lm2 = survreg(Surv(time=times, event=delta)~.,
    ##                        data=subset(df, delta<2), dist=dist)

    ##     ##sigquanto=.95
    ##     sigesto = lm2$scale
    ##     qchio = qchisq(1-sigquanto, nuo)
    ##     lambdao = (sigesto^2)*qchio/nuo
    ##     betao = nuo*lambdao/2
    ## }

    fmu=coef(lm1)[1]
    attr(fmu, 'names')=NULL
    if(take.logs) {
        y=log(times)-fmu
        if(K>0) events=log(events)
    } else {
        ##  fmu=mean(times)
        y=times-fmu
    }

    if(is.na(overalllambda)) overalllambda = var(y)

    ## if(is.na(tauinit)) { #get tauinit from sstart
    ##     sstart = lm2$scale
    ##     tauinit = 1/sstart^2
    ## } else { #get sstart from tauinit
    ##     sstart = 1/sqrt(tauinit)
    ## }
    ## mstart=mean(y)

    ## tau
    if(is.na(sigmaf)) {
        tau=(max(y)-min(y))/(2*k*sqrt(ntree[1]))
    } else {
        tau = sigmaf/sqrt(ntree[1])
    }

    ## LIO.a=1/lm1$scale
    ## LIO.b=exp(coef(lm1)[1])
    ## LIO.Q2=log(qweibull(0.5, LIO.a, LIO.b))
    ## LIO.Q95=log(qweibull(0.95, LIO.a, LIO.b))
    ## LIO.scale=0.5*(LIO.Q95-LIO.Q2)
    ## LIO.tau=1/(LIO.scale^2)
    ## if(drawMuTau==2) {
    ##     prior=list(m=as.integer(neal.m), constrain=as.integer(constrain),
    ##                a0=a0, b0=b0,
    ##                alpha.a=alpha.a, alpha.b=alpha.b)
    ##     phi=matrix(c(0, 1), nrow=n, ncol=2, byrow=TRUE)
    ##     hyper=list(alpha=alpha, alpha.draw=alpha.draw)
    ## } else 
    if(drawMuTau==1) {
        C=rep(1, n)
        C=as.integer(C-1) ## C/C++ indexing
        states=as.integer(c(n, rep(0, n-1)))
            prior=list(m=as.integer(neal.m), constrain=as.integer(constrain),
                       m0=0, k0.a=1.5, k0.b=7.5,
                       a0=3, b0.a=2, b0.b=1,
                       ##a0=1.5, b0.a=0.5, b0.b=1,
                       alpha.a=alpha.a, alpha.b=alpha.b)
        ## else ##unstandardized
        ##     prior=list(m=as.integer(neal.m), constrain=as.integer(constrain),
        ##                ##m0=LIO.Q2/LIO.scale,
        ##                m0=0, k0.a=1.5, k0.b=7.5*LIO.tau,
        ##                a0=2, b0.a=0.5, b0.b=LIO.tau,
        ##                alpha.a=alpha.a, alpha.b=alpha.b)
        phi=matrix(c(0, prior$b0.b), nrow=n, ncol=2, byrow=TRUE)
        ##phi=matrix(c(LIO.Q2/LIO.scale, LIO.tau), nrow=n, ncol=2, byrow=TRUE)
        hyper=list(alpha=alpha, alpha.draw=alpha.draw,
               k0=0, b0=0, k0.draw=k0.draw, b0.draw=b0.draw)
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
              overalllambda,
              overallnu,
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
              drawMuTau,
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
    
    if(drawMuTau>0) {
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
    ##if(draws) {
        res$s.train.mean=apply(res$s.train, 2, mean)
        res$ssd=sqrt(mean(res$s.train.mean^2))
    ## }
    ## else if(drawMuTau==0) {
    ##     if(burn>0) sigma=res$sigma[-(1:burn)]
    ##     else sigma=res$sigma
    ## }

    if(np>0) {
        res$f.test=res$f.test+fmu
        ##res$f.test.mask=is.finite(apply(abs(res$f.test), 1, max))
            res$s.test.mean=apply(res$s.test, 2, mean)
            ##res$s.test.mask=is.finite(apply(res$s.test, 1, max))
            ##res$f.test.mask=(res$f.test.mask&res$s.test.mask)
            res$s.test.mean =apply(res$s.test, 2, mean)
            ##res$s.test.mean =apply(res$s.test[res$f.test.mask, ], 2, mean)
        res$f.test.mean =apply(res$f.test, 2, mean)
        ##res$f.test.mean =apply(res$f.test[res$f.test.mask, ], 2, mean)

        res$x.test=t(xp)
    }

    res$times=times
    res$delta=delta
    res$x.train=t(x)
    res$ntree=ntree[1]
    res$ntree[2]=ntree[2]
    ## if(draws) res$ntree[2]=ntree[2]
    ## else res$scale = lm2$scale
    res$ndpost=nd
    res$xicuts=xicuts
    res$fmu = fmu

    summarystats=(nd>=2)
    if(summarystats) {
        dimnames(res$f.varcount)[[2]]=names(xicuts)
        res$f.varcount[2:ndpost, ]=res$f.varcount[2:ndpost, ]-res$f.varcount[1:(ndpost-1), ]
        res$f.varcount.mean=apply(res$f.varcount, 2, mean)
        res$f.varprob=res$f.varcount.mean/sum(res$f.varcount.mean)
            dimnames(res$s.varcount)[[2]]=names(xicuts)
            res$s.varcount[2:ndpost, ]=res$s.varcount[2:ndpost, ]-res$s.varcount[1:(ndpost-1), ]
            res$s.varcount.mean=apply(res$s.varcount, 2, mean)
            res$s.varprob=res$s.varcount.mean/sum(res$s.varcount.mean)
    }
    
    if(drawMuTau>0) {
        ##res$alphao = alphao
        ##res$betao = betao
        res$prior = prior
        res$hyper = hyper
        ##res$sstart = sstart
        res$C = as.integer(C+1)
        res$states = states
        res$phi = phi

        H=max(res$dpn.)
        dpwt. = matrix(0, nrow=ndpost, ncol=H)
        for(h in 1:H) dpwt.[ , h]=apply(res$dpC==h, 1, sum)
        res$dpwt.=dpwt./n

        if(H<n) { ## it has to be, but just in case 
            res$dpmu.=res$dpmu.[ , -((H+1):n)]
            res$dpsd.=res$dpsd.[ , -((H+1):n)]
            for(i in 1:ndpost) {
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
    cpo = 1/apply(1/dnorm(res$z.train, mu., sd.), 2, mean)
    res$lpml = sum(log(cpo))
    
    z=matrix(log(times), nrow=ndpost, ncol=n, byrow=TRUE)
    delta=matrix(delta, nrow=ndpost, ncol=n, byrow=TRUE)
    cpo=(delta==2)*pnorm(z, mu., sd.)+(delta==0)*pnorm(z, mu., sd., FALSE)+
        (delta==1)*dnorm(z, mu., sd.)
    cpo = 1/apply(1/cpo, 2, mean)
    res$LPML = sum(log(cpo))

    ##if(length(soffset)==0) {
        res$pred=predict(res, res$x.train, tc=tc,
                         XPtr=FALSE, soffset=0) ##(, draws=TRUE)
        ## res$soffset=(res$pred$s.test.mean/res$s.train.mean)^2
        ## res$soffset=0.5*log(mean(res$soffset[is.finite(res$soffset)]))
        res$soffset=0.5*log(mean(res$pred$s.test.mean^2)/
                            mean(res$s.train.mean^2)) ## for stability
        ##if(!is.finite(res$soffset)) res$soffset=NA
    ##} else res$soffset=soffset
    
    res$drawMuTau=drawMuTau
    
    res$elapsed <- (proc.time()-ptm)['elapsed']
    attr(res$elapsed, 'names')=NULL
    return(res)
}
