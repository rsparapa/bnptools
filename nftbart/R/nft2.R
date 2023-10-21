## Copyright (C) 2022 Rodney A. Sparapani

## This file is part of nftbart.
## nft2.R

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

nft2 = function(## data
               xftrain, xstrain, times, delta=NULL, 
               xftest=matrix(nrow=0, ncol=0),
               xstest=matrix(nrow=0, ncol=0),
               rm.const=TRUE, rm.dupe=TRUE,
               ##edraws2=matrix(nrow=0, ncol=0),
               ##zdraws2=matrix(nrow=0, ncol=0),
               ##impute.bin=NULL, impute.prob=NULL,
               ## multi-threading
               tc=getOption("mc.cores", 1), ##OpenMP thread count
               ##MCMC
               nskip=1000, ndpost=2000, 
               nadapt=1000, adaptevery=100, 
               chvf = NULL, chvs = NULL,
               method="spearman", use="pairwise.complete.obs",
               pbd=c(0.7, 0.7), pb=c(0.5, 0.5),
               stepwpert=c(0.1, 0.1), probchv=c(0.1, 0.1),
               minnumbot=c(5, 5),
               ## BART and HBART prior parameters
               ntree=c(50, 10), numcut=100,
               xifcuts=NULL, xiscuts=NULL,
               power=c(2, 2), base=c(0.95, 0.95),
               ## f function
               fmu=NA, k=5, tau=NA, dist='weibull', 
               ## s function
               total.lambda=NA, total.nu=10, mask=NULL,
               ## survival analysis 
               K=100, events=NULL, TSVS=FALSE,
               ## DPM LIO
               drawDPM=1L, 
               alpha=1, alpha.a=1, alpha.b=0.1, alpha.draw=1,
               neal.m=2, constrain=1, 
               m0=0, k0.a=1.5, k0.b=7.5, k0=1, k0.draw=1,
               ##a0=1.5, b0.a=0.5,
               a0=3, b0.a=2, b0.b=1, b0=1, b0.draw=1,
               ## misc
               na.rm=FALSE, probs=c(0.025, 0.975), printevery=100,
               transposed=FALSE, pred=FALSE
               )
{
    n=length(times)
    if(length(delta)==0) delta=rep(1, n)
    ##if(length(sigmav)==0) sigmav=rep(1, n)

    if(!transposed) {
        xftrain=cbind(xftrain)
        if(n!=nrow(xftrain))
            stop('The number of times and rows in xftrain must be the same')
        xstrain=cbind(xstrain)
        if(n!=nrow(xstrain))
            stop('The number of times and rows in xstrain must be the same')

        xftest=cbind(xftest)
        np=nrow(xftest)
        xstest=cbind(xstest)
        if(np!=nrow(xstest))
            stop('The number of rows in xftest and xstest must be the same')

        if(np>0) {
            xf.=bMM(rbind(xftrain, xftest), numcut=numcut, rm.const=rm.const, rm.dupe=rm.dupe,
                   method=method, use=use)
            xftrain=cbind(xf.$X[1:n, ])
            xftest =cbind(xf.$X[n+(1:np), ])
            xs.=bMM(rbind(xstrain, xstest), numcut=numcut, rm.const=rm.const, rm.dupe=rm.dupe,
                   method=method, use=use)
            xstrain=cbind(xs.$X[1:n, ])
            xstest =cbind(xs.$X[n+(1:np), ])
        } else {
            xf.=bMM(xftrain, numcut=numcut, rm.const=rm.const, rm.dupe=rm.dupe,
                   method=method, use=use)
            xftrain=xf.$X
            xs.=bMM(xstrain, numcut=numcut, rm.const=rm.const, rm.dupe=rm.dupe,
                   method=method, use=use)
            xstrain=xs.$X
        }
        xifcuts=xf.$xicuts
        xiscuts=xs.$xicuts
        chvf = xf.$chv
        chvs = xs.$chv
        
        impute=CDimpute(x.train=xftrain, x.test=xftest)
        xftrain=t(impute$x.train)
        xftest=t(impute$x.test)
        impute=CDimpute(x.train=xstrain, x.test=xstest)
        xstrain=t(impute$x.train)
        xstest=t(impute$x.test)
        transposed=TRUE
    } else {
        np=ncol(xftest)
        if(length(chvf)==0) chvf = cor(t(xftrain), method=method, use=use)
        if(length(chvs)==0) chvs = cor(t(xstrain), method=method, use=use)

        if(length(xifcuts)==0) 
            xifcuts=xicuts(xftrain, transposed=transposed, numcut=numcut)
        if(length(xiscuts)==0) 
            xiscuts=xicuts(xstrain, transposed=transposed, numcut=numcut)
    }
    
    pf=nrow(xftrain)
    ps=nrow(xstrain)

    take.logs = (dist!='extreme')
    if(length(events)==0 && K>0) {
        events <- unique(quantile(times, probs=(0:(K-1))/(K-1)))
        attr(events, 'names') <- NULL
    }
    K <- length(events)

    if(length(mask)==1 && 0<mask && mask<1) ndpost=ceiling(ndpost/mask)
    else mask=NULL
        
    ## if(length(edraws2)>0) {
    ##     if(class(edraws2)[1]!='matrix')
    ##         stop('edraws2 must be a matrix')
    ##     if(class(zdraws2)[1]!='matrix')
    ##         stop('zdraws2 must be a matrix')
    ##     ndpost=nrow(edraws2)
    ##     if(ndpost!=nrow(zdraws2))
    ##         stop('number of rows in edraws2 and zdraws2 must be equal')
    ##     if(n!=ncol(edraws2))
    ##         stop(paste0('number of cols in edraws2 must be ', n))
    ##     if(n!=ncol(zdraws2))
    ##         stop(paste0('number of cols in zdraws2 must be ', n))
    ##     ndpost=ndpost-nskip
    ## }
    
    burn=nskip
    nd=ndpost
    nc=numcut

    df = data.frame(times, delta)
    aft = survreg(Surv(time=times, event=delta)~1,
                  data=subset(df, delta<2), dist=dist)

    if(is.na(fmu)) {
        fmu=coef(aft)[1]
        attr(fmu, 'names')=NULL
    }
    
    if(take.logs) {
        y=log(times)-fmu
        if(K>0) events=log(events)
    } else {
        y=times-fmu
    }

    if(is.na(total.lambda)) total.lambda = (aft$scale^2)

    ## tau
    if(is.na(tau)) {
        ##if(is.na(sigmaf)) {
            tau=(max(y)-min(y))/(2*k*sqrt(ntree[1]))
        ## } else {
        ##     tau = sigmaf/sqrt(ntree[1])
        ## }
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

    ## check = length(impute.bin)
    ## if(!(check %in% 0:1))
    ##     stop("The number of imputed binomial columns must be 0 or 1")
    ## if(check==1) {
    ##     impute.flag=TRUE
    ##     if(length(impute.prob)==0) {
    ##         impute.prob=x[impute.bin, ] ## transposed
    ##         impute.prob[is.na(impute.prob)]=
    ##             mean(x[impute.bin, ]==1, na.rm = TRUE) ## transposed
    ##     }
    ##     impute.bin=as.integer(impute.bin-1) ## convert from R index to C/C++
    ## } else {
    ##     impute.flag=FALSE
    ##     impute.bin=-1
    ##     impute.prob=0
    ## }

    ptm <- proc.time()

    res=.Call("cnft2",
              xftrain, xstrain,
              y,
              as.integer(delta),
              xftest, xstest,
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
              ##sigmav, -1
              chvf,
              chvs,
              pbd,
              pb,
              stepwpert,
              probchv,
              minnumbot,
              printevery,
              xifcuts, xiscuts,
              ##summarystats,
              ##alphao, betao,
              ##mstart, sstart,
              hyper, C, states, phi, prior,
              ##draws,
              drawDPM,
              ##edraws2,
              ##zdraws2,
              ##impute.bin,
              ##impute.prob,
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
    if(length(mask)==1) {
        res$s.train[is.na(res$s.train)] = 0
        s.train.min =apply(res$s.train, 1, quantile, probs=0)
        s.train.max =apply(res$s.train, 1, quantile, probs=1)
        s.train.max.=quantile(s.train.max, probs=mask)
        s.train.mask=(0<s.train.min & s.train.max<=s.train.max.)
        nd=sum(s.train.mask)
        s.train.mask=which(s.train.mask)
        res$f.train=res$f.train[s.train.mask, ]
        res$s.train=res$s.train[s.train.mask, ]
    } else s.train.mask=1:nd 
    ##res$z.train=res$z.train[s.train.mask, ]
    s.train.burn=c(1:burn, s.train.mask)

    summarystats=(ndpost>=2)
    if(summarystats) {
        res$f.varcount[2:ndpost, ]=cbind(res$f.varcount[2:ndpost, ]-res$f.varcount[1:(ndpost-1), ])
        res$f.varcount=cbind(res$f.varcount[s.train.mask, ])
        dimnames(res$f.varcount)[[2]]=names(xifcuts)
        res$f.varcount.mean=apply(res$f.varcount, 2, mean)
        res$f.varprob=res$f.varcount.mean/sum(res$f.varcount.mean)
        res$s.varcount[2:ndpost, ]=cbind(res$s.varcount[2:ndpost, ]-res$s.varcount[1:(ndpost-1), ])
        res$s.varcount=cbind(res$s.varcount[s.train.mask, ])
        dimnames(res$s.varcount)[[2]]=names(xiscuts)
        res$s.varcount.mean=apply(res$s.varcount, 2, mean)
        res$s.varprob=res$s.varcount.mean/sum(res$s.varcount.mean)
    }

    if(TSVS) return(res)
   
    if(drawDPM>0) {
        res$dpalpha=res$dpalpha[s.train.mask]
        res$dpn=res$dpn[s.train.burn]
        res$dpmu=res$dpmu[s.train.mask, ]
        res$dpsd=res$dpsd[s.train.mask, ]
        res$dpmu.=res$dpmu.[s.train.mask, ]
        res$dpsd.=res$dpsd.[s.train.mask, ]
        res$dpC=res$dpC[s.train.mask, ]
        if(burn>0) res$dpn.=res$dpn[-(1:burn)]
        else res$dpn.=res$dpn
    }

    res$f.train=res$f.train+fmu
    res$f.train.mean=apply(res$f.train, 2, mean)
    res$f.train.min=apply(res$f.train, 1, quantile, probs=0)
    res$f.train.max=apply(res$f.train, 1, quantile, probs=1)
    ## if(take.logs) {
    ##     z.trunc=apply(res$z.train, 2, pmin, max(log(times)))
    ##     res$z.train.mean=log(apply(exp(z.trunc), 2, mean))
    ## } else {
    ##     z.trunc=apply(res$z.train, 2, pmin, max(times))
    ##     res$z.train.mean=apply(z.trunc, 2, mean)
    ## }
    res$s.train.mean=apply(res$s.train, 2, mean)
    res$ssd=sqrt(mean(res$s.train.mean^2))
    res$mask=mask
    res$s.train.mask=s.train.mask
    res$s.train.max.=s.train.max.
    res$s.train.min=apply(res$s.train, 1, quantile, probs=0)
    res$s.train.max=apply(res$s.train, 1, quantile, probs=1)
        
    if(np>0) {
        res$f.test=res$f.test[s.train.mask, ]+fmu
        res$f.test.mean=apply(res$f.test, 2, mean)
        res$s.test=res$s.test[s.train.mask, ]
        res$s.test.mean=apply(res$s.test, 2, mean)
        res$xftest=t(xftest)
        res$xstest=t(xstest)
    }

    res$times=times
    res$delta=delta
    res$xftrain=t(xftrain)
    res$xstrain=t(xstrain)
    res$ntree=ntree[1]
    res$ntree[2]=ntree[2]
    res$ndpost=ndpost
    res$xifcuts=xifcuts
    res$xiscuts=xiscuts
    ##res$fmu = fmu
    res$NFT=list(total.lambda = total.lambda, total.nu = total.nu,
                 fmu=fmu, tau=tau, k=k) ## , sigmaf=sigmaf)

    if(drawDPM>0) {
        ## res$prior = prior
        ## res$hyper = hyper
        ## res$C = as.integer(C+1)
        ## res$states = states
        ## res$phi = phi
        res$LIO=list(prior=prior, hyper=hyper, C = as.integer(C+1),
                     states=states, phi=phi)

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
    attr(res, 'class') <- 'nft2'
    
    mu. = res$dpmu*res$s.train+res$f.train
    sd. = res$dpsd*res$s.train
    ## cpo = 1/apply(1/dnorm(res$z.train, mu., sd.), 2, mean)
    ## res$lpml = sum(log(cpo))
    
    ## z=matrix(log(times), nrow=nd, ncol=n, byrow=TRUE)
    ## delta=matrix(delta, nrow=nd, ncol=n, byrow=TRUE)
    ## tmp=(delta==2)*pnorm(z, mu., sd.)+(delta==0)*pnorm(z, mu., sd., FALSE)+
    ##     (delta==1)*dnorm(z, mu., sd.)
    ## res$CPO = 1/apply(1/tmp, 2, mean)
    ## res$LPML = sum(log(res$CPO))

    ##return(res)
    res$pred=predict(res, res$xftrain, res$xstrain, tc=tc, XPtr=FALSE,
                     soffset=0, probs=probs, na.rm=na.rm) 
    res$soffset=0.5*log(mean(res$pred$s.test.mean^2)/
                        mean(res$s.train.mean^2)) ## for stability
    if(pred) {
        res$pred$soffset=res$soffset
        res$pred$s.test=res$pred$s.test-res$soffset
        res$pred$s.test.mean=res$pred$s.test.mean-res$soffset
        res$pred$s.test.lower=res$pred$s.test.lower-res$soffset
        res$pred$s.test.upper=res$pred$s.test.upper-res$soffset    
    } else {
        res$pred=NULL
    }

    res$drawDPM=drawDPM
    res$aft=aft
    ##res$total.lambda=total.lambda
    res$elapsed <- (proc.time()-ptm)['elapsed']
    attr(res$elapsed, 'names')=NULL
    return(res)
}
