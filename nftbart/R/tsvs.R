## Copyright (C) 2021-2022 Rodney A. Sparapani

## This file is part of nftbart.
## tsvs.R

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

tsvs = function(
               ## data
               x.train, times, delta=NULL,
               rm.const=TRUE, rm.dupe=TRUE,
               ##tsvs args
               K=20, a.=1, b.=0.5, C=0.5,
               rds.file='tsvs.rds', pdf.file='tsvs.pdf',
               ## multi-threading
               tc=getOption("mc.cores", 1), ##OpenMP thread count
               ##MCMC
               nskip=1000, ndpost=2000, 
               nadapt=1000, adaptevery=100, 
               chv = NULL,
               method="spearman", use="pairwise.complete.obs",
               pbd=c(0.7, 0.7), pb=c(0.5, 0.5),
               stepwpert=c(0.1, 0.1), probchv=c(0.1, 0.1),
               minnumbot=c(5, 5),
               ## BART and HBART prior parameters
               ntree=c(10, 2), numcut=100, xicuts=NULL,
               power=c(2, 2), base=c(0.95, 0.95),
               ## f function
               k=5, sigmaf=NA,
               dist='weibull', 
               ## s function
               sigmav=NULL, total.lambda=NA, total.nu=10, mask=0.95,
               ## survival analysis 
               ##K=100, events=NULL, 
               ## DPM LIO
               drawDPM=1L, 
               alpha=1, alpha.a=1, alpha.b=0.1, alpha.draw=1,
               neal.m=2, constrain=1, 
               m0=0, k0.a=1.5, k0.b=7.5, k0=1, k0.draw=1,
               a0=3, b0.a=2, b0.b=1, b0=1, b0.draw=1,
               ## misc
               na.rm=FALSE, probs=c(0.025, 0.975), printevery=100,
               transposed=FALSE
               )
{
    if(K==0) return(K)

    if(transposed) 
        stop('tsvs is run with x.train untransposed, i.e., prior to bMM processing')

    x.=bMM(x.train, numcut=numcut, rm.const=rm.const, rm.dupe=rm.dupe,
           method=method, use=use)
    xicuts=x.$xicuts
    chv   =x.$chv
    dummy =x.$dummy
    impute=CDimpute(x.train=x.$X)
    x.train=impute$x.train

    namesX=dimnames(x.train)[[2]] 
    P=ncol(x.train)
    a=matrix(a., nrow=K, ncol=P)
    b=matrix(b., nrow=K, ncol=P)
    S=matrix(0, nrow=K, ncol=P)
    dimnames(S)[[2]]=namesX
    theta=matrix(nrow=K, ncol=P)
    dimnames(theta)[[2]]=namesX
    gamma=matrix(0, nrow=K, ncol=P)
    dimnames(gamma)[[2]]=namesX
    prob=matrix(nrow=K, ncol=P)
    dimnames(prob)[[2]]=namesX
    varcount=matrix(0, nrow=K, ncol=P)
    dimnames(varcount)[[2]]=namesX
    for(i in 1:K) {
        set.seed(i)
        print(paste('Step:', i))
        if(i>1) for(j in 1:P) {
                a[i, j]=a[i-1, j]
                b[i, j]=b[i-1, j]
        }
        theta[i, ]=rbeta(P, a[i, ], b[i, ])
        S[i, which(theta[i, ]>=C)]=1
        j=sum(S[i, ])
        if(j==0) S[i, sample.int(P, 2)]=1
        else if(j==1) S[i, sample(which(S[i, ]==0), 1)]=1
        
        for(j in 1:P)
            if(S[i, j]==1) S[i, dummy[1, j]:dummy[2, j] ]=1
        
        set.seed(K+i*K)
        pick=(S[i, ]==1)
        x.train.=cbind(x.train[ , pick])
        dimnames(x.train.)[[2]]=namesX[pick]
        chv.=cbind(chv[pick, pick])
        xicuts.=xicuts
        for(j in P:1) if(!pick[j]) xicuts.[[j]]=NULL
        post=nft(x.train=t(x.train.), times=times, delta=delta,
                 ##rm.const=TRUE, rm.dupe=TRUE,
                 ## multi-threading
                 tc=tc, ##OpenMP thread count
                 ##MCMC
                 nskip=nskip, ndpost=ndpost, 
                 nadapt=nadapt, adaptevery=adaptevery, 
                 chv=chv., ##method=method, use=use,
                 pbd=pbd, pb=pb,
                 stepwpert=stepwpert, probchv=probchv,
                 minnumbot=minnumbot,
                 ## BART and HBART prior parameters
                 ntree=ntree, numcut=numcut, xicuts=xicuts.,
                 power=power, base=base,
                 ## f function
                 k=k, sigmaf=sigmaf,
                 dist=dist, 
                 ## s function
                 sigmav=sigmav, total.lambda=total.lambda,
                 total.nu=total.nu, mask=mask,
                 ## survival analysis 
                 K=0, events=NULL, TSVS=TRUE, 
                 ## DPM LIO
                 drawDPM=drawDPM,
                 alpha=alpha, alpha.a=alpha.a,
                 alpha.b=alpha.b, alpha.draw=alpha.draw,
                 neal.m=neal.m, constrain=constrain, 
                 m0=m0, k0.a=k0.a, k0.b=k0.b, k0=k0, k0.draw=k0.draw,
                 a0=a0, b0.a=b0.a, b0.b=b0.b, b0=b0, b0.draw=b0.draw,
                 ## misc
                 na.rm=na.rm, probs=probs, printevery=printevery,
                 transposed=TRUE
                 )

        namesf=dimnames(post$f.varcount)[[2]]
        M=nrow(post$f.varcount)
        for(j in 1:P) {
            if(S[i, j]==1) {
                h=which(namesX[j]==namesf)
                l=post$f.varcount[M, h]+post$s.varcount[M, h]
                if(l>0) {
                    varcount[i, j]=l
                    gamma[i, j]=1
                }
                a[i, j]=a[i, j]+gamma[i, j]
                b[i, j]=b[i, j]+1-gamma[i, j]
            } else {
                b[i, j]=b[i, j]+1
            }
            prob[i, j]=a[i, j]/(a[i, j]+b[i, j])
        }
        if(length(warnings())>0) print(warnings())
        res=list(step=i, prob=prob, S=S, a=a, b=b, gamma=gamma,
                 theta=theta, varcount=varcount, dummy=dummy, impute=impute)
        saveRDS(res, rds.file)

        pdf(file=pdf.file)
        plot(1:i, prob[1:i, 1], type='n', ylim=c(0, 1), xlim=c(0, K),
             xlab='Steps', ylab='Inclusion Probability')
        abline(h=0:1, v=c(0, K))
        abline(h=0.5, col=8, lty=2)
        for(j in 1:P)
            if(prob[i, j]>0.5) {
                if(i==1) points(i, theta[i, j], col=j)
                else lines(1:i, theta[1:i, j], col=j)
                h=sample(1:i, 1)
                text(h, theta[h, j], namesX[j], col=j, pos=1)
            }
        dev.off()
    }

    return(res)
}
