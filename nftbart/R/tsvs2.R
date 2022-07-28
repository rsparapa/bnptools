## Copyright (C) 2021-2022 Rodney A. Sparapani

## This file is part of nftbart.
## tsvs2.R

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

tsvs2 = function(
               ## data
               xftrain, xstrain, times, delta=NULL, 
               ##tsvs args
               K=20, a.=1, b.=0.5, C=0.5,
               rds.file='tsvs2.rds', pdf.file='tsvs2.pdf',
               ## multi-threading
               tc=getOption("mc.cores", 1), ##OpenMP thread count
               ##MCMC
               nskip=1000, ndpost=2000, 
               nadapt=1000, adaptevery=100,
               method="spearman",
               ##chvf = cor(xftrain, method="spearman"),
               ##chvs = cor(xstrain, method="spearman"),
               pbd=c(0.7, 0.7), pb=c(0.5, 0.5),
               stepwpert=c(0.1, 0.1), probchv=c(0.1, 0.1),
               minnumbot=c(5, 5),
               ## BART and HBART prior parameters
               ntree=c(7, 3), numcut=100,
               xifcuts=NULL, xiscuts=NULL,
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
               na.rm=FALSE, probs=c(0.025, 0.975), printevery=100
               )
{
    if(K==0) return(K)

    xf.=bMM(xftrain, numcut=numcut, xicuts=xifcuts)
    xftrain=xf.$X
    dummyf=xf.$dummy

    xs.=bMM(xstrain, numcut=numcut, xicuts=xiscuts)
    xstrain=xs.$X
    dummys=xs.$dummy
    
    Namesf=dimnames(xftrain)[[2]] 
    Pf=ncol(xftrain)
    Af=matrix(0, nrow=K, ncol=Pf)
    Bf=matrix(0, nrow=K, ncol=Pf)
    Sf=matrix(0, nrow=K, ncol=Pf)
    dimnames(Sf)[[2]]=Namesf
    thetaf=matrix(nrow=K, ncol=Pf)
    dimnames(thetaf)[[2]]=Namesf
    gammaf=matrix(0, nrow=K, ncol=Pf)
    dimnames(gammaf)[[2]]=Namesf
    probf=matrix(nrow=K, ncol=Pf)
    dimnames(probf)[[2]]=Namesf
    varcountf=matrix(0, nrow=K, ncol=Pf)
    dimnames(varcountf)[[2]]=Namesf
    Namess=dimnames(xstrain)[[2]] 
    Ps=ncol(xstrain)
    As=matrix(0, nrow=K, ncol=Ps)
    Bs=matrix(0, nrow=K, ncol=Ps)
    Ss=matrix(0, nrow=K, ncol=Ps)
    dimnames(Ss)[[2]]=Namess
    thetas=matrix(nrow=K, ncol=Ps)
    dimnames(thetas)[[2]]=Namess
    gammas=matrix(0, nrow=K, ncol=Ps)
    dimnames(gammas)[[2]]=Namess
    probs=matrix(nrow=K, ncol=Ps)
    dimnames(probs)[[2]]=Namess
    varcounts=matrix(0, nrow=K, ncol=Ps)
    dimnames(varcounts)[[2]]=Namess
    for(i in 1:K) {
        set.seed(i)
        print(paste('Step:', i))
        for(j in 1:Pf) {
            if(i==1) {
                Af[i, j]=a.
                Bf[i, j]=b.
            } else {
                Af[i, j]=Af[i-1, j]
                Bf[i, j]=Bf[i-1, j]
            }
        }
        for(j in 1:Ps) {
            if(i==1) {
                As[i, j]=a.
                Bs[i, j]=b.
            } else {
                As[i, j]=As[i-1, j]
                Bs[i, j]=Bs[i-1, j]
            }
        }
        while(sum(Sf[i, ])==0) { 
            thetaf[i, ]=rbeta(Pf, Af[i, ], Bf[i, ])
            Sf[i, which(thetaf[i, ]>=C)]=1
        }
        for(j in 1:Pf)
            if(Sf[i, j]==1) Sf[i, dummyf[1, j]:dummyf[2, j] ]=1
        while(sum(Ss[i, ])==0) { 
            thetas[i, ]=rbeta(Ps, As[i, ], Bs[i, ])
            Ss[i, which(thetas[i, ]>=C)]=1
        }
        for(j in 1:Ps)
            if(Ss[i, j]==1) Ss[i, dummys[1, j]:dummys[2, j] ]=1
        set.seed(K+i*K)
        xftrain.=cbind(xftrain[ , Sf[i, ]==1])
        dimnames(xftrain.)[[2]]=Namesf[Sf[i, ]==1]
        xstrain.=cbind(xstrain[ , Ss[i, ]==1])
        dimnames(xstrain.)[[2]]=Namess[Ss[i, ]==1]

        post=nft2(xftrain=xftrain., xstrain=xstrain.,
                  times=times, delta=delta, 
                 ## multi-threading
                 tc=tc, ##OpenMP thread count
                 ##MCMC
                 nskip=nskip, ndpost=ndpost, 
                 nadapt=nadapt, adaptevery=adaptevery, 
                 method=method,
                 ##chvf=chvf, chvs=chvs,
                 pbd=pbd, pb=pb,
                 stepwpert=stepwpert, probchv=probchv,
                 minnumbot=minnumbot,
                 ## BART and HBART prior parameters
                 ntree=ntree, numcut=numcut,
                 xifcuts=xifcuts, xiscuts=xiscuts,
                 power=power, base=base,
                 ## f function
                 k=k, sigmaf=sigmaf, dist=dist, 
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
                 na.rm=na.rm, probs=probs, printevery=printevery
                 )

        namesf=dimnames(post$f.varcount)[[2]]
        M=nrow(post$f.varcount)
        namess=dimnames(post$s.varcount)[[2]]
        for(j in 1:Pf) {
            if(Sf[i, j]==1) {
                h=which(Namesf[j]==namesf)
                l=post$f.varcount[M, h]
                ## if(length(l)==0) print(str(post))
                ## else
                    if(l>0) {
                    varcountf[i, j]=l
                    gammaf[i, j]=1
                }
                Af[i, j]=Af[i, j]+gammaf[i, j]
                Bf[i, j]=Bf[i, j]+1-gammaf[i, j]
            } else {
                Bf[i, j]=Bf[i, j]+1
            }
            probf[i, j]=Af[i, j]/(Af[i, j]+Bf[i, j])
        }
        for(j in 1:Ps) {
            if(Ss[i, j]==1) {
                h=which(Namess[j]==namess)
                l=post$s.varcount[M, h]
                ## if(length(l)==0) print(str(post))
                ## else
                    if(l>0) {
                    varcounts[i, j]=l
                    gammas[i, j]=1
                }
                As[i, j]=As[i, j]+gammas[i, j]
                Bs[i, j]=Bs[i, j]+1-gammas[i, j]
            } else {
                Bs[i, j]=Bs[i, j]+1
            }
            probs[i, j]=As[i, j]/(As[i, j]+Bs[i, j])
        }
        if(length(warnings())>0) print(warnings())
        res=list(step=i,
                 probf=probf, Sf=Sf, af=Af, bf=Bf, gammaf=gammaf,
                 thetaf=thetaf, varcountf=varcountf,
                 probs=probs, Ss=Ss, as=As, bs=Bs, gammas=gammas,
                 thetas=thetas, varcounts=varcounts
                 )
        saveRDS(res, rds.file)

        pdf(file=pdf.file)
        par(mfrow=c(2, 1))
        plot(1:i, probf[1:i, 1], type='n', ylim=c(0, 1), xlim=c(0, K),
             xlab='Steps: f(x)', ylab='Inclusion Probability')
        abline(h=0:1, v=c(0, K))
        abline(h=0.5, col=8, lty=3)
        for(j in 1:Pf) 
            if(probf[i, j]>0.5) {
                if(i==1) points(i, thetaf[i, j], col=j)
                else lines(1:i, thetaf[1:i, j], col=j)
                h=sample(1:i, 1)
                text(h, thetaf[h, j], Namesf[j], col=j, pos=1)
            }
        plot(1:i, probs[1:i, 1], type='n', ylim=c(0, 1), xlim=c(0, K),
             xlab='Steps: s(x)', ylab='Inclusion Probability')
        abline(h=0:1, v=c(0, K))
        abline(h=0.5, col=8, lty=3)
        for(j in 1:Ps) 
            if(probs[i, j]>0.5) {
                if(i==1) points(i, thetas[i, j], col=j)
                else lines(1:i, thetas[1:i, j], col=j)
                h=sample(1:i, 1)
                text(h, thetas[h, j], Namess[j], col=j, pos=1)
            }
        par(mfrow=c(1, 1))
        dev.off()
    }

    return(res)
}
